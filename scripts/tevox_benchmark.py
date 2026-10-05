#!/usr/bin/env python3
"""Evaluate TEvoX runs against run-independent semantic truth.

This program is intentionally an evaluator, not a trainer or calibrator.  A
single dataset truth contract can be reused by many runs, methods, conditions
and seeds without copying truth rows or relying on TEvoX-generated IDs.
"""

from __future__ import annotations

import gzip

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


BENCHMARK_SCHEMA_VERSION = "1.0.0"
REPORT_STATUS = "EVALUATION_ONLY_NOT_A_CALIBRATED_MODEL"
UNKNOWN = "."

TECHNICAL_STATES = {"CALLABLE", "GAP", "AMBIGUOUS", "UNCALLABLE"}
BIOLOGICAL_STATES = {"PRESENT", "EMPTY", "STRUCTURAL_ALTERNATIVE", "UNKNOWN"}
ANNOTATION_STATES = {
    "MATCHED", "MISSING", "FAMILY_CONFLICT", "NOT_APPLICABLE", "UNKNOWN"
}
LEGACY_STATES = {
    "PRESENT_ANNOTATED", "PRESENT_UNANNOTATED", "EMPTY_SITE_CONFIRMED",
    "STRUCTURAL_ALTERNATIVE", "FAMILY_OR_BOUNDARY_DISCORDANCE",
    "ASSEMBLY_GAP", "PROJECTION_AMBIGUOUS", "UNCALLABLE",
}
CANDIDATE_LABELS = {"SAME_LOCUS", "DIFFERENT_LOCUS", "UNKNOWN"}


class ContractError(ValueError):
    """Raised when a benchmark input violates its public contract."""


@dataclass(frozen=True)
class Dataset:
    dataset_id: str
    genome_truth: Path
    candidate_truth: Path | None
    locus_truth: Path | None
    instance_truth: Path | None
    candidate_truth_scope: str
    candidate_sampling_design: str
    candidate_inclusion_probability: float | None
    locus_truth_scope: str


@dataclass(frozen=True)
class Run:
    run_id: str
    dataset_id: str
    method_id: str
    condition_id: str
    replicate_id: str
    group_id: str
    prediction_format: str
    prefix: Path


def table_path(path):
    path = Path(path)
    compressed = Path(str(path) + ".gz")
    if path.suffix == ".tsv" and compressed.is_file():
        if path.is_file():
            raise ValueError("ambiguous plain/compressed table: " + str(path))
        return compressed
    return path


def open_table(path):
    path = table_path(path)
    return gzip.open(path, "rt", newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(newline="", encoding="utf-8")


def read_tsv(
    path: Path, required: set[str], label: str, *, allow_empty: bool = False
) -> list[dict[str, str]]:
    try:
        with open_table(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ContractError(f"{label} {path} has no header")
            if any(
                name is None or not name or name != name.strip()
                for name in reader.fieldnames
            ):
                raise ContractError(f"{label} {path} has a malformed header")
            duplicate_columns = [
                name for name, count in Counter(reader.fieldnames).items()
                if count > 1
            ]
            if duplicate_columns:
                raise ContractError(
                    f"{label} {path} has duplicate columns: "
                    + ", ".join(sorted(duplicate_columns))
                )
            missing = required - set(reader.fieldnames)
            if missing:
                raise ContractError(
                    f"{label} {path} is missing columns: "
                    + ", ".join(sorted(missing))
                )
            rows = []
            for number, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise ContractError(
                        f"{label} {path} row {number} has the wrong number of columns"
                    )
                if not any(value.strip() for value in row.values()):
                    raise ContractError(f"{label} {path} row {number} is empty")
                empty_fields = [
                    field for field, value in row.items() if value == ""
                ]
                if empty_fields:
                    raise ContractError(
                        f"{label} {path} row {number} has empty field(s) "
                        f"{', '.join(empty_fields)}; use '.' for missing values"
                    )
                rows.append(row)
    except OSError as error:
        raise ContractError(f"cannot read {label} {path}: {error}") from error
    if not rows and not allow_empty:
        raise ContractError(f"{label} {path} has no data rows")
    return rows


def resolve_path(base: Path, text: str, *, optional: bool = False) -> Path | None:
    value = text.strip()
    if optional and (not value or value == UNKNOWN):
        return None
    if not value or value == UNKNOWN:
        raise ContractError("required path is missing")
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def unique_nonempty(rows: Sequence[dict[str, str]], field: str, label: str) -> None:
    seen: set[str] = set()
    for number, row in enumerate(rows, start=2):
        value = row[field].strip()
        if not value or value == UNKNOWN:
            raise ContractError(f"{label} row {number} has no {field}")
        if value in seen:
            raise ContractError(f"{label} has duplicate {field} {value!r}")
        seen.add(value)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            while block := handle.read(1024 * 1024):
                digest.update(block)
    except OSError as error:
        raise ContractError(f"cannot hash {path}: {error}") from error
    return digest.hexdigest()


def load_datasets(path: Path) -> dict[str, Dataset]:
    required = {
        "benchmark_schema_version", "dataset_id", "genome_truth",
        "candidate_truth",
        "locus_truth", "instance_truth", "candidate_truth_scope",
        "candidate_sampling_design", "candidate_inclusion_probability",
        "locus_truth_scope",
    }
    rows = read_tsv(path, required, "dataset manifest")
    unique_nonempty(rows, "dataset_id", "dataset manifest")
    datasets: dict[str, Dataset] = {}
    for number, row in enumerate(rows, start=2):
        if row["benchmark_schema_version"] != BENCHMARK_SCHEMA_VERSION:
            raise ContractError(
                f"dataset manifest row {number} has unsupported schema "
                f"{row['benchmark_schema_version']!r}"
            )
        locus_scope = row["locus_truth_scope"].strip()
        if locus_scope not in {"EXHAUSTIVE", "PARTIAL", "NOT_EVALUATED"}:
            raise ContractError(
                f"dataset manifest row {number} has invalid locus_truth_scope "
                f"{locus_scope!r}"
            )
        candidate_scope = row["candidate_truth_scope"].strip()
        if candidate_scope not in {"EXHAUSTIVE", "PARTIAL", "NOT_EVALUATED"}:
            raise ContractError(
                f"dataset manifest row {number} has invalid candidate_truth_scope "
                f"{candidate_scope!r}"
            )
        sampling_design = row["candidate_sampling_design"].strip()
        if sampling_design not in {
            "CENSUS", "CASE_CONTROL", "STRATIFIED", "UNKNOWN",
            "NOT_APPLICABLE",
        }:
            raise ContractError(
                f"dataset manifest row {number} has invalid "
                f"candidate_sampling_design {sampling_design!r}"
            )
        inclusion_text = row["candidate_inclusion_probability"].strip()
        inclusion_probability = None
        if inclusion_text != UNKNOWN:
            inclusion_probability = finite_unit(
                inclusion_text, "candidate_inclusion_probability"
            )
            if inclusion_probability == 0.0:
                raise ContractError(
                    "candidate_inclusion_probability must be greater than zero"
                )
        dataset_id = row["dataset_id"].strip()
        genome_truth = resolve_path(path.parent, row["genome_truth"])
        assert genome_truth is not None
        candidate_truth = resolve_path(
            path.parent, row["candidate_truth"], optional=True
        )
        locus_truth = resolve_path(path.parent, row["locus_truth"], optional=True)
        instance_truth = resolve_path(
            path.parent, row["instance_truth"], optional=True
        )
        datasets[dataset_id] = Dataset(
            dataset_id=dataset_id,
            genome_truth=genome_truth,
            candidate_truth=candidate_truth,
            locus_truth=locus_truth,
            instance_truth=instance_truth,
            candidate_truth_scope=candidate_scope,
            candidate_sampling_design=sampling_design,
            candidate_inclusion_probability=inclusion_probability,
            locus_truth_scope=locus_scope,
        )
        if not any((candidate_truth, locus_truth, instance_truth)):
            raise ContractError(
                f"dataset {dataset_id!r} has no truth source"
            )
        if (locus_scope == "NOT_EVALUATED") != (locus_truth is None):
            raise ContractError(
                f"dataset {dataset_id!r} locus truth path/scope are inconsistent"
            )
        if (candidate_scope == "NOT_EVALUATED") != (candidate_truth is None):
            raise ContractError(
                f"dataset {dataset_id!r} candidate truth path/scope are inconsistent"
            )
        if candidate_scope == "EXHAUSTIVE" and (
            sampling_design != "CENSUS" or inclusion_probability != 1.0
        ):
            raise ContractError(
                f"dataset {dataset_id!r} EXHAUSTIVE candidate truth requires "
                "CENSUS sampling with inclusion probability 1"
            )
        if candidate_scope == "NOT_EVALUATED" and (
            sampling_design != "NOT_APPLICABLE" or inclusion_probability is not None
        ):
            raise ContractError(
                f"dataset {dataset_id!r} without candidate truth requires "
                "NOT_APPLICABLE sampling and '.' inclusion probability"
            )
    return datasets


def valid_sha256(value: object, label: str) -> str:
    if (
        not isinstance(value, str) or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ContractError(f"{label} must be a lowercase SHA-256 digest")
    return value


def load_genome_truth(
    path: Path, dataset_id: str
) -> dict[str, dict[str, str]]:
    required = {
        "benchmark_schema_version", "dataset_id", "genome_id",
        "assembly_sha256", "annotation_sha256", "taxon",
    }
    rows = read_tsv(path, required, "genome truth")
    validate_truth_identity(rows, dataset_id, "genome truth")
    unique_nonempty(rows, "genome_id", "genome truth")
    result: dict[str, dict[str, str]] = {}
    for number, row in enumerate(rows, start=2):
        valid_sha256(row["assembly_sha256"], f"genome truth row {number} assembly")
        valid_sha256(row["annotation_sha256"], f"genome truth row {number} annotation")
        if not row["taxon"].strip() or row["taxon"] == UNKNOWN:
            raise ContractError(f"genome truth row {number} has no taxon")
        result[row["genome_id"]] = row
    return result


def validate_run_input_contract(
    run_json: dict[str, object], genome_truth: dict[str, dict[str, str]],
    run_id: str,
) -> None:
    input_files = run_json.get("input_files")
    if not isinstance(input_files, list) or not input_files:
        raise ContractError(f"run {run_id!r} has no frozen input_files inventory")
    frozen_inputs: dict[str, str] = {}
    for index, item in enumerate(input_files):
        if not isinstance(item, dict):
            raise ContractError(
                f"run {run_id!r} input_files[{index}] is not an object"
            )
        role = item.get("role")
        path = item.get("path")
        if not isinstance(role, str) or not role:
            raise ContractError(f"run {run_id!r} input_files[{index}] has no role")
        if not isinstance(path, str) or not path:
            raise ContractError(f"run {run_id!r} input_files[{index}] has no path")
        digest = valid_sha256(
            item.get("sha256"), f"run {run_id!r} input_files[{index}]"
        )
        previous = frozen_inputs.setdefault(path, digest)
        if previous != digest:
            raise ContractError(
                f"run {run_id!r} records conflicting hashes for input {path!r}"
            )
    counts = run_json.get("counts")
    if not isinstance(counts, dict) or counts.get("input_files") != len(input_files):
        raise ContractError(f"run {run_id!r} input_files count is inconsistent")
    genomes = run_json.get("genomes")
    if not isinstance(genomes, list):
        raise ContractError(f"run {run_id!r} has no genomes array")
    observed: dict[str, dict[str, object]] = {}
    for index, item in enumerate(genomes):
        if not isinstance(item, dict):
            raise ContractError(f"run {run_id!r} genomes[{index}] is not an object")
        genome_id = item.get("genome_id")
        if not isinstance(genome_id, str) or not genome_id or genome_id in observed:
            raise ContractError(f"run {run_id!r} has missing/duplicate genome ID")
        observed[genome_id] = item
    if set(observed) != set(genome_truth):
        raise ContractError(
            f"run {run_id!r} genome IDs disagree with frozen dataset genome truth"
        )
    for genome_id, frozen in genome_truth.items():
        item = observed[genome_id]
        assembly = valid_sha256(
            item.get("fasta_sha256"),
            f"run {run_id!r} genome {genome_id!r} FASTA",
        )
        annotation = valid_sha256(
            item.get("te_annotation_sha256"),
            f"run {run_id!r} genome {genome_id!r} annotation",
        )
        if assembly != frozen["assembly_sha256"]:
            raise ContractError(
                f"run {run_id!r} genome {genome_id!r} assembly hash disagrees "
                "with dataset truth"
            )
        if annotation != frozen["annotation_sha256"]:
            raise ContractError(
                f"run {run_id!r} genome {genome_id!r} annotation hash disagrees "
                "with dataset truth"
            )
        if frozen_inputs.get(item.get("fasta")) != assembly:
            raise ContractError(
                f"run {run_id!r} genome {genome_id!r} FASTA is not bound "
                "to input_files"
            )
        if frozen_inputs.get(item.get("te_annotation")) != annotation:
            raise ContractError(
                f"run {run_id!r} genome {genome_id!r} annotation is not "
                "bound to input_files"
            )
    for section in ("alignment_evidence", "synteny_evidence"):
        values = run_json.get(section)
        if not isinstance(values, list):
            raise ContractError(f"run {run_id!r} has no {section} array")
        for index, item in enumerate(values):
            if not isinstance(item, dict):
                raise ContractError(
                    f"run {run_id!r} {section}[{index}] is not an object"
                )
            digest = valid_sha256(
                item.get("path_sha256"),
                f"run {run_id!r} {section}[{index}] input",
            )
            if frozen_inputs.get(item.get("path")) != digest:
                raise ContractError(
                    f"run {run_id!r} {section}[{index}] is not bound to "
                    "input_files"
                )
    synteny_inputs = run_json.get("synteny_inputs")
    if not isinstance(synteny_inputs, list):
        raise ContractError(f"run {run_id!r} has no synteny_inputs array")
    for index, item in enumerate(synteny_inputs):
        if not isinstance(item, dict):
            raise ContractError(
                f"run {run_id!r} synteny_inputs[{index}] is not an object"
            )
        if item.get("provider") != "MCScanX":
            raise ContractError(
                f"run {run_id!r} synteny_inputs[{index}] has unsupported provider"
            )
        for path_field, hash_field in (
            ("collinearity_path", "collinearity_sha256"),
            ("gene_table_path", "gene_table_sha256"),
        ):
            path_value = item.get(path_field)
            if not isinstance(path_value, str) or not path_value:
                raise ContractError(
                    f"run {run_id!r} synteny_inputs[{index}] has no {path_field}"
                )
            digest = valid_sha256(
                item.get(hash_field),
                f"run {run_id!r} synteny_inputs[{index}] {path_field}",
            )
            if frozen_inputs.get(path_value) != digest:
                raise ContractError(
                    f"run {run_id!r} synteny_inputs[{index}] {path_field} "
                    "is not bound to input_files"
                )


def validate_dataset_truth_genomes(
    dataset: Dataset, genome_truth: dict[str, dict[str, str]]
) -> None:
    if dataset.candidate_truth is not None:
        candidates = load_candidate_truth(
            dataset.candidate_truth, dataset.dataset_id
        )
        for row in candidates.values():
            for genome_field, taxon_field in (
                ("source_genome_id", "taxon_a"),
                ("target_genome_id", "taxon_b"),
            ):
                genome_id = row[genome_field]
                if genome_id not in genome_truth:
                    raise ContractError(
                        f"candidate truth references unknown genome {genome_id!r}"
                    )
                if row[taxon_field] != genome_truth[genome_id]["taxon"]:
                    raise ContractError(
                        f"candidate truth {taxon_field} disagrees with genome truth"
                    )
    if dataset.locus_truth is not None:
        assignment, _ = load_locus_truth(
            dataset.locus_truth, dataset.dataset_id
        )
        unknown_genomes = {genome for genome, _ in assignment} - set(genome_truth)
        if unknown_genomes:
            raise ContractError(
                f"locus truth references unknown genome(s) {sorted(unknown_genomes)!r}"
            )
    if dataset.instance_truth is not None:
        instances = load_instance_truth(
            dataset.instance_truth, dataset.dataset_id
        )
        unknown_genomes = {
            row["genome_id"] for row in instances
        } - set(genome_truth)
        if unknown_genomes:
            raise ContractError(
                f"instance truth references unknown genome(s) {sorted(unknown_genomes)!r}"
            )


def load_runs(path: Path, datasets: dict[str, Dataset]) -> list[Run]:
    required = {
        "benchmark_schema_version", "run_id", "dataset_id", "method_id",
        "condition_id", "replicate_id", "group_id", "prediction_format",
        "prefix",
    }
    rows = read_tsv(path, required, "run manifest")
    unique_nonempty(rows, "run_id", "run manifest")
    runs: list[Run] = []
    for number, row in enumerate(rows, start=2):
        if row["benchmark_schema_version"] != BENCHMARK_SCHEMA_VERSION:
            raise ContractError(
                f"run manifest row {number} has unsupported schema "
                f"{row['benchmark_schema_version']!r}"
            )
        dataset_id = row["dataset_id"].strip()
        if dataset_id not in datasets:
            raise ContractError(
                f"run manifest row {number} references unknown dataset {dataset_id!r}"
            )
        for field in (
            "method_id", "condition_id", "replicate_id", "group_id",
            "prediction_format",
        ):
            if not row[field].strip() or row[field].strip() == UNKNOWN:
                raise ContractError(f"run manifest row {number} has no {field}")
        if row["prediction_format"] != "tevox-1.2":
            raise ContractError(
                f"run manifest row {number} has unsupported prediction_format "
                f"{row['prediction_format']!r}; alpha.2 accepts tevox-1.2"
            )
        prefix = resolve_path(path.parent, row["prefix"])
        assert prefix is not None
        runs.append(Run(
            run_id=row["run_id"].strip(), dataset_id=dataset_id,
            method_id=row["method_id"].strip(),
            condition_id=row["condition_id"].strip(),
            replicate_id=row["replicate_id"].strip(),
            group_id=row["group_id"].strip(),
            prediction_format=row["prediction_format"].strip(), prefix=prefix,
        ))
    return runs


def canonical_pair(
    genome_a: str, te_a: str, genome_b: str, te_b: str
) -> tuple[tuple[str, str], tuple[str, str]]:
    a = (genome_a.strip(), te_a.strip())
    b = (genome_b.strip(), te_b.strip())
    if any(not value or value == UNKNOWN for item in (a, b) for value in item):
        raise ContractError("candidate semantic key contains a missing genome or TE ID")
    if a == b:
        raise ContractError(f"candidate semantic pair repeats the same TE {a!r}")
    return (a, b) if a < b else (b, a)


def finite_unit(text: str, field: str) -> float:
    try:
        value = float(text)
    except ValueError as error:
        raise ContractError(f"invalid {field} {text!r}") from error
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ContractError(f"{field} must be finite and in [0,1], found {text!r}")
    return value


def optional_coordinate(text: str, field: str) -> int | None:
    if not text or text == UNKNOWN:
        return None
    try:
        value = int(text)
    except ValueError as error:
        raise ContractError(f"invalid {field} coordinate {text!r}") from error
    if value < 0:
        raise ContractError(f"{field} coordinate must be non-negative")
    return value


def percentile(values: Sequence[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def auroc(scores: Sequence[tuple[float, int]]) -> float | None:
    positives = sum(label for _, label in scores)
    negatives = len(scores) - positives
    if positives == 0 or negatives == 0:
        return None
    ordered = sorted(scores, key=lambda item: item[0])
    positive_rank_sum = 0.0
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][0] == ordered[start][0]:
            end += 1
        rank = ((start + 1) + end) / 2.0
        positive_rank_sum += rank * sum(label for _, label in ordered[start:end])
        start = end
    return (
        positive_rank_sum - positives * (positives + 1) / 2.0
    ) / (positives * negatives)


def auprc(scores: Sequence[tuple[float, int]]) -> float | None:
    positives = sum(label for _, label in scores)
    if positives == 0:
        return None
    ordered = sorted(scores, key=lambda item: (-item[0], -item[1]))
    true_positive = 0
    false_positive = 0
    previous_recall = 0.0
    area = 0.0
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][0] == ordered[index][0]:
            end += 1
        for _, label in ordered[index:end]:
            true_positive += label
            false_positive += 1 - label
        recall = true_positive / positives
        precision = true_positive / (true_positive + false_positive)
        area += (recall - previous_recall) * precision
        previous_recall = recall
        index = end
    return area


def equal_mass_ece(
    scores: Sequence[tuple[float, int]], requested_bins: int = 10
) -> tuple[float | None, list[dict[str, object]]]:
    if not scores:
        return None, []
    ordered = sorted(scores, key=lambda item: item[0])
    bins = min(requested_bins, len(ordered))
    # Quantile boundaries are score values, not row offsets: all observations
    # with the same score must remain in one bin.  Otherwise ECE depends on the
    # input order of truth labels within a tie and can be arbitrarily inflated.
    boundaries: list[float] = []
    maximum = ordered[-1][0]
    for index in range(1, bins):
        rank = math.ceil(index * len(ordered) / bins) - 1
        boundary = ordered[rank][0]
        if boundary < maximum and (
            not boundaries or boundary > boundaries[-1]
        ):
            boundaries.append(boundary)
    members_by_bin: list[list[tuple[float, int]]] = [
        [] for _ in range(len(boundaries) + 1)
    ]
    for score, label in ordered:
        bin_index = 0
        while bin_index < len(boundaries) and score > boundaries[bin_index]:
            bin_index += 1
        members_by_bin[bin_index].append((score, label))
    rows: list[dict[str, object]] = []
    ece = 0.0
    for index, members in enumerate(members_by_bin, start=1):
        mean_score = sum(score for score, _ in members) / len(members)
        observed = sum(label for _, label in members) / len(members)
        gap = abs(mean_score - observed)
        ece += len(members) / len(ordered) * gap
        rows.append({
            "bin": index,
            "count": len(members),
            "minimum_score": members[0][0],
            "maximum_score": members[-1][0],
            "mean_score": mean_score,
            "observed_fraction": observed,
            "absolute_gap": gap,
        })
    return ece, rows


def binary_metrics(scores: Sequence[tuple[float, int]]) -> dict[str, object]:
    if not scores:
        return {
            "n": 0, "positives": 0, "negatives": 0, "brier_score": None,
            "log_loss": None, "auroc": None, "auprc": None,
            "equal_mass_ece": None,
        }
    epsilon = 1e-15
    brier = sum((score - label) ** 2 for score, label in scores) / len(scores)
    log_loss = -sum(
        label * math.log(min(1.0 - epsilon, max(epsilon, score)))
        + (1 - label) * math.log(
            min(1.0 - epsilon, max(epsilon, 1.0 - score))
        )
        for score, label in scores
    ) / len(scores)
    ece, _ = equal_mass_ece(scores)
    positives = sum(label for _, label in scores)
    return {
        "n": len(scores), "positives": positives,
        "negatives": len(scores) - positives, "brier_score": brier,
        "log_loss": log_loss, "auroc": auroc(scores),
        "auprc": auprc(scores), "equal_mass_ece": ece,
    }


def validate_truth_identity(
    rows: Sequence[dict[str, str]], dataset_id: str, label: str
) -> None:
    for number, row in enumerate(rows, start=2):
        if row["benchmark_schema_version"] != BENCHMARK_SCHEMA_VERSION:
            raise ContractError(
                f"{label} row {number} has unsupported benchmark schema "
                f"{row['benchmark_schema_version']!r}"
            )
        if row["dataset_id"] != dataset_id:
            raise ContractError(
                f"{label} row {number} belongs to dataset {row['dataset_id']!r}, "
                f"expected {dataset_id!r}"
            )


def reject_generated_truth_keys(rows: Sequence[dict[str, str]], label: str) -> None:
    forbidden = {
        "run_id", "candidate_id", "evidence_id", "evidence_group_id",
        "decision_id", "edge_id", "locus_id", "instance_id",
        "solver_component_id", "relation_id",
    }
    present = forbidden & set(rows[0])
    if present:
        raise ContractError(
            f"{label} contains generated prediction key(s): "
            + ", ".join(sorted(present))
        )


def load_candidate_truth(
    path: Path, dataset_id: str
) -> dict[tuple[tuple[str, str], tuple[str, str]], dict[str, str]]:
    required = {
        "benchmark_schema_version", "dataset_id", "truth_record_id",
        "source_genome_id", "source_te_id",
        "target_genome_id", "target_te_id", "label", "ancestral_event_id",
        "source_truth_locus_id", "target_truth_locus_id", "truth_confidence",
        "validation_method",
        "validation_source", "validation_batch_id", "curation_blinded",
        "taxon_a", "taxon_b", "clade_holdout_id",
    }
    rows = read_tsv(path, required, "candidate truth")
    reject_generated_truth_keys(rows, "candidate truth")
    validate_truth_identity(rows, dataset_id, "candidate truth")
    unique_nonempty(rows, "truth_record_id", "candidate truth")
    result: dict[tuple[tuple[str, str], tuple[str, str]], dict[str, str]] = {}
    for number, row in enumerate(rows, start=2):
        for field in (
            "source_genome_id", "source_te_id", "target_genome_id",
            "target_te_id", "source_truth_locus_id",
            "target_truth_locus_id", "validation_method",
            "validation_source", "validation_batch_id",
        ):
            if not row[field].strip() or row[field].strip() == UNKNOWN:
                raise ContractError(
                    f"candidate truth row {number} has no {field}"
                )
        for field in ("taxon_a", "taxon_b", "clade_holdout_id"):
            if not row[field].strip():
                raise ContractError(
                    f"candidate truth row {number} has empty {field}; use '.'"
                )
        if row["label"] not in CANDIDATE_LABELS:
            raise ContractError(
                f"candidate truth row {number} has invalid label {row['label']!r}"
            )
        if row["truth_confidence"] not in {"HIGH", "MEDIUM", "LOW"}:
            raise ContractError(
                f"candidate truth row {number} has invalid truth_confidence"
            )
        if row["curation_blinded"] not in {"true", "false", "unknown"}:
            raise ContractError(
                f"candidate truth row {number} has invalid curation_blinded"
            )
        if row["label"] == "SAME_LOCUS":
            if (
                not row["ancestral_event_id"].strip()
                or row["ancestral_event_id"] == UNKNOWN
            ):
                raise ContractError(
                    f"candidate truth row {number} SAME_LOCUS has no ancestral event"
                )
            if row["source_truth_locus_id"] != row["target_truth_locus_id"]:
                raise ContractError(
                    f"candidate truth row {number} SAME_LOCUS has different locus IDs"
                )
        elif row["label"] == "DIFFERENT_LOCUS":
            if (
                UNKNOWN in {
                    row["source_truth_locus_id"], row["target_truth_locus_id"]
                }
                or row["source_truth_locus_id"] == row["target_truth_locus_id"]
            ):
                raise ContractError(
                    f"candidate truth row {number} DIFFERENT_LOCUS requires "
                    "two known distinct locus IDs"
                )
        key = canonical_pair(
            row["source_genome_id"], row["source_te_id"],
            row["target_genome_id"], row["target_te_id"],
        )
        if key in result:
            raise ContractError(f"candidate truth repeats semantic pair {key!r}")
        result[key] = row
    return result


def predicted_candidate_scores(
    prefix: Path, expected_model_id: str, expected_calibration: str,
) -> tuple[
    dict[tuple[tuple[str, str], tuple[str, str]], float], dict[str, object]
]:
    evidence_path = table_path(Path(f"{prefix}.evidence.tsv"))
    candidates_path = table_path(Path(f"{prefix}.candidates.tsv"))
    features_path = table_path(Path(f"{prefix}.candidate_features.tsv"))
    evidence_rows = read_tsv(
        evidence_path,
        {"schema_version", "evidence_id", "query_genome_id", "source_te_id",
         "target_genome_id"},
        "evidence output", allow_empty=True,
    )
    candidate_rows = read_tsv(
        candidates_path,
        {"schema_version", "candidate_id", "evidence_id", "target_genome_id", "target_te_id"},
        "candidate output", allow_empty=True,
    )
    feature_rows = read_tsv(
        features_path,
        {"schema_version", "candidate_id", "evidence_id", "membership_score",
         "model_id", "calibration_status"},
        "candidate feature output", allow_empty=True,
    )
    evidence: dict[str, dict[str, str]] = {}
    for row in evidence_rows:
        if row["schema_version"] != "1.2.0":
            raise ContractError("evidence output row has a non-1.2.0 schema")
        key = row["evidence_id"]
        if not key or key in evidence:
            raise ContractError(f"missing or duplicate evidence_id {key!r}")
        evidence[key] = row
    candidates: dict[str, dict[str, str]] = {}
    for row in candidate_rows:
        if row["schema_version"] != "1.2.0":
            raise ContractError("candidate output row has a non-1.2.0 schema")
        key = row["candidate_id"]
        if not key or key in candidates:
            raise ContractError(f"missing or duplicate candidate_id {key!r}")
        if row["evidence_id"] not in evidence:
            raise ContractError(f"candidate {key!r} references unknown evidence")
        if row["target_genome_id"] != evidence[row["evidence_id"]]["target_genome_id"]:
            raise ContractError(
                f"candidate {key!r} target genome disagrees with evidence"
            )
        candidates[key] = row
    semantic_scores: dict[tuple[tuple[str, str], tuple[str, str]], float] = {}
    semantic_counts: Counter[tuple[tuple[str, str], tuple[str, str]]] = Counter()
    model_ids: set[str] = set()
    calibration_states: set[str] = set()
    seen_features: set[str] = set()
    for row in feature_rows:
        if row["schema_version"] != "1.2.0":
            raise ContractError("candidate feature row has a non-1.2.0 schema")
        candidate_id = row["candidate_id"]
        if not candidate_id or candidate_id in seen_features:
            raise ContractError(
                f"missing or duplicate feature candidate_id {candidate_id!r}"
            )
        seen_features.add(candidate_id)
        if candidate_id not in candidates:
            raise ContractError(f"feature {candidate_id!r} has no candidate row")
        candidate = candidates[candidate_id]
        if row["evidence_id"] != candidate["evidence_id"]:
            raise ContractError(
                f"feature {candidate_id!r} evidence_id disagrees with candidate"
            )
        observation = evidence[candidate["evidence_id"]]
        key = canonical_pair(
            observation["query_genome_id"], observation["source_te_id"],
            candidate["target_genome_id"], candidate["target_te_id"],
        )
        score = finite_unit(row["membership_score"], "membership_score")
        semantic_scores[key] = max(score, semantic_scores.get(key, -math.inf))
        semantic_counts[key] += 1
        model_ids.add(row["model_id"])
        calibration_states.add(row["calibration_status"])
    if seen_features != set(candidates):
        raise ContractError(
            "candidate and candidate-feature outputs are not one-to-one"
        )
    if feature_rows:
        if model_ids != {expected_model_id}:
            raise ContractError(
                "candidate feature model ID disagrees with run metadata"
            )
        if calibration_states != {expected_calibration}:
            raise ContractError(
                "candidate feature calibration status disagrees with run metadata"
            )
    return semantic_scores, {
        "model_ids": sorted(model_ids),
        "calibration_states": sorted(calibration_states),
        "candidate_rows": len(feature_rows),
        "semantic_candidate_pairs": len(semantic_scores),
        "duplicate_direction_or_evidence_rows": sum(
            count - 1 for count in semantic_counts.values()
        ),
        "input_files": [evidence_path, candidates_path, features_path],
    }


def evaluate_candidates(
    prefix: Path, truth_path: Path, dataset_id: str,
    locus_assignment: dict[tuple[str, str], str] | None,
    truth_scope: str, sampling_design: str,
    inclusion_probability: float | None,
    expected_model_id: str, expected_calibration: str,
) -> dict[str, object]:
    truth = load_candidate_truth(truth_path, dataset_id)
    predicted, metadata = predicted_candidate_scores(
        prefix, expected_model_id, expected_calibration
    )
    predicted_only = set(predicted) - set(truth)
    if truth_scope == "EXHAUSTIVE" and predicted_only:
        examples = sorted(predicted_only)[:3]
        raise ContractError(
            "EXHAUSTIVE candidate truth does not cover predicted semantic "
            f"pair(s), including {examples!r}"
        )
    records: list[tuple[float, int]] = []
    end_to_end_records: list[tuple[float, int]] = []
    generated_positive = 0
    generated_negative = 0
    labelled_positive = 0
    labelled_negative = 0
    unknown = 0
    confidence = Counter()
    for key, row in truth.items():
        label = row["label"]
        confidence[row["truth_confidence"]] += 1
        if locus_assignment is not None:
            source = (row["source_genome_id"], row["source_te_id"])
            target = (row["target_genome_id"], row["target_te_id"])
            if source not in locus_assignment or target not in locus_assignment:
                raise ContractError(
                    f"candidate truth {row['truth_record_id']!r} references a TE "
                    "absent from locus truth"
                )
            source_locus = locus_assignment[source]
            target_locus = locus_assignment[target]
            if row["source_truth_locus_id"] != source_locus:
                raise ContractError(
                    f"candidate truth {row['truth_record_id']!r} source locus "
                    "foreign key disagrees with locus truth"
                )
            if row["target_truth_locus_id"] != target_locus:
                raise ContractError(
                    f"candidate truth {row['truth_record_id']!r} target locus "
                    "foreign key disagrees with locus truth"
                )
            if label != "UNKNOWN" and (
                (source_locus == target_locus) != (label == "SAME_LOCUS")
            ):
                raise ContractError(
                    f"candidate truth {row['truth_record_id']!r} label disagrees "
                    "with locus truth membership"
                )
        if label == "UNKNOWN":
            unknown += 1
            continue
        binary = 1 if label == "SAME_LOCUS" else 0
        labelled_positive += binary
        labelled_negative += 1 - binary
        if key in predicted:
            records.append((predicted[key], binary))
            end_to_end_records.append((predicted[key], binary))
            generated_positive += binary
            generated_negative += 1 - binary
        else:
            end_to_end_records.append((0.0, binary))
    generated_metrics = binary_metrics(records)
    end_to_end_metrics = binary_metrics(end_to_end_records)
    _, bins = equal_mass_ece(records)
    metrics = {
        "truth_scope": truth_scope,
        "sampling_design": sampling_design,
        "uniform_inclusion_probability": inclusion_probability,
        "estimand": "MICRO_AVERAGE_OVER_ASSESSED_SEMANTIC_TE_PAIRS",
        "pair_dependence": (
            "ANCESTRAL_EVENTS_AND_VALIDATION_BATCHES_NOT_INDEPENDENT;NO_CI"
        ),
        "auprc_definition": "TIED_THRESHOLD_STEP_AVERAGE_PRECISION",
        "population_prevalence_or_calibration_claimable": (
            truth_scope == "EXHAUSTIVE" and sampling_design == "CENSUS"
        ),
        "truth_rows": len(truth),
        "labelled_truth_rows": labelled_positive + labelled_negative,
        "unknown_truth_rows": unknown,
        "truth_confidence_counts": dict(sorted(confidence.items())),
        "generated_labelled_rows": len(records),
        "truth_coverage": len(records) / (labelled_positive + labelled_negative)
        if labelled_positive + labelled_negative else None,
        "positive_candidate_generation_recall": generated_positive / labelled_positive
        if labelled_positive else None,
        "negative_truth_generation_fraction": generated_negative / labelled_negative
        if labelled_negative else None,
        "unmatched_truth_rows": labelled_positive + labelled_negative - len(records),
        "unassessed_predicted_semantic_pairs": len(predicted_only),
        "predicted_semantic_pair_truth_coverage": (
            (len(predicted) - len(predicted_only)) / len(predicted)
            if predicted else None
        ),
        "generated_only_score_metrics": generated_metrics,
        "end_to_end_zero_for_ungenerated_metrics": end_to_end_metrics,
        "semantic_view_reduction": "MAX_FIXED_SCORE_ACROSS_REPORTED_ALIGNMENT_VIEWS",
        "calibration_bins": bins,
        **{key: value for key, value in metadata.items() if key != "input_files"},
    }
    return {"metrics": metrics, "input_files": [truth_path, *metadata["input_files"]]}


def load_locus_truth(path: Path, dataset_id: str) -> tuple[
    dict[tuple[str, str], str], dict[str, set[tuple[str, str]]]
]:
    required = {
        "benchmark_schema_version", "dataset_id", "truth_record_id",
        "truth_locus_id", "genome_id", "te_id",
        "ancestral_event_id", "homology_group_id", "validation_batch_id",
    }
    rows = read_tsv(path, required, "locus truth")
    reject_generated_truth_keys(rows, "locus truth")
    validate_truth_identity(rows, dataset_id, "locus truth")
    unique_nonempty(rows, "truth_record_id", "locus truth")
    item_to_locus: dict[tuple[str, str], str] = {}
    loci: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for number, row in enumerate(rows, start=2):
        item = (row["genome_id"].strip(), row["te_id"].strip())
        locus = row["truth_locus_id"].strip()
        if any(not value or value == UNKNOWN for value in (*item, locus)):
            raise ContractError(f"locus truth row {number} has a missing semantic key")
        if not row["ancestral_event_id"].strip() or row["ancestral_event_id"] == UNKNOWN:
            raise ContractError(
                f"locus truth row {number} has no ancestral_event_id"
            )
        if not row["validation_batch_id"].strip() or row["validation_batch_id"] == UNKNOWN:
            raise ContractError(
                f"locus truth row {number} has no validation_batch_id"
            )
        if item in item_to_locus:
            raise ContractError(f"locus truth repeats TE member {item!r}")
        item_to_locus[item] = locus
        loci[locus].add(item)
    return item_to_locus, dict(loci)


def predicted_loci(prefix: Path) -> tuple[
    dict[tuple[str, str], str], dict[str, set[tuple[str, str]]], Path
]:
    path = table_path(Path(f"{prefix}.instances.tsv"))
    rows = read_tsv(
        path, {"schema_version", "locus_id", "genome_id", "copy_count", "member_ids"},
        "instance output", allow_empty=True,
    )
    item_to_locus: dict[tuple[str, str], str] = {}
    loci: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for number, row in enumerate(rows, start=2):
        if row["schema_version"] != "1.2.0":
            raise ContractError(
                f"instance row {number} has non-1.2.0 schema"
            )
        if any(
            not row[field].strip() or row[field] == UNKNOWN
            for field in ("locus_id", "genome_id")
        ):
            raise ContractError(
                f"instance row {number} has a missing locus/genome ID"
            )
        try:
            copy_count = int(row["copy_count"])
        except ValueError as error:
            raise ContractError(f"instance row {number} has invalid copy_count") from error
        if copy_count < 0:
            raise ContractError(f"instance row {number} has negative copy_count")
        members = [] if row["member_ids"] == UNKNOWN else row["member_ids"].split(",")
        if any(not member or member == UNKNOWN for member in members) or len(
            members
        ) != len(set(members)):
            raise ContractError(
                f"instance row {number} has malformed/duplicate member IDs"
            )
        if copy_count != len(members):
            raise ContractError(
                f"instance row {number} copy_count/member_ids disagree"
            )
        for te_id in members:
            item = (row["genome_id"], te_id)
            if item in item_to_locus:
                raise ContractError(f"prediction repeats TE member {item!r}")
            item_to_locus[item] = row["locus_id"]
            loci[row["locus_id"]].add(item)
    return item_to_locus, dict(loci), path


def combinations2(count: int) -> int:
    return count * (count - 1) // 2


def clustering_metrics(
    truth_assignment: dict[tuple[str, str], str],
    predicted_assignment: dict[tuple[str, str], str],
    exhaustive: bool,
) -> dict[str, object]:
    truth_items = set(truth_assignment)
    predicted_items = set(predicted_assignment)
    extra_predicted = predicted_items - truth_items
    if exhaustive and extra_predicted:
        raise ContractError(
            "EXHAUSTIVE locus truth does not contain predicted TE member(s), "
            f"including {sorted(extra_predicted)[:3]!r}"
        )
    # The assessed node universe is always declared by truth.  Prediction-only
    # nodes under PARTIAL truth are unassessed, never invented truth singletons.
    universe = truth_items
    truth_labels: dict[tuple[str, str], str] = dict(truth_assignment)
    predicted_labels: dict[tuple[str, str], str] = {}
    for item in universe:
        predicted_labels[item] = predicted_assignment.get(
            item, f"__MISSING_PREDICTION__:{item[0]}:{item[1]}"
        )
    truth_clusters: dict[str, set[tuple[str, str]]] = defaultdict(set)
    predicted_clusters: dict[str, set[tuple[str, str]]] = defaultdict(set)
    contingency: Counter[tuple[str, str]] = Counter()
    for item in universe:
        truth_clusters[truth_labels[item]].add(item)
        predicted_clusters[predicted_labels[item]].add(item)
        contingency[(truth_labels[item], predicted_labels[item])] += 1
    if not universe:
        raise ContractError("locus evaluation has no TE members")
    b_precision = 0.0
    b_recall = 0.0
    for item in universe:
        if item not in predicted_items:
            # Missing-node detection is part of end-to-end clustering.  A
            # biological singleton that was never emitted must not score 1.
            continue
        overlap = len(
            truth_clusters[truth_labels[item]]
            & predicted_clusters[predicted_labels[item]]
        )
        b_precision += overlap / len(predicted_clusters[predicted_labels[item]])
        b_recall += overlap / len(truth_clusters[truth_labels[item]])
    b_precision /= len(universe)
    b_recall /= len(universe)
    b_f1 = (
        2 * b_precision * b_recall / (b_precision + b_recall)
        if b_precision + b_recall else 0.0
    )

    intersection = truth_items & predicted_items
    intersection_b_precision = 0.0
    intersection_b_recall = 0.0
    if intersection:
        intersection_truth: dict[str, set[tuple[str, str]]] = defaultdict(set)
        intersection_predicted: dict[str, set[tuple[str, str]]] = defaultdict(set)
        for item in intersection:
            intersection_truth[truth_assignment[item]].add(item)
            intersection_predicted[predicted_assignment[item]].add(item)
        for item in intersection:
            overlap = len(
                intersection_truth[truth_assignment[item]]
                & intersection_predicted[predicted_assignment[item]]
            )
            intersection_b_precision += overlap / len(
                intersection_predicted[predicted_assignment[item]]
            )
            intersection_b_recall += overlap / len(
                intersection_truth[truth_assignment[item]]
            )
        intersection_b_precision /= len(intersection)
        intersection_b_recall /= len(intersection)
        intersection_b_f1 = (
            2 * intersection_b_precision * intersection_b_recall
            / (intersection_b_precision + intersection_b_recall)
            if intersection_b_precision + intersection_b_recall else 0.0
        )
    else:
        intersection_b_precision = None
        intersection_b_recall = None
        intersection_b_f1 = None

    true_pairs = sum(combinations2(len(items)) for items in truth_clusters.values())
    predicted_pairs = sum(
        combinations2(len(items)) for items in predicted_clusters.values()
    )
    true_positive_pairs = sum(combinations2(count) for count in contingency.values())
    pair_precision = true_positive_pairs / predicted_pairs if predicted_pairs else (
        1.0 if true_pairs == 0 else 0.0
    )
    pair_recall = true_positive_pairs / true_pairs if true_pairs else 1.0
    pair_f1 = (
        2 * pair_precision * pair_recall / (pair_precision + pair_recall)
        if pair_precision + pair_recall else 0.0
    )

    n_pairs = combinations2(len(universe))
    if n_pairs == 0:
        ari = 1.0
    else:
        expected = true_pairs * predicted_pairs / n_pairs
        maximum = 0.5 * (true_pairs + predicted_pairs)
        denominator = maximum - expected
        ari = (true_positive_pairs - expected) / denominator if denominator else 1.0
    return {
        "truth_scope": "EXHAUSTIVE" if exhaustive else "PARTIAL",
        "evaluated_members": len(universe),
        "truth_members": len(truth_items),
        "predicted_members": len(predicted_items),
        "unassessed_predicted_members": len(extra_predicted),
        "truth_member_prediction_coverage": len(truth_items & predicted_items) / len(truth_items),
        "missing_truth_member_policy": "ZERO_B_CUBED_CONTRIBUTION",
        "b_cubed_precision": b_precision, "b_cubed_recall": b_recall,
        "b_cubed_f1": b_f1, "adjusted_rand_index": ari,
        "intersection_only_b_cubed": {
            "members": len(intersection),
            "precision": intersection_b_precision,
            "recall": intersection_b_recall,
            "f1": intersection_b_f1,
        },
        "pairwise_precision": pair_precision, "pairwise_recall": pair_recall,
        "pairwise_f1": pair_f1, "true_positive_pairs": true_positive_pairs,
        "truth_pairs": true_pairs, "predicted_pairs": predicted_pairs,
        "pairwise_and_ari_detection_warning": (
            "MISSING_TRUTH_MEMBERS_IMPUTED_AS_UNIQUE_SINGLETONS;"
            "INTERPRET_ONLY_WITH_MEMBER_COVERAGE"
        ),
    }


def hungarian_maximum(weights: list[list[int]]) -> list[int]:
    """Return the maximum-weight column for each row of a square matrix."""
    size = len(weights)
    if size == 0:
        return []
    maximum = max(max(row) for row in weights)
    potentials_left = [0] * (size + 1)
    potentials_right = [0] * (size + 1)
    matched_row = [0] * (size + 1)
    predecessor = [0] * (size + 1)
    for row in range(1, size + 1):
        matched_row[0] = row
        minimum = [math.inf] * (size + 1)
        used = [False] * (size + 1)
        column = 0
        while True:
            used[column] = True
            active_row = matched_row[column]
            delta = math.inf
            next_column = 0
            for candidate_column in range(1, size + 1):
                if used[candidate_column]:
                    continue
                cost = maximum - weights[active_row - 1][candidate_column - 1]
                reduced = (
                    cost - potentials_left[active_row]
                    - potentials_right[candidate_column]
                )
                if reduced < minimum[candidate_column]:
                    minimum[candidate_column] = reduced
                    predecessor[candidate_column] = column
                if minimum[candidate_column] < delta:
                    delta = minimum[candidate_column]
                    next_column = candidate_column
            for candidate_column in range(size + 1):
                if used[candidate_column]:
                    potentials_left[matched_row[candidate_column]] += int(delta)
                    potentials_right[candidate_column] -= int(delta)
                else:
                    minimum[candidate_column] -= delta
            column = next_column
            if matched_row[column] == 0:
                break
        while True:
            previous = predecessor[column]
            matched_row[column] = matched_row[previous]
            column = previous
            if column == 0:
                break
    assignment = [-1] * size
    for column in range(1, size + 1):
        if matched_row[column] != 0:
            assignment[matched_row[column] - 1] = column - 1
    return assignment


def maximum_overlap_locus_mapping(
    truth_assignment: dict[tuple[str, str], str],
    predicted_assignment: dict[tuple[str, str], str],
) -> tuple[dict[str, str], dict[str, object]]:
    """Match sparse overlap components exactly, avoiding a global dense matrix."""
    overlap: Counter[tuple[str, str]] = Counter()
    for item, truth_locus in truth_assignment.items():
        predicted_locus = predicted_assignment.get(item)
        if predicted_locus is not None:
            overlap[(truth_locus, predicted_locus)] += 1
    truth_neighbors: dict[str, set[str]] = defaultdict(set)
    predicted_neighbors: dict[str, set[str]] = defaultdict(set)
    for truth_locus, predicted_locus in overlap:
        truth_neighbors[truth_locus].add(predicted_locus)
        predicted_neighbors[predicted_locus].add(truth_locus)
    raw_mapping: dict[str, str] = {}
    visited_truth: set[str] = set()
    visited_predicted: set[str] = set()
    component_sizes: list[tuple[int, int]] = []
    for seed in sorted(truth_neighbors):
        if seed in visited_truth:
            continue
        truth_component: set[str] = set()
        predicted_component: set[str] = set()
        pending: list[tuple[str, str]] = [("truth", seed)]
        while pending:
            side, locus = pending.pop()
            if side == "truth":
                if locus in visited_truth:
                    continue
                visited_truth.add(locus)
                truth_component.add(locus)
                pending.extend(
                    ("predicted", neighbor)
                    for neighbor in truth_neighbors[locus]
                    if neighbor not in visited_predicted
                )
            else:
                if locus in visited_predicted:
                    continue
                visited_predicted.add(locus)
                predicted_component.add(locus)
                pending.extend(
                    ("truth", neighbor)
                    for neighbor in predicted_neighbors[locus]
                    if neighbor not in visited_truth
                )
        truth_loci = sorted(truth_component)
        predicted_loci = sorted(predicted_component)
        if len(truth_loci) == 1:
            weights = [[
                overlap[(truth_loci[0], predicted_locus)]
                for predicted_locus in predicted_loci
            ]]
            assignment = [max(
                range(len(predicted_loci)), key=lambda index: weights[0][index]
            )]
        elif len(predicted_loci) == 1:
            weights = [
                [overlap[(truth_locus, predicted_loci[0])]]
                for truth_locus in truth_loci
            ]
            winner = max(
                range(len(truth_loci)), key=lambda index: weights[index][0]
            )
            assignment = [-1] * len(truth_loci)
            assignment[winner] = 0
        else:
            size = max(len(truth_loci), len(predicted_loci))
            if size > 256:
                raise ContractError(
                    "locus overlap component exceeds the alpha.2 exact mapping "
                    f"cap (truth={len(truth_loci)}, predicted={len(predicted_loci)}, "
                    "maximum dimension=256)"
                )
            weights = [[0 for _ in range(size)] for _ in range(size)]
            for truth_index, truth_locus in enumerate(truth_loci):
                for predicted_index, predicted_locus in enumerate(predicted_loci):
                    weights[truth_index][predicted_index] = overlap[
                        (truth_locus, predicted_locus)
                    ]
            assignment = hungarian_maximum(weights)
        for truth_index, predicted_index in enumerate(assignment[:len(truth_loci)]):
            if (
                0 <= predicted_index < len(predicted_loci)
                and weights[truth_index][predicted_index] > 0
            ):
                raw_mapping[truth_loci[truth_index]] = predicted_loci[predicted_index]
        component_sizes.append((len(truth_loci), len(predicted_loci)))
    all_truth_loci = set(truth_assignment.values())
    all_predicted_loci = set(predicted_assignment.values())
    truth_sizes = Counter(truth_assignment.values())
    predicted_sizes = Counter(predicted_assignment.values())
    mapping: dict[str, str] = {}
    mapping_rows: list[dict[str, object]] = []
    for truth_locus in sorted(all_truth_loci):
        predicted_locus = raw_mapping.get(truth_locus)
        if predicted_locus is None:
            mapping_rows.append({
                "truth_locus_id": truth_locus,
                "predicted_locus_id": None,
                "overlap": 0, "truth_locus_size": truth_sizes[truth_locus],
                "predicted_locus_size": None, "completeness": 0.0,
                "purity": None, "unique_truth_best": False,
                "unique_predicted_best": False, "status": "NO_OVERLAP",
            })
            continue
        weight = overlap[(truth_locus, predicted_locus)]
        truth_best = max(
            (overlap[(truth_locus, neighbor)]
             for neighbor in truth_neighbors[truth_locus]), default=0
        )
        predicted_best = max(
            (overlap[(neighbor, predicted_locus)]
             for neighbor in predicted_neighbors[predicted_locus]), default=0
        )
        unique_truth_best = sum(
            overlap[(truth_locus, neighbor)] == truth_best
            for neighbor in truth_neighbors[truth_locus]
        ) == 1
        unique_predicted_best = sum(
            overlap[(neighbor, predicted_locus)] == predicted_best
            for neighbor in predicted_neighbors[predicted_locus]
        ) == 1
        completeness = weight / truth_sizes[truth_locus]
        purity = weight / predicted_sizes[predicted_locus]
        accepted = (
            weight == truth_best == predicted_best
            and unique_truth_best and unique_predicted_best
            and completeness > 0.5 and purity > 0.5
        )
        if accepted:
            mapping[truth_locus] = predicted_locus
        mapping_rows.append({
            "truth_locus_id": truth_locus,
            "predicted_locus_id": predicted_locus,
            "overlap": weight, "truth_locus_size": truth_sizes[truth_locus],
            "predicted_locus_size": predicted_sizes[predicted_locus],
            "completeness": completeness, "purity": purity,
            "unique_truth_best": unique_truth_best,
            "unique_predicted_best": unique_predicted_best,
            "status": "ACCEPTED" if accepted else "AMBIGUOUS_OR_WEAK",
        })
    return mapping, {
        "method": (
            "EXACT_HUNGARIAN_PER_SPARSE_COMPONENT_THEN_"
            "RECIPROCAL_UNIQUE_MAJORITY_GATE"
        ),
        "matched_truth_loci": len(mapping),
        "truth_loci": len(all_truth_loci),
        "predicted_loci": len(all_predicted_loci),
        "unmatched_truth_loci": len(all_truth_loci) - len(mapping),
        "unmatched_predicted_loci": len(all_predicted_loci) - len(set(mapping.values())),
        "overlap_components": len(component_sizes),
        "largest_component_truth_loci": max((left for left, _ in component_sizes), default=0),
        "largest_component_predicted_loci": max((right for _, right in component_sizes), default=0),
        "mapping_rows": mapping_rows,
    }


def evaluate_loci(
    prefix: Path, truth_path: Path, scope: str, dataset_id: str
) -> dict[str, object]:
    truth_assignment, _ = load_locus_truth(truth_path, dataset_id)
    predicted_assignment, _, predicted_path = predicted_loci(prefix)
    metrics = clustering_metrics(
        truth_assignment, predicted_assignment, scope == "EXHAUSTIVE"
    )
    return {"metrics": metrics, "input_files": [truth_path, predicted_path]}


def bool_value(text: str, field: str) -> bool:
    if text not in {"true", "false"}:
        raise ContractError(f"invalid {field} boolean {text!r}")
    return text == "true"


def validate_state_contract(
    row: dict[str, str], label: str, number: int
) -> None:
    """Validate known legacy/axis combinations without imputing truth '.'."""
    legacy = row["legacy_state"]
    technical = row["technical_state"]
    biological = row["biological_state"]
    annotation = row["annotation_state"]
    claimable = row["claimable"]
    expected: dict[str, dict[str, str]] = {
        "PRESENT_ANNOTATED": {
            "technical_state": "CALLABLE", "biological_state": "PRESENT",
            "annotation_state": "MATCHED", "claimable": "true",
        },
        "PRESENT_UNANNOTATED": {
            "technical_state": "CALLABLE", "biological_state": "PRESENT",
            "annotation_state": "MISSING", "claimable": "true",
        },
        "EMPTY_SITE_CONFIRMED": {
            "technical_state": "CALLABLE", "biological_state": "EMPTY",
            "annotation_state": "NOT_APPLICABLE", "claimable": "true",
        },
        "STRUCTURAL_ALTERNATIVE": {
            "technical_state": "CALLABLE",
            "biological_state": "STRUCTURAL_ALTERNATIVE",
            "annotation_state": "NOT_APPLICABLE", "claimable": "true",
        },
        "FAMILY_OR_BOUNDARY_DISCORDANCE": {
            "technical_state": "CALLABLE", "biological_state": "PRESENT",
            "annotation_state": "FAMILY_CONFLICT", "claimable": "true",
        },
        "ASSEMBLY_GAP": {"technical_state": "GAP", "claimable": "false"},
        "PROJECTION_AMBIGUOUS": {
            "technical_state": "AMBIGUOUS", "claimable": "false",
        },
        # UNCALLABLE is also the legacy fallback for technically CALLABLE but
        # biologically insufficient observations, so only its claim gate is fixed.
        "UNCALLABLE": {"claimable": "false"},
    }
    if legacy != UNKNOWN:
        values = {
            "technical_state": technical, "biological_state": biological,
            "annotation_state": annotation, "claimable": claimable,
        }
        for field, required in expected[legacy].items():
            if values[field] != UNKNOWN and values[field] != required:
                raise ContractError(
                    f"{label} row {number} {field}={values[field]!r} is "
                    f"incompatible with legacy_state={legacy!r}"
                )
    if claimable == "true" and technical not in {"CALLABLE", UNKNOWN}:
        raise ContractError(
            f"{label} row {number} claimable=true requires CALLABLE technical state"
        )


def validate_predicted_coordinates(
    row: dict[str, str], copy_count: int, number: int
) -> tuple[int | None, int | None]:
    if copy_count > 1:
        raise ContractError(
            "alpha.2 state/breakpoint benchmark supports copy_count <= 1; "
            f"instance row {number} has copy_count={copy_count}"
        )
    fields = (row["contig"], row["start"], row["end"])
    if any(value == "" for value in fields):
        raise ContractError(
            f"predicted instance row {number} has an empty coordinate field; use '.'"
        )
    present = tuple(value != UNKNOWN for value in fields)
    if any(present) and not all(present):
        raise ContractError(
            f"predicted instance row {number} must provide contig/start/end together"
        )
    if not any(present):
        return None, None
    if any("," in value for value in fields):
        raise ContractError(
            f"predicted instance row {number} contains unsupported multi-copy coordinates"
        )
    start = optional_coordinate(row["start"], "predicted start")
    end = optional_coordinate(row["end"], "predicted end")
    assert start is not None and end is not None
    if end < start:
        raise ContractError(
            f"predicted instance row {number} has end before start"
        )
    return start, end


def confusion_metrics(pairs: Sequence[tuple[str, str]]) -> dict[str, object]:
    labels = sorted({value for pair in pairs for value in pair})
    matrix: dict[str, dict[str, int]] = {
        truth: {prediction: 0 for prediction in labels} for truth in labels
    }
    for truth, prediction in pairs:
        matrix[truth][prediction] += 1
    per_class: dict[str, dict[str, object]] = {}
    for label in labels:
        true_positive = matrix[label][label]
        false_positive = sum(matrix[truth][label] for truth in labels if truth != label)
        false_negative = sum(matrix[label][prediction] for prediction in labels if prediction != label)
        support = sum(matrix[label].values())
        precision = true_positive / (true_positive + false_positive) if true_positive + false_positive else None
        recall = true_positive / support if support else None
        # For every class present in truth, F1 is defined by counts and is zero
        # when the class is completely missed.  Prediction-only sentinel labels
        # (for example NO_PREDICTION) have no truth support and stay undefined.
        f1_denominator = 2 * true_positive + false_positive + false_negative
        f1 = 2 * true_positive / f1_denominator if support else None
        per_class[label] = {
            "support": support, "precision": precision,
            "recall": recall, "f1": f1,
        }
    f1_values = [
        row["f1"] for row in per_class.values() if row["support"] > 0
    ]
    return {
        "n": len(pairs),
        "accuracy": sum(truth == prediction for truth, prediction in pairs) / len(pairs)
        if pairs else None,
        "macro_f1": sum(f1_values) / len(f1_values) if f1_values else None,
        "labels": labels, "confusion": matrix, "per_class": per_class,
    }


def load_instance_truth(path: Path, dataset_id: str) -> list[dict[str, str]]:
    required = {
        "benchmark_schema_version", "dataset_id", "truth_record_id",
        "truth_locus_id", "genome_id", "technical_state", "biological_state",
        "annotation_state", "legacy_state", "claimable", "contig", "start",
        "end", "truth_confidence", "validation_method", "validation_source",
        "validation_batch_id",
    }
    rows = read_tsv(path, required, "instance truth")
    reject_generated_truth_keys(rows, "instance truth")
    validate_truth_identity(rows, dataset_id, "instance truth")
    unique_nonempty(rows, "truth_record_id", "instance truth")
    semantic_seen: set[tuple[str, str]] = set()
    for number, row in enumerate(rows, start=2):
        key = (row["truth_locus_id"], row["genome_id"])
        if any(not value or value == UNKNOWN for value in key):
            raise ContractError(f"instance truth row {number} has no locus/genome key")
        if key in semantic_seen:
            raise ContractError(f"instance truth repeats locus/genome {key!r}")
        semantic_seen.add(key)
        checks = (
            ("technical_state", TECHNICAL_STATES),
            ("biological_state", BIOLOGICAL_STATES),
            ("annotation_state", ANNOTATION_STATES),
            ("legacy_state", LEGACY_STATES),
        )
        for field, values in checks:
            if row[field] != UNKNOWN and row[field] not in values:
                raise ContractError(
                    f"instance truth row {number} has invalid {field} {row[field]!r}"
                )
        if row["claimable"] != UNKNOWN:
            bool_value(row["claimable"], "claimable")
        if row["truth_confidence"] not in {"HIGH", "MEDIUM", "LOW"}:
            raise ContractError(
                f"instance truth row {number} has invalid truth_confidence"
            )
        for field in ("validation_method", "validation_source", "validation_batch_id"):
            if not row[field].strip() or row[field] == UNKNOWN:
                raise ContractError(
                    f"instance truth row {number} has no {field}"
                )
        start = optional_coordinate(row["start"], "start")
        end = optional_coordinate(row["end"], "end")
        coordinate_parts = (row["contig"] != UNKNOWN, start is not None, end is not None)
        if any(coordinate_parts) and not all(coordinate_parts):
            raise ContractError(
                f"instance truth row {number} must provide contig/start/end together"
            )
        if start is not None and end is not None and end < start:
            raise ContractError(f"instance truth row {number} has end before start")
        validate_state_contract(row, "instance truth", number)
    return rows


def evaluate_instances(
    prefix: Path, truth_path: Path, locus_truth_path: Path, dataset_id: str
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]], list[Path]]:
    truth = load_instance_truth(truth_path, dataset_id)
    truth_assignment, _ = load_locus_truth(locus_truth_path, dataset_id)
    known_truth_loci = set(truth_assignment.values())
    for row in truth:
        if row["truth_locus_id"] not in known_truth_loci:
            raise ContractError(
                f"instance truth {row['truth_record_id']!r} references unknown "
                f"truth locus {row['truth_locus_id']!r}"
            )
    predicted_assignment, _, locus_path = predicted_loci(prefix)
    locus_mapping, mapping_summary = maximum_overlap_locus_mapping(
        truth_assignment, predicted_assignment
    )
    instance_path = table_path(Path(f"{prefix}.instances.tsv"))
    predicted_rows = read_tsv(
        instance_path,
        {"schema_version", "locus_id", "genome_id", "technical_state", "biological_state",
         "annotation_state", "legacy_state", "claimable", "quality", "contig",
         "start", "end", "copy_count", "member_ids"},
        "instance output", allow_empty=True,
    )
    predicted: dict[tuple[str, str], dict[str, str]] = {}
    for number, row in enumerate(predicted_rows, start=2):
        if row["schema_version"] != "1.2.0":
            raise ContractError(
                f"predicted instance row {number} has non-1.2.0 schema"
            )
        key = (row["locus_id"], row["genome_id"])
        if key in predicted:
            raise ContractError(f"prediction repeats instance {key!r}")
        checks = (
            ("technical_state", TECHNICAL_STATES),
            ("biological_state", BIOLOGICAL_STATES),
            ("annotation_state", ANNOTATION_STATES),
            ("legacy_state", LEGACY_STATES),
        )
        for field, values in checks:
            if row[field] not in values:
                raise ContractError(
                    f"predicted instance row {number} has invalid {field}"
                )
        bool_value(row["claimable"], "predicted claimable")
        try:
            copy_count = int(row["copy_count"])
        except ValueError as error:
            raise ContractError(
                f"predicted instance row {number} has invalid copy_count"
            ) from error
        if copy_count < 0:
            raise ContractError(
                f"predicted instance row {number} has negative copy_count"
            )
        members = [] if row["member_ids"] == UNKNOWN else row["member_ids"].split(",")
        if len(members) != copy_count:
            raise ContractError(
                f"predicted instance row {number} copy_count/member_ids disagree"
            )
        predicted_start, predicted_end = validate_predicted_coordinates(
            row, copy_count, number
        )
        if copy_count == 1 and (
            predicted_start is None or predicted_end is None
            or row["technical_state"] != "CALLABLE"
            or row["biological_state"] != "PRESENT"
            or row["annotation_state"] != "MATCHED"
            or row["legacy_state"] != "PRESENT_ANNOTATED"
            or row["claimable"] != "true"
        ):
            raise ContractError(
                f"predicted instance row {number} annotated member/state contract disagrees"
            )
        try:
            quality = float(row["quality"])
        except ValueError as error:
            raise ContractError(
                f"predicted instance row {number} has invalid quality"
            ) from error
        if not math.isfinite(quality) or not 0.0 <= quality <= 100.0:
            raise ContractError(
                f"predicted instance row {number} quality is outside [0,100]"
            )
        validate_state_contract(row, "predicted instance", number)
        predicted[key] = row

    axes: dict[str, list[tuple[str, str]]] = {
        "technical_state": [], "biological_state": [],
        "annotation_state": [], "legacy_state": [], "claimable": [],
    }
    missing_locus_mapping = 0
    missing_instance = 0
    breakpoint_errors: list[dict[str, object]] = []
    quality_records: list[tuple[float, bool, bool, bool]] = []
    empty_claims = 0
    false_empty_claims = 0
    true_empty_claims = 0
    assessed_prediction_keys: set[tuple[str, str]] = set()
    truth_empty = sum(row["biological_state"] == "EMPTY" for row in truth)
    for row in truth:
        truth_start = optional_coordinate(row["start"], "start")
        truth_end = optional_coordinate(row["end"], "end")
        predicted_locus = locus_mapping.get(row["truth_locus_id"])
        if predicted_locus is None:
            missing_locus_mapping += 1
            for field in (
                "technical_state", "biological_state", "annotation_state",
                "legacy_state",
            ):
                if row[field] != UNKNOWN:
                    axes[field].append((row[field], "NO_PREDICTION"))
            if row["claimable"] != UNKNOWN:
                axes["claimable"].append((row["claimable"], "NO_PREDICTION"))
            if row["biological_state"] in {
                "EMPTY", "PRESENT", "STRUCTURAL_ALTERNATIVE"
            }:
                quality_records.append((0.0, False, False, False))
            if truth_start is not None and truth_end is not None:
                breakpoint_errors.append({
                    "truth_record_id": row["truth_record_id"],
                    "truth_locus_id": row["truth_locus_id"],
                    "genome_id": row["genome_id"],
                    "prediction_status": "NO_LOCUS_MAPPING",
                    "same_contig": False,
                    "start_absolute_error": None, "end_absolute_error": None,
                    "mean_boundary_error": None, "max_boundary_error": None,
                })
            continue
        prediction = predicted.get((predicted_locus, row["genome_id"]))
        if prediction is None:
            missing_instance += 1
            for field in (
                "technical_state", "biological_state", "annotation_state",
                "legacy_state",
            ):
                if row[field] != UNKNOWN:
                    axes[field].append((row[field], "NO_PREDICTION"))
            if row["claimable"] != UNKNOWN:
                axes["claimable"].append((row["claimable"], "NO_PREDICTION"))
            if row["biological_state"] in {
                "EMPTY", "PRESENT", "STRUCTURAL_ALTERNATIVE"
            }:
                quality_records.append((0.0, False, False, False))
            if truth_start is not None and truth_end is not None:
                breakpoint_errors.append({
                    "truth_record_id": row["truth_record_id"],
                    "truth_locus_id": row["truth_locus_id"],
                    "genome_id": row["genome_id"],
                    "prediction_status": "NO_INSTANCE",
                    "same_contig": False,
                    "start_absolute_error": None, "end_absolute_error": None,
                    "mean_boundary_error": None, "max_boundary_error": None,
                })
            continue
        biological_assessed = row["biological_state"] in {
            "EMPTY", "PRESENT", "STRUCTURAL_ALTERNATIVE"
        }
        if biological_assessed:
            assessed_prediction_keys.add((predicted_locus, row["genome_id"]))
        for field in ("technical_state", "biological_state", "annotation_state", "legacy_state"):
            if row[field] != UNKNOWN:
                axes[field].append((row[field], prediction[field]))
        if row["claimable"] != UNKNOWN:
            axes["claimable"].append((row["claimable"], prediction["claimable"]))
        try:
            quality = float(prediction["quality"])
        except ValueError as error:
            raise ContractError(f"invalid prediction quality {prediction['quality']!r}") from error
        if not math.isfinite(quality) or not 0.0 <= quality <= 100.0:
            raise ContractError(f"prediction quality outside [0,100]: {quality}")
        if biological_assessed:
            quality_records.append((
                quality,
                prediction["claimable"] == "true",
                prediction["biological_state"] != "UNKNOWN",
                row["biological_state"] == prediction["biological_state"],
            ))
        if (
            prediction["legacy_state"] == "EMPTY_SITE_CONFIRMED"
            and prediction["biological_state"] == "EMPTY"
            and prediction["claimable"] == "true"
            and biological_assessed
        ):
            empty_claims += 1
            if row["biological_state"] == "EMPTY":
                true_empty_claims += 1
            else:
                false_empty_claims += 1
        if truth_start is not None and truth_end is not None:
            predicted_start = optional_coordinate(prediction["start"], "predicted start")
            predicted_end = optional_coordinate(prediction["end"], "predicted end")
            have_coordinates = predicted_start is not None and predicted_end is not None
            same_contig = have_coordinates and prediction["contig"] == row["contig"]
            start_error = abs(predicted_start - truth_start) \
                if same_contig and predicted_start is not None else None
            end_error = abs(predicted_end - truth_end) \
                if same_contig and predicted_end is not None else None
            breakpoint_errors.append({
                "truth_record_id": row["truth_record_id"],
                "truth_locus_id": row["truth_locus_id"],
                "genome_id": row["genome_id"],
                "prediction_status": (
                    "NO_COORDINATES" if not have_coordinates
                    else "MATCHED" if same_contig else "CONTIG_MISMATCH"
                ),
                "same_contig": same_contig,
                "start_absolute_error": start_error,
                "end_absolute_error": end_error,
                "mean_boundary_error": (start_error + end_error) / 2.0
                if start_error is not None and end_error is not None else None,
                "max_boundary_error": max(start_error, end_error)
                if start_error is not None and end_error is not None else None,
            })
    state_metrics = {axis: confusion_metrics(records) for axis, records in axes.items()}
    unassessed_empty_claims = sum(
        key not in assessed_prediction_keys
        and row["legacy_state"] == "EMPTY_SITE_CONFIRMED"
        and row["biological_state"] == "EMPTY"
        and row["claimable"] == "true"
        for key, row in predicted.items()
    )
    all_errors = [
        float(row[field]) for row in breakpoint_errors
        if row["start_absolute_error"] is not None
        and row["end_absolute_error"] is not None
        for field in ("start_absolute_error", "end_absolute_error")
    ]
    callable_breakpoints = sum(
        row["start_absolute_error"] is not None and row["end_absolute_error"] is not None
        for row in breakpoint_errors
    )
    breakpoint_summary = {
        "truth_rows_with_coordinates": len(breakpoint_errors),
        "callable_same_contig_rows": callable_breakpoints,
        "coordinate_coverage": callable_breakpoints / len(breakpoint_errors)
        if breakpoint_errors else None,
        "mean_absolute_error": sum(all_errors) / len(all_errors) if all_errors else None,
        "median_absolute_error": percentile(all_errors, 0.5),
        "rmse": math.sqrt(sum(value * value for value in all_errors) / len(all_errors))
        if all_errors else None,
        "p90_absolute_error": percentile(all_errors, 0.90),
        "p95_absolute_error": percentile(all_errors, 0.95),
        "within_threshold_recall": {
            str(threshold): sum(
                row["max_boundary_error"] is not None
                and float(row["max_boundary_error"]) <= threshold
                for row in breakpoint_errors
            ) / len(breakpoint_errors)
            if breakpoint_errors else None
            for threshold in (0, 5, 10, 25, 50, 100)
        },
    }
    callability_rows: list[dict[str, object]] = []
    for threshold in range(0, 101, 10):
        retained = [
            record for record in quality_records
            if record[0] >= threshold and record[1] and record[2]
        ]
        accuracy = sum(record[3] for record in retained) / len(retained) \
            if retained else None
        callability_rows.append({
            "quality_threshold": threshold,
            "retained": len(retained),
            "total": len(quality_records),
            "coverage": len(retained) / len(quality_records) if quality_records else None,
            "biological_state_accuracy": accuracy,
            "selective_risk": 1.0 - accuracy if accuracy is not None else None,
        })
    metrics = {
        "truth_rows": len(truth),
        "evaluated_rows": len(truth) - missing_locus_mapping - missing_instance,
        "missing_locus_mapping_rows": missing_locus_mapping,
        "missing_predicted_instance_rows": missing_instance,
        "instance_coverage": (
            len(truth) - missing_locus_mapping - missing_instance
        ) / len(truth),
        "locus_mapping": mapping_summary,
        "axes": state_metrics,
        "empty_site_claims_evaluated": empty_claims,
        "unassessed_empty_site_claims": unassessed_empty_claims,
        "true_empty_site_claims": true_empty_claims,
        "false_empty_site_claims": false_empty_claims,
        "truth_empty_sites": truth_empty,
        "assessed_empty_site_false_discovery_rate": false_empty_claims / empty_claims
        if empty_claims else None,
        "empty_claim_assessment_coverage": empty_claims / (
            empty_claims + unassessed_empty_claims
        ) if empty_claims + unassessed_empty_claims else None,
        "empty_site_recall": true_empty_claims / truth_empty if truth_empty else None,
        "breakpoints": breakpoint_summary,
    }
    return metrics, callability_rows, breakpoint_errors, [
        truth_path, locus_truth_path, locus_path, instance_path
    ]


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            output: dict[str, object] = {}
            for field in fieldnames:
                value = row.get(field)
                if value is None:
                    output[field] = UNKNOWN
                elif isinstance(value, bool):
                    output[field] = "true" if value else "false"
                else:
                    output[field] = value
            writer.writerow(output)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate multiple TEvoX runs against semantic, run-independent truth"
    )
    parser.add_argument("--datasets", required=True, type=Path)
    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    try:
        datasets = load_datasets(args.datasets.resolve())
        runs = load_runs(args.runs.resolve(), datasets)
        genome_truth_by_dataset: dict[
            str, dict[str, dict[str, str]]
        ] = {}
        for dataset in datasets.values():
            genome_truth = load_genome_truth(
                dataset.genome_truth, dataset.dataset_id
            )
            validate_dataset_truth_genomes(dataset, genome_truth)
            genome_truth_by_dataset[dataset.dataset_id] = genome_truth
        reports: list[dict[str, object]] = []
        callability_output: list[dict[str, object]] = []
        breakpoint_output: list[dict[str, object]] = []
        input_files: set[Path] = {
            args.datasets.resolve(), args.runs.resolve(),
            *(dataset.genome_truth for dataset in datasets.values()),
        }
        for run in runs:
            dataset = datasets[run.dataset_id]
            run_json_path = Path(f"{run.prefix}.run.json")
            try:
                with run_json_path.open(encoding="utf-8") as handle:
                    run_json = json.load(handle)
            except (OSError, json.JSONDecodeError) as error:
                raise ContractError(
                    f"cannot read valid run metadata {run_json_path}: {error}"
                ) from error
            input_files.add(run_json_path)
            report: dict[str, object] = {
                "run_id": run.run_id, "dataset_id": run.dataset_id,
                "method_id": run.method_id, "condition_id": run.condition_id,
                "replicate_id": run.replicate_id, "group_id": run.group_id,
                "prediction_format": run.prediction_format,
                "prefix": str(run.prefix),
                "software_version": run_json.get("version"),
                "evidence_schema_version": run_json.get("schema_version"),
                "model_id": run_json.get("inference_model", {}).get("model_id"),
                "calibration_status": run_json.get("inference_model", {}).get(
                    "calibration_status"
                ),
            }
            if run_json.get("schema_version") != "1.2.0":
                raise ContractError(
                    f"run {run.run_id!r} uses evidence schema "
                    f"{run_json.get('schema_version')!r}; expected '1.2.0'"
                )
            validate_run_input_contract(
                run_json, genome_truth_by_dataset[run.dataset_id], run.run_id
            )
            model = run_json.get("inference_model")
            if not isinstance(model, dict):
                raise ContractError(
                    f"run {run.run_id!r} has no inference_model metadata"
                )
            model_id = model.get("model_id")
            calibration_status = model.get("calibration_status")
            if not isinstance(model_id, str) or not model_id:
                raise ContractError(f"run {run.run_id!r} has no model_id")
            if calibration_status != "UNCALIBRATED":
                raise ContractError(
                    f"run {run.run_id!r} must remain explicitly UNCALIBRATED"
                )
            locus_assignment = None
            if dataset.locus_truth is not None:
                locus_assignment, _ = load_locus_truth(
                    dataset.locus_truth, dataset.dataset_id
                )
            if dataset.candidate_truth is not None:
                if run_json.get("config", {}).get("max_candidates") != 0:
                    raise ContractError(
                        f"run {run.run_id!r} has max_candidates != 0; semantic "
                        "candidate evaluation requires an untruncated report"
                    )
                result = evaluate_candidates(
                    run.prefix, dataset.candidate_truth, dataset.dataset_id,
                    locus_assignment, dataset.candidate_truth_scope,
                    dataset.candidate_sampling_design,
                    dataset.candidate_inclusion_probability,
                    model_id, calibration_status,
                )
                report["candidate_membership"] = result["metrics"]
                input_files.update(result["input_files"])
            else:
                report["candidate_membership"] = None
            if dataset.locus_truth is not None and dataset.locus_truth_scope != "NOT_EVALUATED":
                result = evaluate_loci(
                    run.prefix, dataset.locus_truth, dataset.locus_truth_scope,
                    dataset.dataset_id,
                )
                report["locus_clustering"] = result["metrics"]
                input_files.update(result["input_files"])
            else:
                report["locus_clustering"] = None
            if dataset.instance_truth is not None:
                if dataset.locus_truth is None:
                    raise ContractError(
                        f"dataset {dataset.dataset_id!r} has instance truth but "
                        "no locus truth for semantic locus matching"
                    )
                metrics, curves, breakpoints, files = evaluate_instances(
                    run.prefix, dataset.instance_truth, dataset.locus_truth,
                    dataset.dataset_id,
                )
                report["state_and_breakpoint"] = metrics
                input_files.update(files)
                for row in curves:
                    callability_output.append({
                        "run_id": run.run_id, "dataset_id": run.dataset_id,
                        "method_id": run.method_id,
                        "condition_id": run.condition_id,
                        "replicate_id": run.replicate_id,
                        "group_id": run.group_id, **row,
                    })
                for row in breakpoints:
                    breakpoint_output.append({
                        "run_id": run.run_id, "dataset_id": run.dataset_id,
                        "method_id": run.method_id,
                        "condition_id": run.condition_id,
                        "replicate_id": run.replicate_id,
                        "group_id": run.group_id, **row,
                    })
            else:
                report["state_and_breakpoint"] = None
            reports.append(report)
        checksums = {
            str(path): sha256_file(path) for path in sorted(input_files, key=str)
        }
        output_prefix = args.output.resolve()
        metrics_path = Path(f"{output_prefix}.metrics.json")
        callability_path = Path(f"{output_prefix}.callability_accuracy.tsv")
        breakpoint_path = Path(f"{output_prefix}.breakpoints.tsv")
        metrics_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "benchmark_schema_version": BENCHMARK_SCHEMA_VERSION,
            "report_status": REPORT_STATUS,
            "truth_identity": "DATASET_SEMANTIC_KEYS_INDEPENDENT_OF_RUN_ID",
            "score_semantics": "EVALUATED_AS_FIXED_SCORES_NOT_CALIBRATED_PROBABILITIES",
            "run_count": len(reports), "dataset_count": len(datasets),
            "runs": reports, "input_sha256": checksums,
        }
        with metrics_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        write_tsv(
            callability_path,
            ["run_id", "dataset_id", "method_id", "condition_id",
             "replicate_id", "group_id",
             "quality_threshold", "retained", "total", "coverage",
             "biological_state_accuracy", "selective_risk"],
            callability_output,
        )
        write_tsv(
            breakpoint_path,
            ["run_id", "dataset_id", "method_id", "condition_id",
             "replicate_id", "group_id",
             "truth_record_id", "truth_locus_id", "genome_id",
             "prediction_status", "same_contig",
             "start_absolute_error", "end_absolute_error",
             "mean_boundary_error", "max_boundary_error"],
            breakpoint_output,
        )
        print(f"Wrote {metrics_path}, {callability_path} and {breakpoint_path}")
    except ContractError as error:
        raise SystemExit(f"benchmark contract error: {error}") from error


if __name__ == "__main__":
    main()
