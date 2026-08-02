#!/usr/bin/env python3
"""Audit benchmark split manifests for information leakage.

This command validates split assignments only.  It never fits, calibrates, or
selects a model.  Records that share a biological or technical binding entity
must remain in the same partition and folds.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence


BENCHMARK_SCHEMA_VERSION = "1.0.0"
SPLIT_SCHEMA_VERSION = "1.0.0"
PARTITIONS = {"DEVELOPMENT", "CALIBRATION", "EXTERNAL_TEST"}
CANDIDATE_LABELS = {"SAME_LOCUS", "DIFFERENT_LOCUS", "UNKNOWN"}
TRUTH_CONFIDENCE = {"HIGH", "MEDIUM", "LOW"}
CURATION_BLINDED = {"true", "false", "unknown"}
UNKNOWN = "."

TRUTH_COLUMNS = {
    "benchmark_schema_version",
    "dataset_id",
    "truth_record_id",
    "source_genome_id",
    "source_te_id",
    "target_genome_id",
    "target_te_id",
    "label",
    "ancestral_event_id",
    "source_truth_locus_id",
    "target_truth_locus_id",
    "truth_confidence",
    "validation_method",
    "validation_source",
    "validation_batch_id",
    "curation_blinded",
    "taxon_a",
    "taxon_b",
    "clade_holdout_id",
}

LOCUS_TRUTH_COLUMNS = {
    "benchmark_schema_version",
    "dataset_id",
    "truth_record_id",
    "truth_locus_id",
    "genome_id",
    "te_id",
    "ancestral_event_id",
    "homology_group_id",
    "validation_batch_id",
}

SPLIT_COLUMNS = {
    "split_schema_version",
    "dataset_id",
    "truth_record_id",
    "partition",
    "outer_fold",
    "inner_fold",
}

class ContractError(ValueError):
    """Raised when an input table violates the split-audit contract."""


def resolved_path(path: Path, label: str) -> Path:
    try:
        return path.resolve(strict=False)
    except (OSError, RuntimeError) as error:
        raise ContractError(f"cannot resolve {label} path {path}: {error}") from error


def require_output_disjoint_from_inputs(
    output_path: Path, input_paths: Sequence[tuple[str, Path]]
) -> None:
    output_resolved = resolved_path(output_path, "output")
    for label, input_path in input_paths:
        if output_resolved == resolved_path(input_path, label):
            raise ContractError(
                f"resolved output path {output_resolved} would overwrite "
                f"the {label} input"
            )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            while block := handle.read(1024 * 1024):
                digest.update(block)
    except OSError as error:
        raise ContractError(f"cannot read {path}: {error}") from error
    return digest.hexdigest()


def read_tsv(path: Path, required: set[str], label: str) -> list[dict[str, str]]:
    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as error:
        raise ContractError(f"cannot open {label} {path}: {error}") from error

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if fields is None:
            raise ContractError(f"{label} {path} has no header")
        if any(field is None or not field.strip() for field in fields):
            raise ContractError(f"{label} {path} contains an empty header field")
        duplicates = sorted(
            field for field, count in Counter(fields).items() if count > 1
        )
        if duplicates:
            raise ContractError(
                f"{label} {path} has duplicate columns: {', '.join(duplicates)}"
            )
        missing = sorted(required - set(fields))
        if missing:
            raise ContractError(
                f"{label} {path} is missing columns: {', '.join(missing)}"
            )

        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ContractError(
                    f"{label} {path} row {line_number} has the wrong number of columns"
                )
            if not any(value.strip() for value in row.values()):
                raise ContractError(f"{label} {path} row {line_number} is empty")
            cleaned = {field: value.strip() for field, value in row.items()}
            rows.append(cleaned)

    if not rows:
        raise ContractError(f"{label} {path} has no data rows")
    return rows


def require_identifier(
    value: str, field: str, line_number: int, label: str = "candidate truth"
) -> str:
    if not value or value == UNKNOWN:
        raise ContractError(
            f"{label} row {line_number} has no {field}"
        )
    return value


def load_truth(path: Path) -> tuple[str, list[dict[str, str]]]:
    rows = read_tsv(path, TRUTH_COLUMNS, "candidate truth TSV")
    forbidden = {
        "run_id", "candidate_id", "evidence_id", "evidence_group_id",
        "decision_id", "edge_id", "locus_id", "instance_id",
        "solver_component_id", "relation_id",
    }
    generated = forbidden & set(rows[0])
    if generated:
        raise ContractError(
            "candidate truth contains generated prediction key(s): "
            + ", ".join(sorted(generated))
        )
    seen_ids: set[str] = set()
    seen_pairs: set[str] = set()
    datasets: set[str] = set()

    for line_number, row in enumerate(rows, start=2):
        if row["benchmark_schema_version"] != BENCHMARK_SCHEMA_VERSION:
            raise ContractError(
                f"candidate truth row {line_number} has unsupported "
                f"benchmark_schema_version "
                f"{row['benchmark_schema_version']!r}"
            )
        dataset_id = require_identifier(row["dataset_id"], "dataset_id", line_number)
        datasets.add(dataset_id)
        record_id = require_identifier(
            row["truth_record_id"], "truth_record_id", line_number, "locus truth"
        )
        if record_id in seen_ids:
            raise ContractError(
                f"candidate truth has duplicate truth_record_id {record_id!r}"
            )
        seen_ids.add(record_id)

        for field in (
            "source_genome_id",
            "source_te_id",
            "target_genome_id",
            "target_te_id",
        ):
            require_identifier(row[field], field, line_number)
        source_key = (row["source_genome_id"], row["source_te_id"])
        target_key = (row["target_genome_id"], row["target_te_id"])
        if source_key == target_key:
            raise ContractError(
                f"candidate truth row {line_number} repeats the same semantic TE "
                f"{source_key!r}"
            )
        if row["label"] not in CANDIDATE_LABELS:
            raise ContractError(
                f"candidate truth row {line_number} has invalid label "
                f"{row['label']!r}"
            )
        pair = canonical_te_pair(row)
        if pair in seen_pairs:
            raise ContractError(
                f"candidate truth repeats semantic pair {pair}"
            )
        seen_pairs.add(pair)
        if row["label"] == "SAME_LOCUS":
            if row["ancestral_event_id"] == UNKNOWN:
                raise ContractError(
                    f"candidate truth row {line_number} SAME_LOCUS has no "
                    "ancestral_event_id"
                )
            if row["source_truth_locus_id"] != row["target_truth_locus_id"]:
                raise ContractError(
                    f"candidate truth row {line_number} SAME_LOCUS has "
                    "different truth locus IDs"
                )
        elif row["label"] == "DIFFERENT_LOCUS":
            if (
                UNKNOWN in {
                    row["source_truth_locus_id"], row["target_truth_locus_id"]
                }
                or row["source_truth_locus_id"] == row["target_truth_locus_id"]
            ):
                raise ContractError(
                    f"candidate truth row {line_number} DIFFERENT_LOCUS "
                    "requires two known distinct truth locus IDs"
                )
        if row["truth_confidence"] not in TRUTH_CONFIDENCE:
            raise ContractError(
                f"candidate truth row {line_number} has invalid truth_confidence "
                f"{row['truth_confidence']!r}"
            )
        if row["curation_blinded"] not in CURATION_BLINDED:
            raise ContractError(
                f"candidate truth row {line_number} has invalid curation_blinded "
                f"{row['curation_blinded']!r}"
            )
        for field in ("validation_method", "validation_source"):
            require_identifier(row[field], field, line_number)
        for field in ("taxon_a", "taxon_b"):
            if not row[field]:
                raise ContractError(
                    f"candidate truth row {line_number} has an empty {field}; "
                    f"use {UNKNOWN!r} for unknown"
                )

    if len(datasets) != 1:
        raise ContractError(
            "candidate truth must contain exactly one dataset_id; found "
            + ", ".join(sorted(datasets))
        )
    return next(iter(datasets)), rows


def load_locus_truth(
    path: Path, dataset_id: str
) -> tuple[list[dict[str, str]], dict[tuple[str, str], dict[str, str]]]:
    """Load the semantic TE-to-locus/context map used for leakage binding."""
    rows = read_tsv(path, LOCUS_TRUTH_COLUMNS, "locus truth TSV")
    forbidden = {
        "run_id", "candidate_id", "evidence_id", "evidence_group_id",
        "decision_id", "edge_id", "locus_id", "instance_id",
        "solver_component_id", "relation_id",
    }
    generated = forbidden & set(rows[0])
    if generated:
        raise ContractError(
            "locus truth contains generated prediction key(s): "
            + ", ".join(sorted(generated))
        )

    seen_record_ids: set[str] = set()
    by_te: dict[tuple[str, str], dict[str, str]] = {}
    locus_contexts: dict[str, tuple[str, str]] = {}
    event_loci: dict[str, str] = {}
    for line_number, row in enumerate(rows, start=2):
        if row["benchmark_schema_version"] != BENCHMARK_SCHEMA_VERSION:
            raise ContractError(
                f"locus truth row {line_number} has unsupported "
                "benchmark_schema_version "
                f"{row['benchmark_schema_version']!r}"
            )
        if row["dataset_id"] != dataset_id:
            raise ContractError(
                f"locus truth row {line_number} belongs to dataset "
                f"{row['dataset_id']!r}, expected {dataset_id!r}"
            )
        record_id = require_identifier(
            row["truth_record_id"], "truth_record_id", line_number
        )
        if record_id in seen_record_ids:
            raise ContractError(
                f"locus truth has duplicate truth_record_id {record_id!r}"
            )
        seen_record_ids.add(record_id)

        for field in ("truth_locus_id", "genome_id", "te_id"):
            require_identifier(row[field], field, line_number, "locus truth")
        for field in (
            "ancestral_event_id", "homology_group_id", "validation_batch_id"
        ):
            if not row[field]:
                raise ContractError(
                    f"locus truth row {line_number} has empty {field}; "
                    f"use {UNKNOWN!r} for unknown"
                )
        if row["ancestral_event_id"] == UNKNOWN:
            raise ContractError(
                f"locus truth row {line_number} has no ancestral_event_id"
            )
        if row["validation_batch_id"] == UNKNOWN:
            raise ContractError(
                f"locus truth row {line_number} has no validation_batch_id"
            )

        te_key = (row["genome_id"], row["te_id"])
        if te_key in by_te:
            raise ContractError(
                f"locus truth repeats semantic TE member {te_key!r}"
            )
        by_te[te_key] = row

        locus_id = row["truth_locus_id"]
        locus_context = (row["ancestral_event_id"], row["homology_group_id"])
        previous_context = locus_contexts.get(locus_id)
        if previous_context is not None and previous_context != locus_context:
            raise ContractError(
                f"locus truth assigns truth_locus_id {locus_id!r} to "
                "inconsistent ancestral_event_id/homology_group_id values"
            )
        locus_contexts[locus_id] = locus_context

        event_id = row["ancestral_event_id"]
        previous_locus = event_loci.get(event_id)
        if previous_locus is not None and previous_locus != locus_id:
            raise ContractError(
                f"locus truth assigns ancestral_event_id {event_id!r} to "
                f"multiple truth loci ({previous_locus!r}, {locus_id!r})"
            )
        event_loci[event_id] = locus_id

    return rows, by_te


def validate_candidate_locus_bindings(
    truth_rows: Sequence[dict[str, str]],
    locus_by_te: dict[tuple[str, str], dict[str, str]],
) -> None:
    """Cross-check candidate endpoint foreign keys against locus truth."""
    for row in truth_rows:
        record_id = row["truth_record_id"]
        endpoint_rows: list[dict[str, str]] = []
        for side in ("source", "target"):
            te_key = (row[f"{side}_genome_id"], row[f"{side}_te_id"])
            locus_row = locus_by_te.get(te_key)
            if locus_row is None:
                raise ContractError(
                    f"candidate truth record {record_id!r} {side} endpoint "
                    f"{te_key!r} is absent from locus truth"
                )
            expected_locus = row[f"{side}_truth_locus_id"]
            if locus_row["truth_locus_id"] != expected_locus:
                raise ContractError(
                    f"candidate truth record {record_id!r} {side} locus "
                    f"foreign key {expected_locus!r} disagrees with locus truth "
                    f"{locus_row['truth_locus_id']!r}"
                )
            endpoint_rows.append(locus_row)

        if row["label"] == "SAME_LOCUS":
            event_ids = {item["ancestral_event_id"] for item in endpoint_rows}
            if event_ids != {row["ancestral_event_id"]}:
                raise ContractError(
                    f"candidate truth record {record_id!r} ancestral_event_id "
                    "disagrees with locus truth"
                )


def load_splits(
    path: Path, dataset_id: str, truth_ids: set[str]
) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    rows = read_tsv(path, SPLIT_COLUMNS, "split TSV")
    by_id: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(rows, start=2):
        if row["split_schema_version"] != SPLIT_SCHEMA_VERSION:
            raise ContractError(
                f"split row {line_number} has unsupported split_schema_version "
                f"{row['split_schema_version']!r}"
            )
        if row["dataset_id"] != dataset_id:
            raise ContractError(
                f"split row {line_number} belongs to dataset "
                f"{row['dataset_id']!r}, expected {dataset_id!r}"
            )
        record_id = row["truth_record_id"]
        if not record_id or record_id == UNKNOWN:
            raise ContractError(f"split row {line_number} has no truth_record_id")
        if record_id in by_id:
            raise ContractError(
                f"split TSV has duplicate truth_record_id {record_id!r}"
            )
        if row["partition"] not in PARTITIONS:
            raise ContractError(
                f"split row {line_number} has invalid partition "
                f"{row['partition']!r}"
            )
        for field in ("outer_fold", "inner_fold"):
            if not row[field]:
                raise ContractError(
                    f"split row {line_number} has empty {field}; "
                    f"use {UNKNOWN!r} when not assigned"
                )
        partition = row["partition"]
        outer_fold = row["outer_fold"]
        inner_fold = row["inner_fold"]
        if partition == "DEVELOPMENT":
            if outer_fold == UNKNOWN or inner_fold == UNKNOWN:
                raise ContractError(
                    f"split row {line_number} DEVELOPMENT requires known "
                    "outer_fold and inner_fold"
                )
        elif partition == "CALIBRATION":
            if outer_fold == UNKNOWN:
                raise ContractError(
                    f"split row {line_number} CALIBRATION requires a known "
                    "outer_fold; inner_fold may be '.'"
                )
        elif outer_fold != UNKNOWN or inner_fold != UNKNOWN:
            raise ContractError(
                f"split row {line_number} EXTERNAL_TEST requires outer_fold "
                "and inner_fold to both be '.'"
            )
        by_id[record_id] = row

    split_ids = set(by_id)
    missing = sorted(truth_ids - split_ids)
    extra = sorted(split_ids - truth_ids)
    if missing or extra:
        details: list[str] = []
        if missing:
            details.append("missing truth_record_id(s): " + ", ".join(missing))
        if extra:
            details.append("unknown truth_record_id(s): " + ", ".join(extra))
        raise ContractError(
            "each candidate truth record must occur exactly once in the split TSV; "
            + "; ".join(details)
        )
    return rows, by_id


def encode_key(parts: Sequence[str]) -> str:
    """Return an unambiguous, deterministic representation of a compound ID."""
    return json.dumps(list(parts), ensure_ascii=False, separators=(",", ":"))


def canonical_te_pair(row: dict[str, str]) -> str:
    first = (row["source_genome_id"], row["source_te_id"])
    second = (row["target_genome_id"], row["target_te_id"])
    ordered = sorted((first, second))
    return encode_key(
        (ordered[0][0], ordered[0][1], ordered[1][0], ordered[1][1])
    )


def add_entity(
    entities: dict[tuple[str, str], set[str]],
    kind: str,
    value: str,
    record_id: str,
) -> None:
    entities[(kind, value)].add(record_id)


def binding_entities(
    truth_rows: Sequence[dict[str, str]],
    strict: bool,
    locus_by_te: dict[tuple[str, str], dict[str, str]] | None,
) -> tuple[dict[tuple[str, str], set[str]], list[dict[str, object]], list[dict[str, object]]]:
    entities: dict[tuple[str, str], set[str]] = defaultdict(set)
    violations: list[dict[str, object]] = []
    unverifiable: list[dict[str, object]] = []

    if locus_by_te is None:
        unverifiable.append({
            "code": "LOCUS_TRUTH_NOT_PROVIDED",
            "entity_type": "homology_group",
            "message": (
                "locus truth was not provided; endpoint homology-group, "
                "ancestral-event, and validation-batch leakage cannot be "
                "independently checked"
            ),
        })

    for row in truth_rows:
        record_id = row["truth_record_id"]
        dataset_id = row["dataset_id"]
        source_te = encode_key(
            (dataset_id, row["source_genome_id"], row["source_te_id"])
        )
        target_te = encode_key(
            (dataset_id, row["target_genome_id"], row["target_te_id"])
        )
        add_entity(entities, "semantic_te", source_te, record_id)
        add_entity(entities, "semantic_te", target_te, record_id)
        add_entity(
            entities,
            "reciprocal_semantic_pair",
            encode_key((dataset_id, canonical_te_pair(row))),
            record_id,
        )

        bindings = (
            ("ancestral_event", row["ancestral_event_id"]),
            ("truth_locus", row["source_truth_locus_id"]),
            ("truth_locus", row["target_truth_locus_id"]),
            ("validation_batch", row["validation_batch_id"]),
            ("clade_holdout", row["clade_holdout_id"]),
        )
        for kind, value in bindings:
            if not value:
                raise ContractError(
                    f"candidate truth record {record_id!r} has empty binding "
                    f"field for {kind}; use {UNKNOWN!r} for unknown"
                )
            if value == UNKNOWN:
                # A negative pair has two different ancestral events, so a
                # single shared event is genuinely not applicable. Endpoint
                # TE and locus bindings still keep both events leakage-safe.
                if kind == "ancestral_event" and row["label"] == "DIFFERENT_LOCUS":
                    continue
                item = {
                    "code": "UNKNOWN_BINDING_ENTITY",
                    "truth_record_id": record_id,
                    "entity_type": kind,
                    "message": f"{kind} is unknown and cannot be leakage-checked",
                }
                if strict:
                    violations.append({
                        **item,
                        "code": "STRICT_UNKNOWN_BINDING_ENTITY",
                    })
                else:
                    unverifiable.append(item)
                continue
            add_entity(entities, kind, encode_key((dataset_id, value)), record_id)

        if locus_by_te is not None:
            for side in ("source", "target"):
                te_key = (row[f"{side}_genome_id"], row[f"{side}_te_id"])
                # validate_candidate_locus_bindings() has already established
                # total coverage, so absence here is an internal contract bug.
                locus_row = locus_by_te[te_key]
                locus_bindings = (
                    ("ancestral_event", locus_row["ancestral_event_id"]),
                    ("homology_group", locus_row["homology_group_id"]),
                    ("validation_batch", locus_row["validation_batch_id"]),
                )
                for kind, value in locus_bindings:
                    if value == UNKNOWN:
                        item = {
                            "code": "UNKNOWN_LOCUS_TRUTH_BINDING_ENTITY",
                            "truth_record_id": record_id,
                            "endpoint": side,
                            "semantic_te": encode_key(
                                (dataset_id, te_key[0], te_key[1])
                            ),
                            "entity_type": kind,
                            "message": (
                                f"locus-truth {kind} is unknown and cannot be "
                                "leakage-checked"
                            ),
                        }
                        if strict:
                            violations.append({
                                **item,
                                "code": "STRICT_UNKNOWN_LOCUS_TRUTH_BINDING_ENTITY",
                            })
                        else:
                            unverifiable.append(item)
                        continue
                    add_entity(
                        entities, kind, encode_key((dataset_id, value)), record_id
                    )

    return entities, violations, unverifiable


def assignment_counts(
    split_rows: Sequence[dict[str, str]], field: str
) -> dict[str, object]:
    by_partition: dict[str, Counter[str]] = defaultdict(Counter)
    for row in split_rows:
        by_partition[row["partition"]][row[field]] += 1
    return {
        partition: dict(sorted(counts.items()))
        for partition, counts in sorted(by_partition.items())
    }


def check_entity_assignments(
    entities: dict[tuple[str, str], set[str]],
    splits: dict[str, dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    violations: list[dict[str, object]] = []
    unverifiable: list[dict[str, object]] = []
    for (kind, value), record_ids in sorted(entities.items()):
        if len(record_ids) < 2:
            continue
        ordered_ids = sorted(record_ids)
        partition_values: dict[str, list[str]] = defaultdict(list)
        for record_id in ordered_ids:
            partition_values[splits[record_id]["partition"]].append(record_id)
        if len(partition_values) > 1:
            violations.append({
                "code": "BINDING_ENTITY_CROSSES_SPLIT",
                "entity_type": kind,
                "entity_id": value,
                "dimension": "partition",
                "assignments": dict(sorted(partition_values.items())),
                "truth_record_ids": ordered_ids,
                "message": f"{kind} {value!r} crosses partition assignments",
            })

        # Fold identifiers are meaningful within a partition.  A dot is a
        # valid not-applicable value for CALIBRATION inner folds and for both
        # EXTERNAL_TEST folds, rather than an unverifiable assignment.
        for partition, partition_ids in sorted(partition_values.items()):
            if partition == "EXTERNAL_TEST":
                dimensions: tuple[str, ...] = ()
            elif partition == "CALIBRATION":
                dimensions = ("outer_fold", "inner_fold")
            else:
                dimensions = ("outer_fold", "inner_fold")
            for dimension in dimensions:
                values: dict[str, list[str]] = defaultdict(list)
                for record_id in partition_ids:
                    values[splits[record_id][dimension]].append(record_id)
                known = {
                    key: ids for key, ids in values.items() if key != UNKNOWN
                }
                unknown_ids = values.get(UNKNOWN, [])
                if len(known) > 1:
                    violations.append({
                        "code": "BINDING_ENTITY_CROSSES_SPLIT",
                        "entity_type": kind,
                        "entity_id": value,
                        "partition": partition,
                        "dimension": dimension,
                        "assignments": dict(sorted(known.items())),
                        "truth_record_ids": sorted(partition_ids),
                        "message": (
                            f"{kind} {value!r} crosses {partition} "
                            f"{dimension} assignments"
                        ),
                    })
                unknown_is_allowed = (
                    partition == "CALIBRATION" and dimension == "inner_fold"
                )
                if unknown_ids and not unknown_is_allowed:
                    unverifiable.append({
                        "code": "UNKNOWN_SPLIT_ASSIGNMENT",
                        "entity_type": kind,
                        "entity_id": value,
                        "partition": partition,
                        "dimension": dimension,
                        "truth_record_ids": sorted(unknown_ids),
                        "message": (
                            f"{kind} {value!r} has unknown {partition} "
                            f"{dimension} assignment"
                        ),
                    })
    return violations, unverifiable


def check_external_taxa(
    truth_rows: Sequence[dict[str, str]],
    splits: dict[str, dict[str, str]],
    allow_shared: bool,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    external: dict[str, set[str]] = defaultdict(set)
    nonexternal: dict[str, set[str]] = defaultdict(set)
    unverifiable: list[dict[str, object]] = []
    for row in truth_rows:
        record_id = row["truth_record_id"]
        destination = (
            external
            if splits[record_id]["partition"] == "EXTERNAL_TEST"
            else nonexternal
        )
        for field in ("taxon_a", "taxon_b"):
            taxon = row[field]
            if taxon == UNKNOWN:
                unverifiable.append({
                    "code": "UNKNOWN_TAXON",
                    "truth_record_id": record_id,
                    "field": field,
                    "message": "taxon is unknown and external isolation cannot be checked",
                })
            else:
                destination[taxon].add(record_id)

    violations: list[dict[str, object]] = []
    if allow_shared:
        return violations, unverifiable
    for taxon in sorted(set(external) & set(nonexternal)):
        violations.append({
            "code": "EXTERNAL_TAXON_LEAKAGE",
            "taxon": taxon,
            "external_truth_record_ids": sorted(external[taxon]),
            "nonexternal_truth_record_ids": sorted(nonexternal[taxon]),
            "message": (
                f"external-test taxon {taxon!r} is shared with a non-external partition"
            ),
        })
    return violations, unverifiable


def input_metadata(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "sha256": sha256_file(path)}


def base_report(
    truth_path: Path,
    splits_path: Path,
    locus_truth_path: Path | None,
    strict: bool,
    allow_shared_external_taxa: bool,
) -> dict[str, object]:
    inputs: dict[str, object] = {}
    input_paths: list[tuple[str, Path]] = [
        ("truth", truth_path), ("splits", splits_path)
    ]
    if locus_truth_path is not None:
        input_paths.append(("locus_truth", locus_truth_path))
    for name, path in input_paths:
        try:
            inputs[name] = input_metadata(path)
        except ContractError as error:
            inputs[name] = {
                "path": str(path.resolve()), "sha256": None,
                "error": str(error),
            }
    if locus_truth_path is None:
        inputs["locus_truth"] = {
            "provided": False,
            "path": None,
            "sha256": None,
        }
    return {
        "status": "EVALUATION_ONLY",
        "model_status": "NO_MODEL_FITTED",
        "benchmark_schema_version": BENCHMARK_SCHEMA_VERSION,
        "split_schema_version": SPLIT_SCHEMA_VERSION,
        "strict": strict,
        "allow_shared_external_taxa": allow_shared_external_taxa,
        "locus_truth_provided": locus_truth_path is not None,
        "inputs": inputs,
        "dataset_id": None,
        "counts": {
            "truth_records": 0,
            "locus_truth_records": 0,
            "split_records": 0,
            "partition_counts": {},
            "outer_fold_counts": {},
            "inner_fold_counts": {},
            "binding_entity_counts": {},
        },
        "violations": [],
        "unverifiable": [],
        "passed": False,
    }


def audit(
    truth_path: Path,
    splits_path: Path,
    locus_truth_path: Path | None,
    *,
    strict: bool,
    allow_shared_external_taxa: bool,
) -> dict[str, object]:
    report = base_report(
        truth_path, splits_path, locus_truth_path, strict,
        allow_shared_external_taxa
    )
    dataset_id, truth_rows = load_truth(truth_path)
    locus_rows: list[dict[str, str]] = []
    locus_by_te: dict[tuple[str, str], dict[str, str]] | None = None
    if locus_truth_path is not None:
        locus_rows, locus_by_te = load_locus_truth(locus_truth_path, dataset_id)
        validate_candidate_locus_bindings(truth_rows, locus_by_te)
    truth_ids = {row["truth_record_id"] for row in truth_rows}
    split_rows, splits = load_splits(splits_path, dataset_id, truth_ids)

    entities, violations, unverifiable = binding_entities(
        truth_rows, strict, locus_by_te
    )
    binding_violations, binding_unverifiable = check_entity_assignments(
        entities, splits
    )
    external_violations, external_unverifiable = check_external_taxa(
        truth_rows, splits, allow_shared_external_taxa
    )
    violations.extend(binding_violations)
    violations.extend(external_violations)
    unverifiable.extend(binding_unverifiable)
    unverifiable.extend(external_unverifiable)

    entity_counts = Counter(kind for kind, _ in entities)
    if locus_by_te is not None:
        # Keep the newly audited WGD binding visible even when every supplied
        # homology-group value is explicitly unknown.
        entity_counts.setdefault("homology_group", 0)
    report.update({
        "dataset_id": dataset_id,
        "counts": {
            "truth_records": len(truth_rows),
            "locus_truth_records": len(locus_rows),
            "split_records": len(split_rows),
            "partition_counts": dict(sorted(Counter(
                row["partition"] for row in split_rows
            ).items())),
            "outer_fold_counts": assignment_counts(split_rows, "outer_fold"),
            "inner_fold_counts": assignment_counts(split_rows, "inner_fold"),
            "binding_entity_counts": dict(sorted(entity_counts.items())),
        },
        "violations": violations,
        "unverifiable": unverifiable,
        "passed": not violations and (not strict or not unverifiable),
    })
    return report


def write_report(path: Path, report: dict[str, object]) -> None:
    temporary_path: Path | None = None
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
        )
        temporary_path = Path(temporary_name)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
        temporary_path = None
        directory_fd = os.open(
            path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
        )
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except OSError as error:
        raise SystemExit(f"cannot write report {path}: {error}") from error
    finally:
        if temporary_path is not None:
            try:
                temporary_path.unlink()
            except FileNotFoundError:
                pass


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "audit candidate-truth split assignments for biological and "
            "technical leakage; no model is fitted"
        )
    )
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument(
        "--locus-truth",
        type=Path,
        help=(
            "optional benchmark locus-member truth used to bind endpoint "
            "events, homology groups, and validation batches; required for "
            "a passing --strict audit"
        ),
    )
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "require locus truth and fail when event, locus, homology-group, "
            "validation-batch, or clade binding IDs are unknown"
        ),
    )
    parser.add_argument(
        "--allow-shared-external-taxa",
        action="store_true",
        help="permit taxa in EXTERNAL_TEST to also occur in other partitions",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        input_paths: list[tuple[str, Path]] = [
            ("truth", args.truth), ("splits", args.splits)
        ]
        if args.locus_truth is not None:
            input_paths.append(("locus_truth", args.locus_truth))
        require_output_disjoint_from_inputs(
            args.output,
            input_paths,
        )
    except ContractError as error:
        print(f"split audit refused unsafe output: {error}", file=sys.stderr)
        return 2
    report = base_report(
        args.truth,
        args.splits,
        args.locus_truth,
        args.strict,
        args.allow_shared_external_taxa,
    )
    try:
        report = audit(
            args.truth,
            args.splits,
            args.locus_truth,
            strict=args.strict,
            allow_shared_external_taxa=args.allow_shared_external_taxa,
        )
    except ContractError as error:
        report["violations"] = [{
            "code": "INPUT_CONTRACT_ERROR",
            "message": str(error),
        }]
        report["passed"] = False
    write_report(args.output, report)
    if not report["passed"]:
        print(
            f"split audit failed; see {args.output}",
            file=sys.stderr,
        )
        return 1
    print(f"split audit passed; report written to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
