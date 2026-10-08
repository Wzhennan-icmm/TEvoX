#!/usr/bin/env python3
"""Export leakage-aware, pre-decision TEvoX candidate features.

This command prepares a frozen evaluation/training table.  It deliberately
does not fit, tune, or calibrate a model, and it never treats the built-in
membership score as a probability.
"""

from __future__ import annotations

import gzip

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Sequence


BENCHMARK_SCHEMA_VERSION = "1.0.0"
EXPORT_SCHEMA_VERSION = "1.1.0"
EVIDENCE_SCHEMA_VERSION = "1.2.0"
EXPORT_STATUS = "EVALUATION_ONLY"
MODEL_FIT_STATUS = "NO_MODEL_FITTED"
CALIBRATION_STATUS = "UNCALIBRATED"
UNKNOWN = "."
SHA256_HEX_LENGTH = 64
PRINTED_FRACTION_TOLERANCE = 0.500001e-6

TRUTH_COLUMNS = (
    "benchmark_schema_version", "dataset_id", "truth_record_id",
    "source_genome_id", "source_te_id", "target_genome_id", "target_te_id",
    "label", "ancestral_event_id", "source_truth_locus_id",
    "target_truth_locus_id", "truth_confidence", "validation_method",
    "validation_source", "validation_batch_id", "curation_blinded",
    "taxon_a", "taxon_b", "clade_holdout_id",
)

EVIDENCE_COLUMNS = (
    "schema_version", "evidence_id", "evidence_group_id", "decision_id",
    "provider", "provider_path", "provider_record", "provider_line", "origin",
    "dependency", "query_genome_id", "target_genome_id", "source_te_id",
    "source_contig", "source_start", "source_end", "alignment_query_contig",
    "alignment_query_start", "alignment_query_end", "alignment_target_contig",
    "alignment_target_start", "alignment_target_end", "projection_contig",
    "projection_start", "projection_end", "alignment_identity",
    "alignment_identity_method", "alignment_error_count",
    "similarity_error_count", "nonalpha_count", "local_identity",
    "identity_method", "mapq", "mapq_status", "mapping_confidence",
    "left_flank", "right_flank", "left_flank_status", "right_flank_status",
    "te_aligned_fraction", "insertion_fraction", "target_n_fraction",
    "nearby_candidate_count", "retained_candidate_count",
    "graph_candidate_count", "eligible_candidate_count",
    "selected_target_te_id", "technical_state", "biological_state",
    "annotation_state", "legacy_state", "claim_type", "claimable", "quality",
    "evidence_completeness", "observation_rank", "primary", "near_best",
    "decision_ambiguous", "decision_code", "claimability_reason",
)

CANDIDATE_COLUMNS = (
    "schema_version", "candidate_id", "evidence_id", "decision_id",
    "target_genome_id", "target_te_id", "target_contig", "target_start",
    "target_end", "score", "reciprocal_overlap", "boundary_score",
    "breakpoint_distance", "family_relation", "context_relation",
    "shared_homology_group_id", "context_compatible", "candidate_rank",
    "graph_retained", "eligible", "selected", "decision_code",
)

FEATURE_COLUMNS = (
    "schema_version", "candidate_id", "evidence_id", "model_id",
    "calibration_status", "observed_feature_mask", "missing_features",
    "flank_min", "local_identity", "aggregate_identity", "mapq_normalized",
    "target_n_fraction", "te_aligned_fraction", "insertion_fraction",
    "reciprocal_overlap", "boundary_score", "family_relation",
    "context_relation", "membership_logit", "membership_score",
    "membership_entropy", "membership_prediction_set", "out_of_domain",
    "eligible",
)

RAW_FEATURE_COLUMNS = (
    "observed_feature_mask", "missing_features", "flank_min",
    "local_identity", "aggregate_identity", "mapq_normalized",
    "target_n_fraction", "te_aligned_fraction", "insertion_fraction",
    "reciprocal_overlap", "boundary_score", "family_relation",
    "context_relation",
)

TRAINING_COLUMNS = (
    "export_schema_version", "export_status", "model_fit_status",
    "calibration_status",
    "benchmark_schema_version", "dataset_id", "run_id", "truth_status",
    "truth_record_id", "label", "truth_confidence", "ancestral_event_id",
    "source_truth_locus_id", "target_truth_locus_id", "validation_method",
    "validation_source", "validation_batch_id", "curation_blinded", "taxon_a",
    "taxon_b", "clade_holdout_id", "candidate_generator_sha256",
    "inclusion_probability", "candidate_id", "evidence_id",
    "evidence_group_id", "provider", "origin", "dependency",
    "source_genome_id", "source_te_id", "source_te_semantic_key",
    "source_assembly_sha256", "source_annotation_sha256", "source_contig",
    "source_start", "source_end", "target_genome_id", "target_te_id",
    "target_te_semantic_key", "target_assembly_sha256",
    "target_annotation_sha256", "target_contig", "target_start", "target_end",
) + RAW_FEATURE_COLUMNS

UNMATCHED_COLUMNS = TRUTH_COLUMNS + ("unmatched_reason",)

LABELS = {"SAME_LOCUS", "DIFFERENT_LOCUS", "UNKNOWN"}
TRUTH_CONFIDENCE = {"HIGH", "MEDIUM", "LOW"}
CURATION_BLINDED = {"true", "false", "unknown"}
FAMILY_RELATIONS = {"MATCH", "CONFLICT", "UNKNOWN"}
CONTEXT_RELATIONS = {"SUPPORTED", "CONFLICT", "AMBIGUOUS", "UNKNOWN"}
BOOL_VALUES = {"true", "false"}
KNOWN_MISSING_FEATURES = {
    "flank_min", "identity", "mapq", "n_fraction", "family", "context"
}

# Public bits from include/tevox.h.  Pinning them here makes malformed or
# schema-incompatible feature sidecars fail closed.
FEATURE_BITS = {
    "flank_min": 1 << 0,
    "local_identity": 1 << 1,
    "aggregate_identity": 1 << 2,
    "mapq": 1 << 3,
    "n_fraction": 1 << 4,
    "te_alignment": 1 << 5,
    "insertion": 1 << 6,
    "reciprocal_overlap": 1 << 7,
    "boundary": 1 << 8,
    "family": 1 << 9,
    "context": 1 << 10,
}
ALL_FEATURE_BITS = sum(FEATURE_BITS.values())


class ContractError(ValueError):
    """An input cannot safely be exported under this contract."""


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


def fail_constant(value: str) -> None:
    raise ContractError(f"JSON contains non-finite constant {value!r}")


def read_json(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(handle, parse_constant=fail_constant)
    except (OSError, json.JSONDecodeError) as error:
        raise ContractError(f"cannot read valid JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ContractError(f"run metadata {path} must be a JSON object")
    validate_json_value(value, "run metadata")
    return value


def validate_json_value(value: Any, label: str) -> None:
    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ContractError(f"{label} contains a non-finite number")
        return
    if isinstance(value, list):
        for item in value:
            validate_json_value(item, label)
        return
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise ContractError(f"{label} contains a non-string object key")
        for item in value.values():
            validate_json_value(item, label)
        return
    raise ContractError(f"{label} contains unsupported JSON value {type(value).__name__}")


def canonical_json(value: Any) -> str:
    validate_json_value(value, "canonical JSON")
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
        allow_nan=False,
    )


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            while block := handle.read(1024 * 1024):
                digest.update(block)
    except OSError as error:
        raise ContractError(f"cannot hash input {path}: {error}") from error
    return digest.hexdigest()


def validate_sha256(value: Any, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != SHA256_HEX_LENGTH
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ContractError(f"{label} must be a 64-character lowercase SHA-256")
    return value


def read_tsv(path: Path, expected: Sequence[str], label: str) -> list[dict[str, str]]:
    try:
        with open_table(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames
            if fields is None:
                raise ContractError(f"{label} {path} has no header")
            if any(field is None or not field or field != field.strip() for field in fields):
                raise ContractError(f"{label} {path} has an empty or malformed header field")
            duplicates = sorted(name for name, count in Counter(fields).items() if count > 1)
            if duplicates:
                raise ContractError(
                    f"{label} {path} has duplicate columns: {', '.join(duplicates)}"
                )
            missing = sorted(set(expected) - set(fields))
            extra = sorted(set(fields) - set(expected))
            if missing or extra:
                details: list[str] = []
                if missing:
                    details.append("missing " + ", ".join(missing))
                if extra:
                    details.append("unexpected " + ", ".join(extra))
                raise ContractError(f"{label} {path} columns do not match contract: " + "; ".join(details))
            rows: list[dict[str, str]] = []
            for line, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise ContractError(f"{label} {path}:{line} has the wrong number of columns")
                if not any(value for value in row.values()):
                    raise ContractError(f"{label} {path}:{line} is empty")
                empty = [field for field, value in row.items() if value == ""]
                if empty:
                    raise ContractError(
                        f"{label} {path}:{line} has empty field(s) "
                        f"{', '.join(empty)}; use '.' for missing values"
                    )
                rows.append(row)
    except OSError as error:
        raise ContractError(f"cannot read {label} {path}: {error}") from error
    if not rows:
        raise ContractError(f"{label} {path} has no data rows")
    return rows


def clean(value: str, field: str, line: int, *, allow_unknown: bool = False) -> str:
    stripped = value.strip()
    if value != stripped:
        raise ContractError(f"row {line} field {field} has leading or trailing whitespace")
    if not stripped or (stripped == UNKNOWN and not allow_unknown):
        raise ContractError(f"row {line} field {field} is missing")
    if any(character in stripped for character in "\t\r\n"):
        raise ContractError(f"row {line} field {field} contains a control delimiter")
    return stripped


def unique_index(
    rows: Sequence[dict[str, str]], field: str, label: str
) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for line, row in enumerate(rows, start=2):
        key = clean(row[field], field, line)
        if key in result:
            raise ContractError(f"{label} repeats {field} {key!r}")
        result[key] = row
    return result


def finite_number(
    value: str, field: str, line: int, *, optional: bool = False,
    minimum: float | None = None, maximum: float | None = None,
) -> float | None:
    if optional and value == UNKNOWN:
        return None
    if value == "":
        raise ContractError(f"row {line} field {field} is empty; use '.' for missing")
    try:
        number = float(value)
    except ValueError as error:
        raise ContractError(f"row {line} field {field} is not numeric: {value!r}") from error
    if not math.isfinite(number):
        raise ContractError(f"row {line} field {field} is not finite")
    if minimum is not None and number < minimum:
        raise ContractError(f"row {line} field {field} is below {minimum}")
    if maximum is not None and number > maximum:
        raise ContractError(f"row {line} field {field} is above {maximum}")
    return number


def optional_fraction(value: str, field: str, line: int) -> float | None:
    return finite_number(
        value, field, line, optional=True, minimum=0.0, maximum=1.0
    )


def require_printed_fraction_equal(
    observed: str, expected: str | float | None, field: str, line: int,
) -> None:
    """Compare values emitted to six decimals, preserving missingness."""
    observed_value = optional_fraction(observed, field, line)
    if isinstance(expected, str):
        expected_value = optional_fraction(expected, field, line)
    else:
        expected_value = expected
        if expected_value is not None and (
            not math.isfinite(expected_value)
            or expected_value < 0.0
            or expected_value > 1.0
        ):
            raise ContractError(f"row {line} reconstructed {field} is out of range")
    if (observed_value is None) != (expected_value is None):
        raise ContractError(
            f"feature row {line} {field} missingness disagrees with source output"
        )
    if (
        observed_value is not None
        and expected_value is not None
        and abs(observed_value - expected_value) > PRINTED_FRACTION_TOLERANCE
    ):
        raise ContractError(
            f"feature row {line} {field}={observed!r} disagrees with "
            f"reconstructed value {expected_value:.9g}"
        )


def integer(value: str, field: str, line: int, *, minimum: int = 0) -> int:
    try:
        number = int(value)
    except ValueError as error:
        raise ContractError(f"row {line} field {field} is not an integer: {value!r}") from error
    if str(number) != value or number < minimum:
        raise ContractError(f"row {line} field {field} is not a canonical integer >= {minimum}")
    return number


def validate_identifier(value: str, label: str) -> str:
    result = value.strip()
    if value != result or not result or result == UNKNOWN:
        raise ContractError(f"{label} must be a non-empty identifier without surrounding whitespace")
    if any(character in result for character in "\t\r\n"):
        raise ContractError(f"{label} contains a control delimiter")
    return result


def canonical_pair(
    genome_a: str, te_a: str, genome_b: str, te_b: str
) -> tuple[tuple[str, str], tuple[str, str]]:
    a = (genome_a, te_a)
    b = (genome_b, te_b)
    if a == b:
        raise ContractError(f"semantic pair repeats the same TE {a!r}")
    return (a, b) if a < b else (b, a)


def orient_truth(
    row: dict[str, str], source_genome: str, source_te: str,
    target_genome: str, target_te: str,
) -> dict[str, str]:
    """Return endpoint metadata in the generated candidate's direction."""
    truth_source = (row["source_genome_id"], row["source_te_id"])
    truth_target = (row["target_genome_id"], row["target_te_id"])
    candidate_source = (source_genome, source_te)
    candidate_target = (target_genome, target_te)
    if (candidate_source, candidate_target) == (truth_source, truth_target):
        return dict(row)
    if (candidate_source, candidate_target) != (truth_target, truth_source):
        raise ContractError(
            "candidate orientation does not match its semantic truth pair"
        )
    oriented = dict(row)
    # These fields describe pair endpoints rather than the unordered event.
    # Keep this list explicit so future endpoint fields cannot be silently
    # emitted in the truth file's direction for a derived-reverse candidate.
    for source_field, target_field in (
        ("source_genome_id", "target_genome_id"),
        ("source_te_id", "target_te_id"),
        ("source_truth_locus_id", "target_truth_locus_id"),
        ("taxon_a", "taxon_b"),
    ):
        oriented[source_field], oriented[target_field] = (
            oriented[target_field], oriented[source_field]
        )
    return oriented


def validate_truth(
    rows: Sequence[dict[str, str]], dataset_id: str
) -> dict[tuple[tuple[str, str], tuple[str, str]], dict[str, str]]:
    by_record = unique_index(rows, "truth_record_id", "truth TSV")
    del by_record
    result: dict[tuple[tuple[str, str], tuple[str, str]], dict[str, str]] = {}
    te_loci: dict[tuple[str, str], str] = {}
    genome_taxa: dict[str, str] = {}
    locus_events: dict[str, str] = {}
    for line, row in enumerate(rows, start=2):
        for field in TRUTH_COLUMNS:
            clean(row[field], field, line, allow_unknown=field in {
                "ancestral_event_id", "source_truth_locus_id",
                "target_truth_locus_id", "taxon_a", "taxon_b",
                "clade_holdout_id",
            })
        if row["benchmark_schema_version"] != BENCHMARK_SCHEMA_VERSION:
            raise ContractError(
                f"truth row {line} uses unsupported benchmark schema "
                f"{row['benchmark_schema_version']!r}"
            )
        if row["dataset_id"] != dataset_id:
            raise ContractError(
                f"truth row {line} belongs to dataset {row['dataset_id']!r}, "
                f"expected {dataset_id!r}"
            )
        if row["label"] not in LABELS:
            raise ContractError(f"truth row {line} has invalid label {row['label']!r}")
        if row["truth_confidence"] not in TRUTH_CONFIDENCE:
            raise ContractError(
                f"truth row {line} has invalid truth_confidence "
                f"{row['truth_confidence']!r}"
            )
        if row["curation_blinded"] not in CURATION_BLINDED:
            raise ContractError(
                f"truth row {line} has invalid curation_blinded "
                f"{row['curation_blinded']!r}"
            )
        source_locus = row["source_truth_locus_id"]
        target_locus = row["target_truth_locus_id"]
        if row["label"] == "SAME_LOCUS":
            if row["ancestral_event_id"] == UNKNOWN:
                raise ContractError(
                    f"truth row {line} SAME_LOCUS requires ancestral_event_id"
                )
            if (
                UNKNOWN in {source_locus, target_locus}
                or source_locus != target_locus
            ):
                raise ContractError(
                    f"truth row {line} SAME_LOCUS requires one known shared truth locus ID"
                )
        if row["label"] == "DIFFERENT_LOCUS":
            if UNKNOWN in {source_locus, target_locus} or source_locus == target_locus:
                raise ContractError(
                    f"truth row {line} DIFFERENT_LOCUS requires two distinct truth locus IDs"
                )

        endpoints = (
            (
                row["source_genome_id"], row["source_te_id"], source_locus,
                row["taxon_a"],
            ),
            (
                row["target_genome_id"], row["target_te_id"], target_locus,
                row["taxon_b"],
            ),
        )
        for genome_id, te_id, locus_id, taxon in endpoints:
            previous_taxon = genome_taxa.setdefault(genome_id, taxon)
            if previous_taxon != taxon:
                raise ContractError(
                    f"truth row {line} maps genome {genome_id!r} to both "
                    f"{previous_taxon!r} and {taxon!r}"
                )
            if locus_id != UNKNOWN:
                te_key = (genome_id, te_id)
                previous_locus = te_loci.setdefault(te_key, locus_id)
                if previous_locus != locus_id:
                    raise ContractError(
                        f"truth row {line} maps TE {te_key!r} to both "
                        f"{previous_locus!r} and {locus_id!r}"
                    )
                event_id = row["ancestral_event_id"]
                if event_id != UNKNOWN:
                    previous_event = locus_events.setdefault(locus_id, event_id)
                    if previous_event != event_id:
                        raise ContractError(
                            f"truth row {line} maps truth locus {locus_id!r} "
                            f"to both {previous_event!r} and {event_id!r}"
                        )
        key = canonical_pair(
            row["source_genome_id"], row["source_te_id"],
            row["target_genome_id"], row["target_te_id"],
        )
        if key in result:
            raise ContractError(f"truth TSV repeats semantic pair {key!r}")
        result[key] = row
    return result


def require_object(parent: dict[str, Any], field: str) -> dict[str, Any]:
    value = parent.get(field)
    if not isinstance(value, dict):
        raise ContractError(f"run metadata field {field!r} must be an object")
    return value


def require_list(parent: dict[str, Any], field: str) -> list[Any]:
    value = parent.get(field)
    if not isinstance(value, list):
        raise ContractError(f"run metadata field {field!r} must be an array")
    return value


def validate_run(run: dict[str, Any]) -> tuple[dict[str, Any], str, str, str]:
    if run.get("software") != "TEvoX":
        raise ContractError("run metadata software must be 'TEvoX'")
    version = run.get("version")
    schema = run.get("schema_version")
    coordinate_system = run.get("coordinate_system")
    if not isinstance(version, str) or not version:
        raise ContractError("run metadata version must be a non-empty string")
    if not isinstance(schema, str) or not schema:
        raise ContractError("run metadata schema_version must be a non-empty string")
    if schema != EVIDENCE_SCHEMA_VERSION:
        raise ContractError(
            f"training exporter supports evidence schema {EVIDENCE_SCHEMA_VERSION}, "
            f"found {schema!r}"
        )
    if coordinate_system != "0-based-half-open":
        raise ContractError("run metadata must use 0-based-half-open coordinates")
    model = require_object(run, "inference_model")
    if model.get("calibration_status") != CALIBRATION_STATUS:
        raise ContractError("training export requires run calibration_status=UNCALIBRATED")
    model_id = model.get("model_id")
    if not isinstance(model_id, str) or not model_id:
        raise ContractError("run metadata inference_model.model_id is missing")
    config = require_object(run, "config")
    candidate_window = config.get("candidate_window")
    if (
        isinstance(candidate_window, bool)
        or not isinstance(candidate_window, int)
        or candidate_window < 0
    ):
        raise ContractError(
            "run config candidate_window must be a non-negative integer"
        )
    for field in ("max_candidates", "max_graph_candidates"):
        value = config.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value != 0:
            raise ContractError(
                f"training export requires run config {field}=0; found {value!r}"
            )
    genomes = require_list(run, "genomes")
    if not genomes:
        raise ContractError("run metadata genomes array is empty")
    return config, version, schema, model_id


def resolve_recorded_path(text: str, metadata_path: Path) -> Path:
    value = Path(text)
    if value.is_absolute():
        candidate = value.resolve()
        if not candidate.is_file():
            raise ContractError(f"recorded input is not a file: {candidate}")
        return candidate
    possibilities = {
        (Path.cwd() / value).resolve(),
        (metadata_path.parent / value).resolve(),
    }
    existing = sorted((path for path in possibilities if path.is_file()), key=str)
    if not existing:
        raise ContractError(
            f"cannot resolve recorded relative input {text!r} from cwd or run directory"
        )
    if len(existing) > 1 and existing[0] != existing[1]:
        raise ContractError(f"recorded relative input {text!r} is ambiguous")
    return existing[0]


def load_genome_hashes(
    run: dict[str, Any], run_path: Path
) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    hashes: dict[str, dict[str, str]] = {}
    inputs: list[dict[str, str]] = []
    for index, item in enumerate(require_list(run, "genomes")):
        if not isinstance(item, dict):
            raise ContractError(f"run metadata genome {index} is not an object")
        genome_id = item.get("genome_id")
        fasta = item.get("fasta")
        annotation = item.get("te_annotation")
        if not all(isinstance(value, str) and value for value in (genome_id, fasta, annotation)):
            raise ContractError(f"run metadata genome {index} has missing identifiers or paths")
        assert isinstance(genome_id, str) and isinstance(fasta, str) and isinstance(annotation, str)
        validate_identifier(genome_id, f"genomes[{index}].genome_id")
        if genome_id in hashes:
            raise ContractError(f"run metadata repeats genome_id {genome_id!r}")
        fasta_path = resolve_recorded_path(fasta, run_path)
        annotation_path = resolve_recorded_path(annotation, run_path)
        fasta_sha = sha256_file(fasta_path)
        annotation_sha = sha256_file(annotation_path)
        recorded_fasta_sha = validate_sha256(
            item.get("fasta_sha256"), f"genomes[{index}].fasta_sha256"
        )
        recorded_annotation_sha = validate_sha256(
            item.get("te_annotation_sha256"),
            f"genomes[{index}].te_annotation_sha256",
        )
        if recorded_fasta_sha != fasta_sha:
            raise ContractError(
                f"genomes[{index}] FASTA hash disagrees with current input {fasta_path}"
            )
        if recorded_annotation_sha != annotation_sha:
            raise ContractError(
                f"genomes[{index}] TE annotation hash disagrees with current input "
                f"{annotation_path}"
            )
        hashes[genome_id] = {
            "assembly_sha256": fasta_sha,
            "annotation_sha256": annotation_sha,
            "fasta_path": str(fasta_path),
            "annotation_path": str(annotation_path),
        }
        inputs.extend((
            {"role": "assembly", "genome_id": genome_id,
             "path": str(fasta_path), "sha256": fasta_sha},
            {"role": "annotation", "genome_id": genome_id,
             "path": str(annotation_path), "sha256": annotation_sha},
        ))
    return hashes, inputs


def load_frozen_input_inventory(
    run: dict[str, Any], run_path: Path
) -> tuple[list[dict[str, str]], dict[str, str]]:
    values = require_list(run, "input_files")
    if not values:
        raise ContractError("run metadata input_files inventory is empty")
    counts = require_object(run, "counts")
    if counts.get("input_files") != len(values):
        raise ContractError("run metadata input_files count is inconsistent")
    inputs: list[dict[str, str]] = []
    by_path: dict[str, str] = {}
    seen: set[tuple[str, str]] = set()
    for index, item in enumerate(values):
        if not isinstance(item, dict):
            raise ContractError(f"run metadata input_files[{index}] is not an object")
        role = item.get("role")
        path_text = item.get("path")
        if not isinstance(role, str) or not role:
            raise ContractError(f"run metadata input_files[{index}] has no role")
        if not isinstance(path_text, str) or not path_text:
            raise ContractError(f"run metadata input_files[{index}] has no path")
        path = resolve_recorded_path(path_text, run_path)
        recorded_sha = validate_sha256(
            item.get("sha256"), f"input_files[{index}].sha256"
        )
        current_sha = sha256_file(path)
        if recorded_sha != current_sha:
            raise ContractError(
                f"run metadata input_files[{index}] hash disagrees with "
                f"current input {path}"
            )
        path_key = str(path)
        previous = by_path.setdefault(path_key, current_sha)
        if previous != current_sha:
            raise ContractError(f"run metadata has conflicting hashes for {path}")
        key = (role, path_key)
        if key in seen:
            raise ContractError(
                f"run metadata repeats input role/path {role!r}, {path_key!r}"
            )
        seen.add(key)
        inputs.append({"role": role, "path": path_key, "sha256": current_sha})
    return inputs, by_path


def provider_contract(
    run: dict[str, Any], evidence_rows: Sequence[dict[str, str]],
    run_path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    providers = {clean(row["provider"], "provider", line)
                 for line, row in enumerate(evidence_rows, start=2)}
    inputs: list[dict[str, str]] = []
    seen_inputs: set[tuple[str, str]] = set()
    path_hashes: dict[Path, str] = {}
    for section in ("alignment_evidence", "synteny_evidence"):
        for index, item in enumerate(require_list(run, section)):
            if not isinstance(item, dict):
                raise ContractError(f"run metadata {section}[{index}] is not an object")
            provider = item.get("provider")
            path_text = item.get("path")
            if not isinstance(provider, str) or not provider:
                raise ContractError(f"run metadata {section}[{index}] has no provider")
            if not isinstance(path_text, str) or not path_text:
                raise ContractError(f"run metadata {section}[{index}] has no path")
            providers.add(provider)
            path = resolve_recorded_path(path_text, run_path)
            recorded_sha = validate_sha256(
                item.get("path_sha256"), f"{section}[{index}].path_sha256"
            )
            current_sha = path_hashes.get(path)
            if current_sha is None:
                current_sha = sha256_file(path)
                path_hashes[path] = current_sha
            if recorded_sha != current_sha:
                raise ContractError(
                    f"run metadata {section}[{index}] hash disagrees with "
                    f"current provider input {path}"
                )
            key = (provider, str(path))
            if key in seen_inputs:
                continue
            seen_inputs.add(key)
            inputs.append({
                "role": "provider_input", "provider": provider,
                "path": str(path), "sha256": current_sha,
            })
    for index, item in enumerate(require_list(run, "synteny_inputs")):
        if not isinstance(item, dict):
            raise ContractError(
                f"run metadata synteny_inputs[{index}] is not an object"
            )
        provider = item.get("provider")
        if provider != "MCScanX":
            raise ContractError(
                f"run metadata synteny_inputs[{index}] has unsupported provider"
            )
        providers.add(provider)
        for path_field, hash_field in (
            ("collinearity_path", "collinearity_sha256"),
            ("gene_table_path", "gene_table_sha256"),
        ):
            path_text = item.get(path_field)
            if not isinstance(path_text, str) or not path_text:
                raise ContractError(
                    f"run metadata synteny_inputs[{index}] has no {path_field}"
                )
            path = resolve_recorded_path(path_text, run_path)
            recorded_sha = validate_sha256(
                item.get(hash_field),
                f"synteny_inputs[{index}].{hash_field}",
            )
            current_sha = path_hashes.get(path)
            if current_sha is None:
                current_sha = sha256_file(path)
                path_hashes[path] = current_sha
            if recorded_sha != current_sha:
                raise ContractError(
                    f"run metadata synteny_inputs[{index}] {path_field} hash "
                    f"disagrees with current provider input {path}"
                )
            key = (provider, str(path))
            if key in seen_inputs:
                continue
            seen_inputs.add(key)
            inputs.append({
                "role": "provider_input", "provider": provider,
                "path": str(path), "sha256": current_sha,
            })
    if not providers:
        raise ContractError("candidate generator provider set is empty")
    inputs.sort(key=lambda item: (item["provider"], item["path"]))
    return sorted(providers), inputs


def validate_schema_versions(
    rows: Sequence[dict[str, str]], expected: str, label: str
) -> None:
    for line, row in enumerate(rows, start=2):
        if row["schema_version"] != expected:
            raise ContractError(
                f"{label} row {line} schema {row['schema_version']!r} "
                f"does not match run schema {expected!r}"
            )


def validate_coordinates(row: dict[str, str], prefix: str, line: int) -> tuple[int, int]:
    start = integer(row[f"{prefix}_start"], f"{prefix}_start", line)
    end = integer(row[f"{prefix}_end"], f"{prefix}_end", line)
    if end <= start:
        raise ContractError(f"row {line} {prefix} interval must have positive length")
    return start, end


def validate_features(row: dict[str, str], line: int) -> None:
    try:
        mask = int(row["observed_feature_mask"], 16)
    except ValueError as error:
        raise ContractError(
            f"feature row {line} has invalid observed_feature_mask "
            f"{row['observed_feature_mask']!r}"
        ) from error
    if row["observed_feature_mask"] != f"0x{mask:08x}" or mask & ~ALL_FEATURE_BITS:
        raise ContractError(f"feature row {line} has non-canonical or unknown feature mask")
    missing_text = row["missing_features"]
    if missing_text == UNKNOWN:
        missing: set[str] = set()
    else:
        parts = missing_text.split(",")
        if len(parts) != len(set(parts)) or set(parts) - KNOWN_MISSING_FEATURES:
            raise ContractError(f"feature row {line} has invalid missing_features")
        missing = set(parts)
    expected_missing = set()
    if not mask & FEATURE_BITS["flank_min"]:
        expected_missing.add("flank_min")
    if not mask & (FEATURE_BITS["local_identity"] | FEATURE_BITS["aggregate_identity"]):
        expected_missing.add("identity")
    if not mask & FEATURE_BITS["mapq"]:
        expected_missing.add("mapq")
    if not mask & FEATURE_BITS["n_fraction"]:
        expected_missing.add("n_fraction")
    if not mask & FEATURE_BITS["family"]:
        expected_missing.add("family")
    if not mask & FEATURE_BITS["context"]:
        expected_missing.add("context")
    if missing != expected_missing:
        raise ContractError(
            f"feature row {line} missing_features disagrees with observed_feature_mask"
        )
    required_core = (
        FEATURE_BITS["te_alignment"] | FEATURE_BITS["insertion"]
        | FEATURE_BITS["reciprocal_overlap"] | FEATURE_BITS["boundary"]
    )
    if mask & required_core != required_core:
        raise ContractError(f"feature row {line} lacks required pre-decision core features")

    fractions = (
        "flank_min", "local_identity", "aggregate_identity", "mapq_normalized",
        "target_n_fraction", "te_aligned_fraction", "insertion_fraction",
        "reciprocal_overlap", "boundary_score",
    )
    values = {
        field: finite_number(row[field], field, line, optional=True, minimum=0.0, maximum=1.0)
        for field in fractions
    }
    expected_observed = {
        "flank_min": bool(mask & FEATURE_BITS["flank_min"]),
        "local_identity": bool(mask & FEATURE_BITS["local_identity"]),
        "mapq_normalized": bool(mask & FEATURE_BITS["mapq"]),
        "target_n_fraction": bool(mask & FEATURE_BITS["n_fraction"]),
    }
    for field, observed in expected_observed.items():
        if observed != (values[field] is not None):
            raise ContractError(f"feature row {line} mask disagrees with {field}")
    local_bit = bool(mask & FEATURE_BITS["local_identity"])
    aggregate_bit = bool(mask & FEATURE_BITS["aggregate_identity"])
    if local_bit and aggregate_bit:
        raise ContractError(f"feature row {line} selects both local and fallback identity")
    # Aggregate identity is printed whenever the provider supplies it, but its
    # bit is selected only when local identity is absent and aggregate identity
    # is therefore the feature actually consumed by the built-in scorer.
    if aggregate_bit and values["aggregate_identity"] is None:
        raise ContractError(f"feature row {line} selects missing aggregate_identity")
    if not local_bit and values["aggregate_identity"] is not None and not aggregate_bit:
        raise ContractError(f"feature row {line} fails to select available fallback identity")
    for field in ("te_aligned_fraction", "insertion_fraction", "reciprocal_overlap", "boundary_score"):
        if values[field] is None:
            raise ContractError(f"feature row {line} required feature {field} is missing")
    if row["family_relation"] not in FAMILY_RELATIONS:
        raise ContractError(f"feature row {line} has invalid family_relation")
    if row["context_relation"] not in CONTEXT_RELATIONS:
        raise ContractError(f"feature row {line} has invalid context_relation")
    if bool(mask & FEATURE_BITS["family"]) != (row["family_relation"] != "UNKNOWN"):
        raise ContractError(f"feature row {line} mask disagrees with family_relation")
    if bool(mask & FEATURE_BITS["context"]) != (row["context_relation"] != "UNKNOWN"):
        raise ContractError(f"feature row {line} mask disagrees with context_relation")
    finite_number(row["membership_logit"], "membership_logit", line)
    finite_number(row["membership_score"], "membership_score", line, minimum=0.0, maximum=1.0)
    finite_number(row["membership_entropy"], "membership_entropy", line, minimum=0.0, maximum=1.0)
    if row["calibration_status"] != CALIBRATION_STATUS:
        raise ContractError(f"feature row {line} is not UNCALIBRATED")
    if row["eligible"] not in BOOL_VALUES or row["out_of_domain"] not in BOOL_VALUES:
        raise ContractError(f"feature row {line} has invalid boolean")


def projection_interval(row: dict[str, str], line: int) -> tuple[int, int]:
    start = integer(row["projection_start"], "projection_start", line)
    end = integer(row["projection_end"], "projection_end", line)
    if end < start:
        raise ContractError(f"evidence row {line} projection interval is reversed")
    return start, end


def expected_candidate_geometry(
    observation: dict[str, str], candidate: dict[str, str],
    candidate_window: int, line: int,
) -> tuple[float, int, float]:
    projection_start, projection_end = projection_interval(observation, line)
    target_start, target_end = validate_coordinates(candidate, "target", line)
    projected_length = projection_end - projection_start
    target_length = target_end - target_start
    shared = max(
        0, min(projection_end, target_end) - max(projection_start, target_start)
    )
    reciprocal = 0.0
    if projected_length > 0 and target_length > 0:
        reciprocal = min(shared / projected_length, shared / target_length)

    start_distance = abs(projection_start - target_start)
    end_distance = abs(projection_end - target_end)
    int32_max = (1 << 31) - 1
    if (
        start_distance >= int32_max
        or end_distance >= int32_max
        or start_distance + end_distance > int32_max
    ):
        distance = int32_max
    else:
        distance = start_distance + end_distance
    boundary = 1.0 - min(distance / (2 * candidate_window + 1), 1.0)
    return reciprocal, distance, boundary


def validate_raw_feature_sources(
    feature: dict[str, str], candidate: dict[str, str],
    observation: dict[str, str], candidate_window: int, line: int,
) -> None:
    """Rebuild every exported raw feature from evidence/candidate outputs."""
    if observation["target_genome_id"] != candidate["target_genome_id"]:
        raise ContractError(
            f"feature row {line} candidate target genome disagrees with evidence"
        )

    flank_statuses = (
        observation["left_flank_status"], observation["right_flank_status"]
    )
    if any(status not in {"OBSERVED", "CONTIG_EDGE"} for status in flank_statuses):
        raise ContractError(f"evidence for feature row {line} has invalid flank status")
    left_flank = optional_fraction(observation["left_flank"], "left_flank", line)
    right_flank = optional_fraction(observation["right_flank"], "right_flank", line)
    if flank_statuses[0] == "OBSERVED" and left_flank is None:
        raise ContractError(f"evidence for feature row {line} has missing observed left flank")
    if flank_statuses[1] == "OBSERVED" and right_flank is None:
        raise ContractError(f"evidence for feature row {line} has missing observed right flank")
    flank_min: float | None = None
    if flank_statuses == ("OBSERVED", "OBSERVED"):
        assert left_flank is not None and right_flank is not None
        flank_min = min(left_flank, right_flank)

    local_identity = optional_fraction(
        observation["local_identity"], "local_identity", line
    )
    aggregate_identity = optional_fraction(
        observation["alignment_identity"], "alignment_identity", line
    )
    if observation["mapq"] == UNKNOWN:
        mapq_normalized: float | None = None
    else:
        mapq = integer(observation["mapq"], "mapq", line)
        if mapq > 255:
            raise ContractError(f"evidence for feature row {line} has mapq above 255")
        mapq_normalized = min(1.0, mapq / 60.0)
    target_n_fraction = optional_fraction(
        observation["target_n_fraction"], "target_n_fraction", line
    )
    te_aligned_fraction = optional_fraction(
        observation["te_aligned_fraction"], "te_aligned_fraction", line
    )
    insertion_fraction = optional_fraction(
        observation["insertion_fraction"], "insertion_fraction", line
    )
    if te_aligned_fraction is None or insertion_fraction is None:
        raise ContractError(
            f"evidence for feature row {line} lacks required TE/insertion fractions"
        )

    reciprocal, breakpoint_distance, boundary = expected_candidate_geometry(
        observation, candidate, candidate_window, line
    )
    reported_distance = integer(
        candidate["breakpoint_distance"], "breakpoint_distance", line
    )
    if reported_distance != breakpoint_distance:
        raise ContractError(
            f"candidate for feature row {line} breakpoint_distance="
            f"{reported_distance} disagrees with reconstructed {breakpoint_distance}"
        )
    require_printed_fraction_equal(
        candidate["reciprocal_overlap"], reciprocal,
        "candidate.reciprocal_overlap", line,
    )
    require_printed_fraction_equal(
        candidate["boundary_score"], boundary, "candidate.boundary_score", line
    )

    expected_values: dict[str, str | float | None] = {
        "flank_min": flank_min,
        "local_identity": local_identity,
        "aggregate_identity": aggregate_identity,
        "mapq_normalized": mapq_normalized,
        "target_n_fraction": target_n_fraction,
        "te_aligned_fraction": te_aligned_fraction,
        "insertion_fraction": insertion_fraction,
        "reciprocal_overlap": reciprocal,
        "boundary_score": boundary,
    }
    for field, expected in expected_values.items():
        require_printed_fraction_equal(feature[field], expected, field, line)

    if feature["family_relation"] != candidate["family_relation"]:
        raise ContractError(
            f"feature row {line} family_relation disagrees with candidate"
        )
    if feature["context_relation"] != candidate["context_relation"]:
        raise ContractError(
            f"feature row {line} context_relation disagrees with candidate"
        )

    mask = (
        FEATURE_BITS["te_alignment"] | FEATURE_BITS["insertion"]
        | FEATURE_BITS["reciprocal_overlap"] | FEATURE_BITS["boundary"]
    )
    if flank_min is not None:
        mask |= FEATURE_BITS["flank_min"]
    if local_identity is not None:
        mask |= FEATURE_BITS["local_identity"]
    elif aggregate_identity is not None:
        mask |= FEATURE_BITS["aggregate_identity"]
    if mapq_normalized is not None:
        mask |= FEATURE_BITS["mapq"]
    if target_n_fraction is not None:
        mask |= FEATURE_BITS["n_fraction"]
    if candidate["family_relation"] != "UNKNOWN":
        mask |= FEATURE_BITS["family"]
    if candidate["context_relation"] != "UNKNOWN":
        mask |= FEATURE_BITS["context"]
    expected_mask = f"0x{mask:08x}"
    if feature["observed_feature_mask"] != expected_mask:
        raise ContractError(
            f"feature row {line} observed_feature_mask disagrees with "
            f"reconstructed {expected_mask}"
        )

    missing = []
    for flag, name in (
        (FEATURE_BITS["flank_min"], "flank_min"),
        (FEATURE_BITS["local_identity"] | FEATURE_BITS["aggregate_identity"],
         "identity"),
        (FEATURE_BITS["mapq"], "mapq"),
        (FEATURE_BITS["n_fraction"], "n_fraction"),
        (FEATURE_BITS["family"], "family"),
        (FEATURE_BITS["context"], "context"),
    ):
        if not mask & flag:
            missing.append(name)
    expected_missing = ",".join(missing) if missing else UNKNOWN
    if feature["missing_features"] != expected_missing:
        raise ContractError(
            f"feature row {line} missing_features disagrees with reconstructed "
            f"{expected_missing!r}"
        )


def te_semantic_key(
    assembly_sha256: str, genome_id: str, contig: str, start: int, end: int,
    te_id: str,
) -> str:
    value = {
        "assembly_sha256": assembly_sha256, "genome_id": genome_id,
        "contig": contig, "start": start, "end": end, "te_id": te_id,
    }
    return "TEKsha256:" + sha256_text(canonical_json(value))


def atomic_write_tsv(path: Path, columns: Sequence[str], rows: Iterable[dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=False, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
            handle.flush()
            os.fsync(handle.fileno())
        return temporary
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_write_json(path: Path, value: Any) -> Path:
    path.parent.mkdir(parents=False, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        return temporary
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def export(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    dataset_id = validate_identifier(args.dataset_id, "--dataset-id")
    run_id = validate_identifier(args.run_id, "--run-id")
    prefix = Path(args.prefix)
    truth_path = Path(args.truth).resolve()
    run_path = Path(f"{prefix}.run.json").resolve()
    evidence_path = table_path(Path(f"{prefix}.evidence.tsv")).resolve()
    candidate_path = table_path(Path(f"{prefix}.candidates.tsv")).resolve()
    feature_path = table_path(Path(f"{prefix}.candidate_features.tsv")).resolve()
    output_prefix = Path(args.output)
    training_path = Path(f"{output_prefix}.training.tsv").resolve()
    unmatched_path = Path(f"{output_prefix}.unmatched_truth.tsv").resolve()
    dataset_path = Path(f"{output_prefix}.dataset.json").resolve()
    output_paths = {training_path, unmatched_path, dataset_path}
    input_paths = {run_path, evidence_path, candidate_path, feature_path, truth_path}
    if output_paths & input_paths:
        raise ContractError("an output path would overwrite an input")
    if len(output_paths) != 3:
        raise ContractError("output prefix does not produce three distinct output paths")
    if not training_path.parent.is_dir():
        raise ContractError(f"output directory does not exist: {training_path.parent}")

    run = read_json(run_path)
    config, software_version, schema_version, model_id = validate_run(run)
    evidence_rows = read_tsv(evidence_path, EVIDENCE_COLUMNS, "evidence TSV")
    candidate_rows = read_tsv(candidate_path, CANDIDATE_COLUMNS, "candidate TSV")
    feature_rows = read_tsv(feature_path, FEATURE_COLUMNS, "candidate feature TSV")
    truth_rows = read_tsv(truth_path, TRUTH_COLUMNS, "truth TSV")
    truth = validate_truth(truth_rows, dataset_id)
    validate_schema_versions(evidence_rows, schema_version, "evidence TSV")
    validate_schema_versions(candidate_rows, schema_version, "candidate TSV")
    validate_schema_versions(feature_rows, schema_version, "candidate feature TSV")

    evidence = unique_index(evidence_rows, "evidence_id", "evidence TSV")
    candidates = unique_index(candidate_rows, "candidate_id", "candidate TSV")
    features = unique_index(feature_rows, "candidate_id", "candidate feature TSV")
    if set(candidates) != set(features):
        missing = len(set(candidates) - set(features))
        extra = len(set(features) - set(candidates))
        raise ContractError(
            f"candidate/features are not one-to-one (missing={missing}, extra={extra})"
        )
    counts = require_object(run, "counts")
    expected_counts = {
        "evidence_observations": len(evidence_rows),
        "internal_candidates": len(candidate_rows),
        "candidates": len(candidate_rows),
        "candidate_feature_rows": len(feature_rows),
    }
    for field, observed in expected_counts.items():
        value = counts.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value != observed:
            raise ContractError(
                f"run count {field}={value!r} disagrees with exported rows {observed}"
            )

    frozen_inputs, frozen_by_path = load_frozen_input_inventory(run, run_path)
    genome_hashes, genome_inputs = load_genome_hashes(run, run_path)
    providers, provider_inputs = provider_contract(run, evidence_rows, run_path)
    for item in genome_inputs + provider_inputs:
        if frozen_by_path.get(item["path"]) != item["sha256"]:
            raise ContractError(
                f"run-specific input {item['path']!r} is not bound to input_files"
            )
    generator = {
        "candidate_generator_contract_version": "1.0.0",
        "software": "TEvoX", "software_version": software_version,
        "evidence_schema_version": schema_version,
        "coordinate_system": run["coordinate_system"],
        "config": config, "provider_set": providers,
    }
    generator_json = canonical_json(generator)
    generator_sha = sha256_text(generator_json)

    te_definitions: dict[tuple[str, str], tuple[str, int, int]] = {}
    for line, row in enumerate(evidence_rows, start=2):
        for field in ("evidence_id", "evidence_group_id", "provider", "origin",
                      "dependency", "query_genome_id", "target_genome_id",
                      "source_te_id", "source_contig"):
            clean(row[field], field, line)
        if row["query_genome_id"] not in genome_hashes or row["target_genome_id"] not in genome_hashes:
            raise ContractError(f"evidence row {line} references an unknown genome")
        start, end = validate_coordinates(row, "source", line)
        key = (row["query_genome_id"], row["source_te_id"])
        definition = (row["source_contig"], start, end)
        if key in te_definitions and te_definitions[key] != definition:
            raise ContractError(f"TE {key!r} has inconsistent semantic coordinates")
        te_definitions[key] = definition

    candidates_per_evidence: Counter[str] = Counter()
    for line, row in enumerate(candidate_rows, start=2):
        for field in ("candidate_id", "evidence_id", "target_genome_id",
                      "target_te_id", "target_contig"):
            clean(row[field], field, line)
        evidence_id = row["evidence_id"]
        if evidence_id not in evidence:
            raise ContractError(f"candidate row {line} references unknown evidence_id")
        candidates_per_evidence[evidence_id] += 1
        observation = evidence[evidence_id]
        if observation["target_genome_id"] != row["target_genome_id"]:
            raise ContractError(f"candidate row {line} target genome disagrees with evidence")
        if row["target_genome_id"] not in genome_hashes:
            raise ContractError(f"candidate row {line} references an unknown target genome")
        start, end = validate_coordinates(row, "target", line)
        key = (row["target_genome_id"], row["target_te_id"])
        definition = (row["target_contig"], start, end)
        if key in te_definitions and te_definitions[key] != definition:
            raise ContractError(f"TE {key!r} has inconsistent semantic coordinates")
        te_definitions[key] = definition
        finite_number(row["score"], "score", line)
        finite_number(row["reciprocal_overlap"], "reciprocal_overlap", line,
                      minimum=0.0, maximum=1.0)
        finite_number(row["boundary_score"], "boundary_score", line,
                      minimum=0.0, maximum=1.0)
        integer(row["breakpoint_distance"], "breakpoint_distance", line)
        integer(row["candidate_rank"], "candidate_rank", line, minimum=1)
        for field in ("context_compatible", "graph_retained", "eligible", "selected"):
            if row[field] not in BOOL_VALUES:
                raise ContractError(f"candidate row {line} field {field} is not boolean")
        if row["graph_retained"] != "true":
            raise ContractError(
                f"candidate row {line} is truncated despite max_graph_candidates=0"
            )

    for line, row in enumerate(evidence_rows, start=2):
        evidence_id = row["evidence_id"]
        actual = candidates_per_evidence[evidence_id]
        for field in (
            "nearby_candidate_count", "retained_candidate_count",
            "graph_candidate_count",
        ):
            reported = integer(row[field], field, line)
            if reported != actual:
                raise ContractError(
                    f"evidence row {line} reports {field}={reported}, but "
                    f"complete candidate output contains {actual} row(s)"
                )

    model_ids: set[str] = set()
    candidate_window = config["candidate_window"]
    assert isinstance(candidate_window, int)
    for line, row in enumerate(feature_rows, start=2):
        validate_features(row, line)
        model_ids.add(clean(row["model_id"], "model_id", line))
        candidate = candidates[row["candidate_id"]]
        if row["evidence_id"] != candidate["evidence_id"]:
            raise ContractError(f"feature row {line} evidence_id disagrees with candidate")
        if row["eligible"] != candidate["eligible"]:
            raise ContractError(f"feature row {line} eligible disagrees with candidate")
        observation = evidence[candidate["evidence_id"]]
        validate_raw_feature_sources(
            row, candidate, observation, candidate_window, line
        )
    if model_ids != {model_id}:
        raise ContractError(
            f"feature model IDs {sorted(model_ids)!r} disagree with run model {model_id!r}"
        )

    semantic_candidates: Counter[
        tuple[tuple[str, str], tuple[str, str]]
    ] = Counter()
    exported: list[dict[str, Any]] = []
    matched_truth_ids: set[str] = set()
    for candidate_id in sorted(candidates):
        candidate = candidates[candidate_id]
        feature = features[candidate_id]
        observation = evidence[candidate["evidence_id"]]
        source_genome = observation["query_genome_id"]
        source_te = observation["source_te_id"]
        target_genome = candidate["target_genome_id"]
        target_te = candidate["target_te_id"]
        pair = canonical_pair(source_genome, source_te, target_genome, target_te)
        semantic_candidates[pair] += 1
        truth_row = truth.get(pair)
        if truth_row is None:
            truth_status = "UNLABELLED"
            truth_values = {field: UNKNOWN for field in TRUTH_COLUMNS}
        else:
            truth_status = (
                "ASSESSED_UNKNOWN" if truth_row["label"] == "UNKNOWN"
                else "ASSESSED_LABEL"
            )
            truth_values = orient_truth(
                truth_row, source_genome, source_te, target_genome, target_te
            )
            matched_truth_ids.add(truth_row["truth_record_id"])
        source_start, source_end = validate_coordinates(observation, "source", 0)
        target_start, target_end = validate_coordinates(candidate, "target", 0)
        source_hash = genome_hashes[source_genome]
        target_hash = genome_hashes[target_genome]
        row: dict[str, Any] = {
            "export_schema_version": EXPORT_SCHEMA_VERSION,
            "export_status": EXPORT_STATUS, "model_fit_status": MODEL_FIT_STATUS,
            "calibration_status": CALIBRATION_STATUS,
            "benchmark_schema_version": BENCHMARK_SCHEMA_VERSION,
            "dataset_id": dataset_id, "run_id": run_id,
            "truth_status": truth_status,
            "truth_record_id": truth_values["truth_record_id"],
            "label": truth_values["label"],
            "truth_confidence": truth_values["truth_confidence"],
            "ancestral_event_id": truth_values["ancestral_event_id"],
            "source_truth_locus_id": truth_values["source_truth_locus_id"],
            "target_truth_locus_id": truth_values["target_truth_locus_id"],
            "validation_method": truth_values["validation_method"],
            "validation_source": truth_values["validation_source"],
            "validation_batch_id": truth_values["validation_batch_id"],
            "curation_blinded": truth_values["curation_blinded"],
            "taxon_a": truth_values["taxon_a"],
            "taxon_b": truth_values["taxon_b"],
            "clade_holdout_id": truth_values["clade_holdout_id"],
            "candidate_generator_sha256": generator_sha,
            "inclusion_probability": "1", "candidate_id": candidate_id,
            "evidence_id": candidate["evidence_id"],
            "evidence_group_id": observation["evidence_group_id"],
            "provider": observation["provider"], "origin": observation["origin"],
            "dependency": observation["dependency"],
            "source_genome_id": source_genome, "source_te_id": source_te,
            "source_te_semantic_key": te_semantic_key(
                source_hash["assembly_sha256"], source_genome,
                observation["source_contig"], source_start, source_end, source_te,
            ),
            "source_assembly_sha256": source_hash["assembly_sha256"],
            "source_annotation_sha256": source_hash["annotation_sha256"],
            "source_contig": observation["source_contig"],
            "source_start": source_start, "source_end": source_end,
            "target_genome_id": target_genome, "target_te_id": target_te,
            "target_te_semantic_key": te_semantic_key(
                target_hash["assembly_sha256"], target_genome,
                candidate["target_contig"], target_start, target_end, target_te,
            ),
            "target_assembly_sha256": target_hash["assembly_sha256"],
            "target_annotation_sha256": target_hash["annotation_sha256"],
            "target_contig": candidate["target_contig"],
            "target_start": target_start, "target_end": target_end,
        }
        row.update({field: feature[field] for field in RAW_FEATURE_COLUMNS})
        exported.append(row)

    unmatched = [
        {**row, "unmatched_reason": "NO_GENERATED_CANDIDATE"}
        for row in sorted(truth_rows, key=lambda item: item["truth_record_id"])
        if row["truth_record_id"] not in matched_truth_ids
    ]
    label_counts = Counter(row["label"] for row in exported)
    truth_status_counts = Counter(row["truth_status"] for row in exported)
    truth_label_counts = Counter(row["label"] for row in truth_rows)
    base_inputs = [
        {"role": "run_metadata", "path": str(run_path), "sha256": sha256_file(run_path)},
        {"role": "evidence_output", "path": str(evidence_path), "sha256": sha256_file(evidence_path)},
        {"role": "candidate_output", "path": str(candidate_path), "sha256": sha256_file(candidate_path)},
        {"role": "candidate_feature_output", "path": str(feature_path), "sha256": sha256_file(feature_path)},
        {"role": "candidate_truth", "path": str(truth_path), "sha256": sha256_file(truth_path)},
    ]
    input_manifest = base_inputs + frozen_inputs
    input_manifest.sort(
        key=lambda item: (
            item["role"], item.get("genome_id", ""), item.get("provider", ""),
            item["path"],
        )
    )
    input_manifest_sha = sha256_text(canonical_json(input_manifest))
    dataset = {
        "export_schema_version": EXPORT_SCHEMA_VERSION,
        "export_status": EXPORT_STATUS, "model_fit_status": MODEL_FIT_STATUS,
        "calibration_status": CALIBRATION_STATUS,
        "score_semantics": "NO_SCORES_EXPORTED_NO_PROBABILITY_CLAIM",
        "benchmark_schema_version": BENCHMARK_SCHEMA_VERSION,
        "dataset_id": dataset_id, "run_id": run_id,
        "coordinate_system": "0-based-half-open",
        "join_contract": {
            "truth_join": "UNORDERED_SEMANTIC_PAIR_OF_GENOME_ID_AND_TE_ID",
            "forbidden_truth_join_fields": [
                "run_id", "candidate_id", "evidence_id", "evidence_group_id",
                "decision_id", "edge_id", "locus_id", "solver_component_id",
            ],
            "candidate_id_semantics": "RUN_LOCAL_FOREIGN_KEY_ONLY",
            "evidence_group_id_semantics": "RUN_LOCAL_LEAKAGE_BINDING_ONLY",
            "te_semantic_key": (
                "TEKsha256:<hex SHA256(canonical_json({assembly_sha256,genome_id,"
                "contig,start,end,te_id}))>"
            ),
        },
        "candidate_generator": generator,
        "candidate_generator_canonical_json": generator_json,
        "candidate_generator_sha256": generator_sha,
        "sampling": {
            "downsampled": False, "inclusion_probability": 1,
            "max_candidates": 0, "max_graph_candidates": 0,
            "scope": "ALL_GENERATED_CANDIDATE_VIEWS",
            "truth_ascertainment_probability": "UNKNOWN_NOT_ESTIMATED",
        },
        "features": {
            "exported": list(RAW_FEATURE_COLUMNS),
            "excluded_post_decision_or_builtin_score_fields": [
                "membership_logit", "membership_score", "membership_entropy",
                "membership_prediction_set", "out_of_domain", "eligible",
                "score", "candidate_rank", "graph_retained", "selected",
                "decision_code", "decision_id", "solver_component_id", "locus_id",
            ],
        },
        "counts": {
            "truth_records": len(truth_rows),
            "truth_label_counts": dict(sorted(truth_label_counts.items())),
            "candidate_rows": len(candidate_rows),
            "semantic_candidate_pairs": len(semantic_candidates),
            "candidate_pairs_without_truth": sum(
                1 for pair in semantic_candidates if pair not in truth
            ),
            "candidate_rows_without_truth": sum(
                count for pair, count in semantic_candidates.items() if pair not in truth
            ),
            "matched_truth_records": len(matched_truth_ids),
            "unmatched_truth_records": len(unmatched),
            "exported_rows": len(exported),
            "exported_label_counts": dict(sorted(label_counts.items())),
            "exported_truth_status_counts": dict(
                sorted(truth_status_counts.items())
            ),
        },
        "input_files": input_manifest,
        "input_manifest_sha256": input_manifest_sha,
        "genomes": [
            {
                "genome_id": genome_id,
                "assembly_sha256": values["assembly_sha256"],
                "annotation_sha256": values["annotation_sha256"],
                "fasta_path": values["fasta_path"],
                "annotation_path": values["annotation_path"],
            }
            for genome_id, values in sorted(genome_hashes.items())
        ],
        "outputs": {
            "training_tsv": str(training_path),
            "unmatched_truth_tsv": str(unmatched_path),
            "dataset_json": str(dataset_path),
        },
    }

    temporaries: list[Path] = []
    try:
        training_temporary = atomic_write_tsv(
            training_path, TRAINING_COLUMNS, exported
        )
        temporaries.append(training_temporary)
        unmatched_temporary = atomic_write_tsv(
            unmatched_path, UNMATCHED_COLUMNS, unmatched
        )
        temporaries.append(unmatched_temporary)
        training_sha = sha256_file(training_temporary)
        unmatched_sha = sha256_file(unmatched_temporary)
        training_size = training_temporary.stat().st_size
        unmatched_size = unmatched_temporary.stat().st_size
        transaction = {
            "protocol": "DATASET_JSON_LAST_V1",
            "export_schema_version": EXPORT_SCHEMA_VERSION,
            "input_manifest_sha256": input_manifest_sha,
            "candidate_generator_sha256": generator_sha,
            "training_tsv_sha256": training_sha,
            "unmatched_truth_tsv_sha256": unmatched_sha,
        }
        dataset["outputs"].update({
            "training_tsv_sha256": training_sha,
            "training_tsv_bytes": training_size,
            "unmatched_truth_tsv_sha256": unmatched_sha,
            "unmatched_truth_tsv_bytes": unmatched_size,
        })
        dataset["commit_marker"] = {
            **transaction,
            "status": "COMMITTED_IF_DATASET_JSON_PRESENT",
            "transaction_sha256": sha256_text(canonical_json(transaction)),
        }
        dataset_temporary = atomic_write_json(dataset_path, dataset)
        temporaries.append(dataset_temporary)

        # A stale marker must never describe a partially replaced pair of TSVs.
        # Removing it first makes absence of dataset.json the explicit
        # interrupted-transaction state.  The new marker is committed last.
        dataset_path.unlink(missing_ok=True)
        fsync_directory(dataset_path.parent)
        os.replace(training_temporary, training_path)
        fsync_directory(training_path.parent)
        os.replace(unmatched_temporary, unmatched_path)
        fsync_directory(unmatched_path.parent)
        os.replace(dataset_temporary, dataset_path)
        fsync_directory(dataset_path.parent)
    finally:
        for temporary in temporaries:
            temporary.unlink(missing_ok=True)
    return training_path, unmatched_path, dataset_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "export complete pre-decision candidate features against independent "
            "semantic truth; no model is fitted or calibrated"
        )
    )
    parser.add_argument("--prefix", required=True, help="TEvoX output prefix")
    parser.add_argument("--truth", required=True, help="candidate truth TSV")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output", required=True, help="output prefix")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        paths = export(args)
    except (ContractError, OSError) as error:
        print(f"tevox-export-training: error: {error}", file=sys.stderr)
        return 2
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
