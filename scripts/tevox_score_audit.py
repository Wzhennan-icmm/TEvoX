#!/usr/bin/env python3
"""Audit uncalibrated TEvoX membership scores against independent truth."""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import os
import random
import tempfile
from pathlib import Path
from typing import Iterable, Sequence


Record = tuple[float, int, str]


def resolved_path(path: Path, label: str) -> Path:
    try:
        return path.resolve(strict=False)
    except (OSError, RuntimeError) as exc:
        raise SystemExit(f"cannot resolve {label} path {path}: {exc}") from exc


def require_outputs_disjoint_from_inputs(
    output_paths: Sequence[tuple[str, Path]],
    input_paths: Sequence[tuple[str, Path]],
) -> None:
    resolved_outputs = [
        (label, resolved_path(path, label)) for label, path in output_paths
    ]
    for output_label, output_path in resolved_outputs:
        for input_label, input_path in input_paths:
            input_resolved = resolved_path(input_path, input_label)
            if output_path == input_resolved:
                raise SystemExit(
                    f"resolved {output_label} output path {output_path} would "
                    f"overwrite the {input_label} input"
                )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            while block := handle.read(1024 * 1024):
                digest.update(block)
    except OSError as exc:
        raise SystemExit(f"cannot hash output {path}: {exc}") from exc
    return digest.hexdigest()


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def read_rows(path: Path, table_name: str) -> tuple[list[str], list[dict[str, str]]]:
    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as exc:
        raise SystemExit(f"cannot open {table_name} TSV {path}: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields:
            raise SystemExit(f"{table_name} TSV has no header")
        if any(field is None or not field.strip() for field in fields):
            raise SystemExit(f"{table_name} TSV contains an empty header field")
        if len(fields) != len(set(fields)):
            raise SystemExit(f"{table_name} TSV contains duplicate header fields")
        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise SystemExit(
                    f"{table_name} TSV row {line_number} has the wrong number of columns"
                )
            if not any(value.strip() for value in row.values()):
                raise SystemExit(f"{table_name} TSV row {line_number} is empty")
            rows.append(row)
    if not rows:
        raise SystemExit(f"{table_name} TSV contains no data rows")
    return fields, rows


def require_columns(fields: list[str], required: set[str], table_name: str) -> None:
    missing = sorted(required - set(fields))
    if missing:
        raise SystemExit(
            f"{table_name} TSV is missing required column(s): {', '.join(missing)}"
        )


def finite_score(text: str) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise ValueError(f"invalid membership_score {text!r}") from exc
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"invalid membership_score {text!r}")
    return value


def binary_label(text: str) -> int:
    if text not in {"0", "1"}:
        raise ValueError(f"truth label must be 0 or 1, found {text!r}")
    return int(text)


def auc(records: list[Record]) -> float | None:
    positives = sum(label for _, label, _ in records)
    negatives = len(records) - positives
    if positives == 0 or negatives == 0:
        return None
    ordered = sorted(records, key=lambda item: item[0])
    positive_rank_sum = 0.0
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][0] == ordered[start][0]:
            end += 1
        average_rank = ((start + 1) + end) / 2.0
        positive_rank_sum += average_rank * sum(
            ordered[index][1] for index in range(start, end)
        )
        start = end
    return (
        positive_rank_sum - positives * (positives + 1) / 2.0
    ) / (positives * negatives)


def auprc(records: list[Record]) -> float | None:
    """Return tied-threshold average precision (step-integrated PR area)."""
    positives = sum(label for _, label, _ in records)
    if positives == 0:
        return None
    ordered = sorted(records, key=lambda item: item[0], reverse=True)
    true_positives = 0
    false_positives = 0
    previous_recall = 0.0
    area = 0.0
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][0] == ordered[start][0]:
            end += 1
        tied_positives = sum(ordered[index][1] for index in range(start, end))
        true_positives += tied_positives
        false_positives += (end - start) - tied_positives
        recall = true_positives / positives
        precision = true_positives / (true_positives + false_positives)
        area += (recall - previous_recall) * precision
        previous_recall = recall
        start = end
    return area


def _sigmoid(value: float) -> float:
    if value >= 0.0:
        factor = math.exp(-value)
        return 1.0 / (1.0 + factor)
    factor = math.exp(value)
    return factor / (1.0 + factor)


def _logit(score: float, epsilon: float = 1e-6) -> float:
    clipped = min(1.0 - epsilon, max(epsilon, score))
    return math.log(clipped / (1.0 - clipped))


def calibration_intercept(records: list[Record]) -> float | None:
    """Fit the diagnostic intercept with membership-score logit as an offset."""
    positives = sum(label for _, label, _ in records)
    if positives == 0 or positives == len(records):
        return None
    offsets = [_logit(score) for score, _, _ in records]
    target = float(positives)
    low = -60.0
    high = 60.0
    for _ in range(200):
        middle = (low + high) / 2.0
        expected = sum(_sigmoid(offset + middle) for offset in offsets)
        if expected < target:
            low = middle
        else:
            high = middle
    result = (low + high) / 2.0
    return result if math.isfinite(result) else None


def calibration_slope(records: list[Record]) -> float | None:
    """Fit an unpenalized diagnostic logistic slope, or return null if unstable."""
    if len(records) < 3:
        return None
    labels = [label for _, label, _ in records]
    positives = sum(labels)
    if positives == 0 or positives == len(labels):
        return None
    logits = [_logit(score) for score, _, _ in records]
    mean_logit = sum(logits) / len(logits)
    variance = sum((value - mean_logit) ** 2 for value in logits) / len(logits)
    if variance <= 1e-14:
        return None
    scale = math.sqrt(variance)
    predictors = [(value - mean_logit) / scale for value in logits]

    intercept = math.log(positives / (len(labels) - positives))
    slope = 0.0

    def log_likelihood(candidate_intercept: float, candidate_slope: float) -> float:
        total = 0.0
        for predictor, label in zip(predictors, labels):
            linear = candidate_intercept + candidate_slope * predictor
            total += label * linear - max(linear, 0.0) - math.log1p(
                math.exp(-abs(linear))
            )
        return total

    converged = False
    current_likelihood = log_likelihood(intercept, slope)
    for _ in range(100):
        gradient_intercept = 0.0
        gradient_slope = 0.0
        information_00 = 0.0
        information_01 = 0.0
        information_11 = 0.0
        for predictor, label in zip(predictors, labels):
            fitted = _sigmoid(intercept + slope * predictor)
            residual = label - fitted
            weight = fitted * (1.0 - fitted)
            gradient_intercept += residual
            gradient_slope += residual * predictor
            information_00 += weight
            information_01 += weight * predictor
            information_11 += weight * predictor * predictor
        determinant = information_00 * information_11 - information_01**2
        if determinant <= 1e-14 or not math.isfinite(determinant):
            return None
        delta_intercept = (
            information_11 * gradient_intercept
            - information_01 * gradient_slope
        ) / determinant
        delta_slope = (
            information_00 * gradient_slope
            - information_01 * gradient_intercept
        ) / determinant
        step = 1.0
        accepted = False
        while step >= 2.0**-20:
            proposed_intercept = intercept + step * delta_intercept
            proposed_slope = slope + step * delta_slope
            proposed_likelihood = log_likelihood(
                proposed_intercept, proposed_slope
            )
            if proposed_likelihood >= current_likelihood - 1e-12:
                intercept = proposed_intercept
                slope = proposed_slope
                current_likelihood = proposed_likelihood
                accepted = True
                break
            step /= 2.0
        if not accepted:
            return None
        if max(abs(step * delta_intercept), abs(step * delta_slope)) < 1e-9:
            converged = True
            break
        if max(abs(intercept), abs(slope)) > 50.0:
            return None
    if not converged:
        return None
    result = slope / scale
    if not math.isfinite(result) or abs(result) > 1e6:
        return None
    return result


def _bin_summary(members: list[tuple[float, int]], index: int) -> dict[str, object]:
    if members:
        mean_score = sum(score for score, _ in members) / len(members)
        observed = sum(label for _, label in members) / len(members)
        gap = abs(mean_score - observed)
        score_min: float | None = min(score for score, _ in members)
        score_max: float | None = max(score for score, _ in members)
    else:
        mean_score = None
        observed = None
        gap = None
        score_min = None
        score_max = None
    return {
        "bin": index,
        "score_min": score_min,
        "score_max": score_max,
        "count": len(members),
        "mean_score": mean_score,
        "observed_fraction": observed,
        "absolute_gap": gap,
    }


def equal_width_bins(records: list[Record], bins: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for index in range(bins):
        low = index / bins
        high = (index + 1) / bins
        members = [
            (score, label)
            for score, label, _ in records
            if low <= score < high or (index == bins - 1 and score == 1.0)
        ]
        rows.append(
            {
                **_bin_summary(members, index + 1),
                "lower": low,
                "upper": high,
            }
        )
    return rows


def equal_mass_bins(records: list[Record], bins: int) -> list[dict[str, object]]:
    """Make approximately equal-count bins without splitting tied scores."""
    scores = sorted(score for score, _, _ in records)
    maximum = scores[-1]
    boundaries: list[float] = []
    for index in range(1, bins):
        rank = math.ceil(index * len(scores) / bins) - 1
        boundary = scores[rank]
        if boundary < maximum and (not boundaries or boundary > boundaries[-1]):
            boundaries.append(boundary)
    members_by_bin: list[list[tuple[float, int]]] = [
        [] for _ in range(len(boundaries) + 1)
    ]
    for score, label, _ in records:
        members_by_bin[bisect.bisect_left(boundaries, score)].append((score, label))
    rows: list[dict[str, object]] = []
    for index, members in enumerate(members_by_bin, start=1):
        rows.append(
            {
                **_bin_summary(members, index),
                "lower": min(score for score, _ in members),
                "upper": max(score for score, _ in members),
            }
        )
    return rows


def expected_calibration_error(
    rows: list[dict[str, object]], total: int
) -> float:
    return sum(
        int(row["count"]) / total * float(row["absolute_gap"])
        for row in rows
        if row["count"]
    )


def summarize(
    records: list[Record], bins: int
) -> tuple[
    dict[str, object], list[dict[str, object]], list[dict[str, object]]
]:
    if not records:
        raise ValueError("cannot audit an empty truth set")
    epsilon = 1e-15
    brier = sum((score - label) ** 2 for score, label, _ in records) / len(records)
    log_loss = -sum(
        label * math.log(min(1.0 - epsilon, max(epsilon, score)))
        + (1 - label)
        * math.log(min(1.0 - epsilon, max(epsilon, 1.0 - score)))
        for score, label, _ in records
    ) / len(records)
    width_rows = equal_width_bins(records, bins)
    mass_rows = equal_mass_bins(records, bins)
    mean_score = sum(score for score, _, _ in records) / len(records)
    prevalence = sum(label for _, label, _ in records) / len(records)
    width_ece = expected_calibration_error(width_rows, len(records))
    metrics: dict[str, object] = {
        "n": len(records),
        "positives": sum(label for _, label, _ in records),
        "negatives": sum(1 - label for _, label, _ in records),
        "mean_membership_score": mean_score,
        "observed_positive_fraction": prevalence,
        "calibration_in_the_large": prevalence - mean_score,
        "calibration_intercept": calibration_intercept(records),
        "calibration_slope": calibration_slope(records),
        "brier_score": brier,
        "log_loss": log_loss,
        # Preserve the v1 field as equal-width ECE for downstream compatibility.
        "ece": width_ece,
        "ece_equal_width": width_ece,
        "ece_equal_mass": expected_calibration_error(mass_rows, len(records)),
        "auroc": auc(records),
        "auprc": auprc(records),
    }
    return metrics, width_rows, mass_rows


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def grouped_bootstrap(
    records: list[Record], bins: int, replicates: int, seed: int
) -> dict[str, object]:
    groups: dict[str, list[Record]] = {}
    for record in records:
        groups.setdefault(record[2], []).append(record)
    group_ids = sorted(groups)
    metric_names = (
        "brier_score",
        "log_loss",
        "ece_equal_width",
        "ece_equal_mass",
        "auroc",
        "auprc",
        "calibration_in_the_large",
        "calibration_intercept",
        "calibration_slope",
    )
    samples: dict[str, list[float]] = {name: [] for name in metric_names}
    generator = random.Random(seed)
    for _ in range(replicates):
        sampled: list[Record] = []
        for _ in group_ids:
            sampled.extend(groups[generator.choice(group_ids)])
        metrics, _, _ = summarize(sampled, bins)
        for name in metric_names:
            value = metrics[name]
            if value is not None and math.isfinite(float(value)):
                samples[name].append(float(value))
    intervals: dict[str, dict[str, object]] = {}
    for name in metric_names:
        values = samples[name]
        intervals[name] = {
            "lower": percentile(values, 0.025) if values else None,
            "upper": percentile(values, 0.975) if values else None,
            "valid_replicates": len(values),
        }
    return {
        "method": "PERCENTILE_GROUP_RESAMPLING",
        "confidence_level": 0.95,
        "replicates_requested": replicates,
        "seed": seed,
        "resampling_groups": len(group_ids),
        "intervals": intervals,
    }


def write_bins(
    path: Path,
    scopes: Iterable[tuple[str, str, list[dict[str, object]]]],
) -> None:
    fields = [
        "scope",
        "binning_method",
        "bin",
        "lower",
        "upper",
        "score_min",
        "score_max",
        "count",
        "mean_score",
        "observed_fraction",
        "absolute_gap",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for scope, method, rows in scopes:
            for row in rows:
                output = {"scope": scope, "binning_method": method, **row}
                for field in (
                    "lower",
                    "upper",
                    "score_min",
                    "score_max",
                    "mean_score",
                    "observed_fraction",
                    "absolute_gap",
                ):
                    if output[field] is None:
                        output[field] = "."
                writer.writerow(output)
        handle.flush()
        os.fsync(handle.fileno())


def write_output_bundle(
    metrics_path: Path,
    bins_path: Path,
    report: dict[str, object],
    bin_scopes: Iterable[tuple[str, str, list[dict[str, object]]]],
) -> None:
    """Publish the TSV first and its JSON manifest last.

    Both files are fully written and fsynced under temporary names in the
    destination directory.  A reader can use the manifest's TSV digest to
    reject an interrupted mixed-version pair.
    """
    temporary_bins: Path | None = None
    temporary_metrics: Path | None = None
    try:
        metrics_path.parent.mkdir(parents=True, exist_ok=True)
        bins_path.parent.mkdir(parents=True, exist_ok=True)

        bins_fd, bins_name = tempfile.mkstemp(
            dir=bins_path.parent,
            prefix=f".{bins_path.name}.",
            suffix=".tmp",
        )
        os.close(bins_fd)
        temporary_bins = Path(bins_name)
        write_bins(temporary_bins, bin_scopes)
        bins_sha256 = sha256_file(temporary_bins)

        report["artifacts"] = {
            "calibration_tsv": {
                "path": str(bins_path),
                "sha256": bins_sha256,
            }
        }
        report["calibration_tsv_sha256"] = bins_sha256

        metrics_fd, metrics_name = tempfile.mkstemp(
            dir=metrics_path.parent,
            prefix=f".{metrics_path.name}.",
            suffix=".tmp",
        )
        temporary_metrics = Path(metrics_name)
        with os.fdopen(metrics_fd, "w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())

        os.replace(temporary_bins, bins_path)
        temporary_bins = None
        fsync_directory(bins_path.parent)
        # The JSON is the commit marker for the pair and must be replaced last.
        os.replace(temporary_metrics, metrics_path)
        temporary_metrics = None
        fsync_directory(metrics_path.parent)
    except OSError as exc:
        raise SystemExit(
            f"cannot atomically write {metrics_path} and {bins_path}: {exc}"
        ) from exc
    finally:
        for temporary_path in (temporary_bins, temporary_metrics):
            if temporary_path is not None:
                try:
                    temporary_path.unlink()
                except FileNotFoundError:
                    pass


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate TEvoX membership scores without converting them into "
            "calibrated probabilities"
        )
    )
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=0,
        help="group-resampling bootstrap replicates (default: disabled)",
    )
    parser.add_argument(
        "--seed", type=int, default=1, help="deterministic bootstrap seed"
    )
    args = parser.parse_args()
    if not 2 <= args.bins <= 100:
        parser.error("--bins must be between 2 and 100")
    if not 0 <= args.bootstrap_replicates <= 100000:
        parser.error("--bootstrap-replicates must be between 0 and 100000")

    metrics_path = Path(f"{args.output}.metrics.json")
    bins_path = Path(f"{args.output}.calibration.tsv")
    require_outputs_disjoint_from_inputs(
        (("metrics JSON", metrics_path), ("calibration TSV", bins_path)),
        (("features", args.features), ("truth", args.truth)),
    )

    feature_fields, features = read_rows(args.features, "features")
    truth_fields, truth = read_rows(args.truth, "truth")
    require_columns(
        feature_fields,
        {"candidate_id", "membership_score", "model_id", "calibration_status"},
        "features",
    )
    require_columns(
        truth_fields, {"candidate_id", "label", "group_id"}, "truth"
    )

    by_id: dict[str, dict[str, str]] = {}
    score_by_id: dict[str, float] = {}
    model_ids: set[str] = set()
    calibration_states: set[str] = set()
    try:
        for row in features:
            candidate_id = row["candidate_id"].strip()
            if not candidate_id or candidate_id == "." or candidate_id in by_id:
                raise SystemExit(
                    f"missing or duplicate feature candidate_id {candidate_id!r}"
                )
            model_id = row["model_id"].strip()
            if not model_id or model_id == ".":
                raise SystemExit(
                    f"feature candidate {candidate_id!r} has no model_id"
                )
            calibration_status = row["calibration_status"].strip()
            by_id[candidate_id] = row
            score_by_id[candidate_id] = finite_score(row["membership_score"].strip())
            model_ids.add(model_id)
            calibration_states.add(calibration_status)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if calibration_states != {"UNCALIBRATED"}:
        raise SystemExit(
            "score audit expects explicitly UNCALIBRATED feature rows"
        )
    if len(model_ids) != 1:
        raise SystemExit(
            "score audit requires exactly one model_id; found "
            + ", ".join(sorted(model_ids))
        )

    records: list[Record] = []
    seen_truth: set[str] = set()
    all_truth_groups: set[str] = set()
    group_coverage: dict[str, dict[str, int]] = {}
    missing_feature_rows = 0
    try:
        for row in truth:
            candidate_id = row["candidate_id"].strip()
            if not candidate_id or candidate_id == ".":
                raise SystemExit("truth TSV contains an empty candidate_id")
            if candidate_id in seen_truth:
                raise SystemExit(f"duplicate truth candidate_id {candidate_id!r}")
            seen_truth.add(candidate_id)
            label = binary_label(row["label"].strip())
            group = row["group_id"].strip()
            if not group or group == ".":
                raise SystemExit(f"truth candidate {candidate_id!r} has no group_id")
            all_truth_groups.add(group)
            coverage = group_coverage.setdefault(group, {"truth_rows": 0, "matched_rows": 0})
            coverage["truth_rows"] += 1
            if candidate_id not in by_id:
                missing_feature_rows += 1
                continue
            coverage["matched_rows"] += 1
            records.append((score_by_id[candidate_id], label, group))
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if not records:
        raise SystemExit("truth TSV has no candidate IDs present in the features TSV")

    groups = sorted({group for _, _, group in records})
    overall, overall_width_bins, overall_mass_bins = summarize(records, args.bins)
    group_metrics: dict[str, dict[str, object]] = {}
    bin_scopes: list[tuple[str, str, list[dict[str, object]]]] = [
        ("overall", "EQUAL_WIDTH", overall_width_bins),
        ("overall", "EQUAL_MASS_TIES_PRESERVED", overall_mass_bins),
    ]
    for group in groups:
        subset = [record for record in records if record[2] == group]
        metrics, group_width_bins, group_mass_bins = summarize(subset, args.bins)
        group_metrics[group] = metrics
        bin_scopes.extend(
            [
                (f"group:{group}", "EQUAL_WIDTH", group_width_bins),
                (
                    f"group:{group}",
                    "EQUAL_MASS_TIES_PRESERVED",
                    group_mass_bins,
                ),
            ]
        )

    for coverage in group_coverage.values():
        coverage["missing_feature_rows"] = (
            coverage["truth_rows"] - coverage["matched_rows"]
        )
    matched_rows = len(records)
    feature_rows_with_truth = len(seen_truth & set(by_id))
    coverage_report = {
        "feature_rows": len(features),
        "truth_rows": len(truth),
        "matched_truth_feature_rows": matched_rows,
        "truth_rows_missing_features": missing_feature_rows,
        "feature_rows_with_truth": feature_rows_with_truth,
        "feature_rows_without_truth": len(features) - feature_rows_with_truth,
        "truth_feature_coverage": matched_rows / len(truth),
        "feature_truth_coverage": feature_rows_with_truth / len(features),
        "by_truth_group": {key: group_coverage[key] for key in sorted(group_coverage)},
    }

    bootstrap_report = grouped_bootstrap(
        records, args.bins, args.bootstrap_replicates, args.seed
    )
    report = {
        "audit_version": "1.1.0",
        "model_ids": sorted(model_ids),
        "input_calibration_status": "UNCALIBRATED",
        "report_status": "EVALUATION_ONLY_NOT_A_CALIBRATED_MODEL",
        "score_semantics": "UNCALIBRATED_MEMBERSHIP_SCORE",
        "truth_groups": sorted(all_truth_groups),
        "evaluated_truth_groups": groups,
        "coverage": coverage_report,
        "overall": overall,
        "group_wise_evaluation": group_metrics,
        "group_bootstrap_95_ci": bootstrap_report,
        "calibration_diagnostics": {
            "legacy_ece_binning": "EQUAL_WIDTH",
            "additional_ece_binning": "EQUAL_MASS_TIES_PRESERVED",
            "auprc_method": "TIED_THRESHOLD_AVERAGE_PRECISION",
            "calibration_in_the_large_definition": (
                "observed_positive_fraction-minus-mean_membership_score"
            ),
            "calibration_intercept_definition": (
                "diagnostic_logistic_intercept_with_membership_score_logit_offset"
            ),
            "calibration_slope_definition": (
                "diagnostic_unpenalized_logistic_slope_or_null_when_unstable"
            ),
            "logit_clipping_for_intercept_and_slope": "[1e-6,1-1e-6]",
        },
        "features": str(args.features),
        "truth": str(args.truth),
    }
    write_output_bundle(metrics_path, bins_path, report, bin_scopes)
    print(f"Wrote {metrics_path} and {bins_path}")


if __name__ == "__main__":
    main()
