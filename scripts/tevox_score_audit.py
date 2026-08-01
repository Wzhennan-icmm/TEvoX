#!/usr/bin/env python3
"""Audit uncalibrated TEvoX membership scores against independent truth."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite_score(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"invalid membership_score {text!r}")
    return value


def binary_label(text: str) -> int:
    if text not in {"0", "1"}:
        raise ValueError(f"truth label must be 0 or 1, found {text!r}")
    return int(text)


def auc(records: list[tuple[float, int, str]]) -> float | None:
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


def summarize(
    records: list[tuple[float, int, str]], bins: int
) -> tuple[dict[str, object], list[dict[str, object]]]:
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
    rows: list[dict[str, object]] = []
    ece = 0.0
    for index in range(bins):
        low = index / bins
        high = (index + 1) / bins
        members = [
            (score, label)
            for score, label, _ in records
            if low <= score < high or (index == bins - 1 and score == 1.0)
        ]
        if members:
            mean_score = sum(score for score, _ in members) / len(members)
            observed = sum(label for _, label in members) / len(members)
            gap = abs(mean_score - observed)
            ece += len(members) / len(records) * gap
        else:
            mean_score = None
            observed = None
            gap = None
        rows.append(
            {
                "bin": index + 1,
                "lower": low,
                "upper": high,
                "count": len(members),
                "mean_score": mean_score,
                "observed_fraction": observed,
                "absolute_gap": gap,
            }
        )
    metrics: dict[str, object] = {
        "n": len(records),
        "positives": sum(label for _, label, _ in records),
        "negatives": sum(1 - label for _, label, _ in records),
        "brier_score": brier,
        "log_loss": log_loss,
        "ece": ece,
        "auroc": auc(records),
    }
    return metrics, rows


def write_bins(
    path: Path,
    scopes: Iterable[tuple[str, list[dict[str, object]]]],
) -> None:
    fields = [
        "scope",
        "bin",
        "lower",
        "upper",
        "count",
        "mean_score",
        "observed_fraction",
        "absolute_gap",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for scope, rows in scopes:
            for row in rows:
                output = {"scope": scope, **row}
                for field in ("mean_score", "observed_fraction", "absolute_gap"):
                    if output[field] is None:
                        output[field] = "."
                writer.writerow(output)


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
    args = parser.parse_args()
    if not 2 <= args.bins <= 100:
        parser.error("--bins must be between 2 and 100")

    features = read_rows(args.features)
    truth = read_rows(args.truth)
    required_truth = {"candidate_id", "label", "group_id"}
    if not truth or not required_truth <= set(truth[0]):
        raise SystemExit(
            "truth TSV must contain candidate_id, label and group_id"
        )
    by_id: dict[str, dict[str, str]] = {}
    model_ids: set[str] = set()
    calibration_states: set[str] = set()
    for row in features:
        candidate_id = row.get("candidate_id", "")
        if not candidate_id or candidate_id in by_id:
            raise SystemExit(f"missing or duplicate feature candidate_id {candidate_id!r}")
        by_id[candidate_id] = row
        model_ids.add(row.get("model_id", ""))
        calibration_states.add(row.get("calibration_status", ""))
    if calibration_states != {"UNCALIBRATED"}:
        raise SystemExit(
            "score audit expects explicitly UNCALIBRATED feature rows"
        )

    records: list[tuple[float, int, str]] = []
    seen_truth: set[str] = set()
    for row in truth:
        candidate_id = row["candidate_id"]
        if candidate_id in seen_truth:
            raise SystemExit(f"duplicate truth candidate_id {candidate_id!r}")
        seen_truth.add(candidate_id)
        if candidate_id not in by_id:
            raise SystemExit(f"truth candidate {candidate_id!r} has no feature row")
        group = row["group_id"].strip()
        if not group or group == ".":
            raise SystemExit(f"truth candidate {candidate_id!r} has no group_id")
        records.append(
            (
                finite_score(by_id[candidate_id]["membership_score"]),
                binary_label(row["label"]),
                group,
            )
        )
    groups = sorted({group for _, _, group in records})
    overall, overall_bins = summarize(records, args.bins)
    group_metrics: dict[str, dict[str, object]] = {}
    bin_scopes: list[tuple[str, list[dict[str, object]]]] = [
        ("overall", overall_bins)
    ]
    for group in groups:
        subset = [record for record in records if record[2] == group]
        metrics, group_bins = summarize(subset, args.bins)
        group_metrics[group] = metrics
        bin_scopes.append((f"group:{group}", group_bins))

    metrics_path = Path(f"{args.output}.metrics.json")
    bins_path = Path(f"{args.output}.calibration.tsv")
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "audit_version": "1.0.0",
        "model_ids": sorted(model_ids),
        "input_calibration_status": "UNCALIBRATED",
        "report_status": "EVALUATION_ONLY_NOT_A_CALIBRATED_MODEL",
        "truth_groups": groups,
        "overall": overall,
        "group_wise_evaluation": group_metrics,
        "features": str(args.features),
        "truth": str(args.truth),
    }
    with metrics_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_bins(bins_path, bin_scopes)
    print(f"Wrote {metrics_path} and {bins_path}")


if __name__ == "__main__":
    main()
