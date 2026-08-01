#!/usr/bin/env python3
"""Validate TEvoX v0.3 evidence-schema foreign keys and missing semantics."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


SCHEMA = "1.0.0"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def ids(value: str) -> list[str]:
    return [] if value in {"", "."} else value.split(",")


def require_schema(table: list[dict[str, str]], name: str) -> None:
    for row in table:
        assert row["schema_version"] == SCHEMA, (name, row)


def main(prefix_text: str) -> None:
    prefix = Path(prefix_text)
    evidence = rows(Path(f"{prefix}.evidence.tsv"))
    candidates = rows(Path(f"{prefix}.candidates.tsv"))
    decisions = rows(Path(f"{prefix}.decisions.tsv"))
    edges = rows(Path(f"{prefix}.edges.tsv"))
    loci = rows(Path(f"{prefix}.loci.tsv"))
    instances = rows(Path(f"{prefix}.instances.tsv"))

    for name, table in {
        "evidence": evidence,
        "candidates": candidates,
        "decisions": decisions,
        "edges": edges,
        "loci": loci,
        "instances": instances,
    }.items():
        require_schema(table, name)

    evidence_ids = {row["evidence_id"] for row in evidence}
    group_ids = {row["evidence_group_id"] for row in evidence}
    decision_ids = {row["decision_id"] for row in decisions}
    locus_ids = {row["locus_id"] for row in loci}

    assert len(evidence_ids) == len(evidence)
    assert len(decision_ids) == len(decisions)

    for row in evidence:
        assert row["decision_id"] in decision_ids
        if row["identity_method"] == "MISSING":
            assert row["local_identity"] == "."
        if row["mapq_status"] == "MISSING_255":
            assert row["mapq"] == "."
        if row["origin"] == "DERIVED_REVERSE":
            assert row["dependency"] == "DERIVED_SAME_GROUP"

    for row in candidates:
        assert row["evidence_id"] in evidence_ids
        assert row["decision_id"] in decision_ids

    for row in decisions:
        winner = row["winner_evidence_id"]
        if winner != ".":
            assert winner in evidence_ids
        linked = ids(row["linked_evidence_ids"])
        assert all(value in evidence_ids for value in linked)
        assert len(linked) == int(row["near_best_count"])
        if row["technical_state"] == "AMBIGUOUS":
            assert len(linked) >= 2 or row["decision_code"] == "MULTIPLE_NEAR_BEST_CANDIDATES"

    for row in edges:
        assert all(value in evidence_ids for value in ids(row["evidence_ids"]))
        assert all(value in group_ids for value in ids(row["evidence_group_ids"]))
        assert int(row["support_count"]) == len(set(ids(row["evidence_group_ids"])))

    for row in instances:
        assert row["locus_id"] in locus_ids
        assert all(value in decision_ids for value in ids(row["supporting_decision_ids"]))
        assert all(value in evidence_ids for value in ids(row["supporting_evidence_ids"]))
        if row["claimable"] == "false":
            assert row["claimability_reason"] not in {"", "."}

    with Path(f"{prefix}.run.json").open(encoding="utf-8") as handle:
        run = json.load(handle)
    assert run["schema_version"] == SCHEMA
    assert run["counts"]["evidence_observations"] == len(evidence)
    assert run["counts"]["candidates"] == len(candidates)
    assert run["counts"]["decisions"] == len(decisions)
    assert run["counts"]["edges"] == len(edges)
    assert run["counts"]["loci"] == len(loci)
    assert run["counts"]["native_evidence_groups"] == len(run["alignment_evidence"])
    assert group_ids <= {
        row["evidence_group_id"] for row in run["alignment_evidence"]
    }


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: validate_schema.py OUTPUT_PREFIX")
    main(sys.argv[1])
