#!/usr/bin/env python3
"""Validate TEvoX schema 1.1 provenance, keys, and missing-value semantics."""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


SCHEMA = "1.1.0"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def ids(value: str) -> list[str]:
    return [] if value in {"", "."} else value.split(",")


def require_schema(table: list[dict[str, str]], name: str) -> None:
    for row in table:
        assert row["schema_version"] == SCHEMA, (name, row)


def unique(table: list[dict[str, str]], field: str) -> set[str]:
    values = [row[field] for row in table]
    assert len(values) == len(set(values)), (field, values)
    return set(values)


def main(prefix_text: str) -> None:
    prefix = Path(prefix_text)
    evidence = rows(Path(f"{prefix}.evidence.tsv"))
    candidates = rows(Path(f"{prefix}.candidates.tsv"))
    candidate_contexts = rows(Path(f"{prefix}.candidate_contexts.tsv"))
    decisions = rows(Path(f"{prefix}.decisions.tsv"))
    edges = rows(Path(f"{prefix}.edges.tsv"))
    loci = rows(Path(f"{prefix}.loci.tsv"))
    instances = rows(Path(f"{prefix}.instances.tsv"))
    blocks = rows(Path(f"{prefix}.synteny.blocks.tsv"))
    anchors = rows(Path(f"{prefix}.synteny.anchors.tsv"))
    contexts = rows(Path(f"{prefix}.contexts.tsv"))
    te_contexts = rows(Path(f"{prefix}.te_contexts.tsv"))

    tables = {
        "evidence": evidence,
        "candidates": candidates,
        "candidate_contexts": candidate_contexts,
        "decisions": decisions,
        "edges": edges,
        "loci": loci,
        "instances": instances,
        "synteny.blocks": blocks,
        "synteny.anchors": anchors,
        "contexts": contexts,
        "te_contexts": te_contexts,
    }
    for name, table in tables.items():
        require_schema(table, name)

    evidence_ids = unique(evidence, "evidence_id")
    candidate_ids = unique(candidates, "candidate_id")
    decision_ids = unique(decisions, "decision_id")
    locus_ids = unique(loci, "locus_id")
    block_ids = unique(blocks, "block_id")
    anchor_ids = unique(anchors, "anchor_id")
    context_ids = unique(contexts, "context_id")
    te_context_ids = unique(te_contexts, "te_context_id")
    del anchor_ids, te_context_ids
    alignment_group_ids = {row["evidence_group_id"] for row in evidence}
    synteny_group_ids = {row["evidence_group_id"] for row in blocks}
    evidence_by_id = {row["evidence_id"]: row for row in evidence}
    context_by_id = {row["context_id"]: row for row in contexts}
    te_contexts_by_node: dict[tuple[str, str], list[dict[str, str]]] = (
        defaultdict(list)
    )
    for row in te_contexts:
        te_contexts_by_node[(row["genome_id"], row["te_id"])].append(row)
    homology_group_ids = {
        row["homology_group_id"] for row in contexts
        if row["homology_group_id"] != "."
    }

    for row in evidence:
        assert row["decision_id"] in decision_ids
        assert row["provider"] in {"PAF", "MUMMER_DELTA"}
        assert int(row["provider_record"]) >= 1
        assert int(row["provider_line"]) >= 1
        if row["identity_method"] == "MISSING":
            assert row["local_identity"] == "."
        else:
            assert row["local_identity"] != "."
        if row["alignment_identity_method"] == "MISSING":
            assert row["alignment_identity"] == "."
        else:
            assert row["alignment_identity"] != "."
        if row["mapq_status"] in {"MISSING_255", "NOT_PROVIDED"}:
            assert row["mapq"] == "."
        else:
            assert row["mapq_status"] == "OBSERVED"
            assert row["mapq"] != "."
        if row["provider"] == "MUMMER_DELTA":
            assert row["mapq_status"] == "NOT_PROVIDED"
            assert row["identity_method"] == "MISSING"
            assert row["alignment_identity_method"] == "DELTA_ERROR_COUNT"
            assert row["similarity_error_count"] != "."
            assert row["nonalpha_count"] != "."
        else:
            assert row["mapq_status"] in {"OBSERVED", "MISSING_255"}
            assert row["similarity_error_count"] == "."
            assert row["nonalpha_count"] == "."
        if row["origin"] == "DERIVED_REVERSE":
            assert row["dependency"] == "DERIVED_SAME_GROUP"
        else:
            assert row["origin"] == "NATIVE"
            assert row["dependency"] == "INDEPENDENT"
        if row["legacy_state"] == "EMPTY_SITE_CONFIRMED":
            assert row["provider"] == "PAF"
            assert row["local_identity"] != "."
            assert row["mapping_confidence"] != "NOT_ESTABLISHED"
        assert int(row["nearby_candidate_count"]) >= int(
            row["retained_candidate_count"]
        )

    assert unique(candidate_contexts, "candidate_id") == candidate_ids
    context_row_by_candidate = {
        row["candidate_id"]: row for row in candidate_contexts
    }
    ranks: dict[str, list[int]] = defaultdict(list)
    for row in candidates:
        assert row["evidence_id"] in evidence_ids
        assert row["decision_id"] in decision_ids
        assert row["candidate_rank"].isdigit()
        ranks[row["evidence_id"]].append(int(row["candidate_rank"]))
        context_row = context_row_by_candidate[row["candidate_id"]]
        assert context_row["evidence_id"] == row["evidence_id"]
        evidence_row = evidence_by_id[row["evidence_id"]]
        assert context_row["source_te_id"] == evidence_row["source_te_id"]
        assert context_row["target_te_id"] == row["target_te_id"]
        assert context_row["context_relation"] == row["context_relation"]
        assert context_row["context_compatible"] == row["context_compatible"]
        assert (
            context_row["shared_homology_group_id"]
            == row["shared_homology_group_id"]
        )
        shared = row["shared_homology_group_id"]
        if shared != ".":
            assert shared in homology_group_ids
        source_key = (
            evidence_row["query_genome_id"],
            context_row["source_te_id"],
        )
        target_key = (row["target_genome_id"], context_row["target_te_id"])
        source_assignments = te_contexts_by_node[source_key]
        target_assignments = te_contexts_by_node[target_key]
        source_context_ids = set(ids(context_row["source_context_ids"]))
        target_context_ids = set(ids(context_row["target_context_ids"]))
        assert source_context_ids == {
            item["context_id"] for item in source_assignments
        }
        assert target_context_ids == {
            item["context_id"] for item in target_assignments
        }
        assert source_context_ids <= context_ids
        assert target_context_ids <= context_ids

        source_hmgs = {
            context_by_id[value]["homology_group_id"]
            for value in source_context_ids
        }
        target_hmgs = {
            context_by_id[value]["homology_group_id"]
            for value in target_context_ids
        }
        source_strong_hmgs = {
            context_by_id[item["context_id"]]["homology_group_id"]
            for item in source_assignments
            if item["assignment"] == "BRACKETED"
            and context_by_id[item["context_id"]]["status"] == "PASS"
        }
        target_strong_hmgs = {
            context_by_id[item["context_id"]]["homology_group_id"]
            for item in target_assignments
            if item["assignment"] == "BRACKETED"
            and context_by_id[item["context_id"]]["status"] == "PASS"
        }
        source_ambiguous = any(
            item["assignment"] == "AMBIGUOUS"
            or (
                item["assignment"] == "BRACKETED"
                and context_by_id[item["context_id"]]["status"] != "PASS"
            )
            for item in source_assignments
        )
        target_ambiguous = any(
            item["assignment"] == "AMBIGUOUS"
            or (
                item["assignment"] == "BRACKETED"
                and context_by_id[item["context_id"]]["status"] != "PASS"
            )
            for item in target_assignments
        )
        relation = row["context_relation"]
        compatible = row["context_compatible"]
        assert relation in {"SUPPORTED", "CONFLICT", "AMBIGUOUS", "UNKNOWN"}
        assert compatible in {"true", "false"}
        if relation == "SUPPORTED":
            assert compatible == "true" and shared != "."
            assert shared in source_hmgs and shared in target_hmgs
            assert shared in source_strong_hmgs
            assert shared in target_strong_hmgs
        elif relation == "CONFLICT":
            assert compatible == "false" and shared == "."
            assert source_strong_hmgs and target_strong_hmgs
            assert source_strong_hmgs.isdisjoint(target_strong_hmgs)
        elif relation == "AMBIGUOUS":
            assert compatible == "false" and shared == "."
            assert source_ambiguous or target_ambiguous
        else:
            assert compatible == "true" and shared == "."
            assert not source_strong_hmgs or not target_strong_hmgs
    for values in ranks.values():
        assert sorted(values) == list(range(1, len(values) + 1))

    for row in decisions:
        winner = row["winner_evidence_id"]
        if winner != ".":
            assert winner in evidence_ids
        linked = ids(row["linked_evidence_ids"])
        assert all(value in evidence_ids for value in linked)
        assert len(linked) == int(row["near_best_count"])
        if row["technical_state"] == "AMBIGUOUS":
            assert linked or row["decision_code"] == "MULTIPLE_NEAR_BEST_CANDIDATES"

    for row in edges:
        assert all(value in evidence_ids for value in ids(row["evidence_ids"]))
        assert all(
            value in alignment_group_ids
            for value in ids(row["evidence_group_ids"])
        )
        assert int(row["support_count"]) == len(
            set(ids(row["evidence_group_ids"]))
        )

    for row in instances:
        assert row["locus_id"] in locus_ids
        assert all(
            value in decision_ids for value in ids(row["supporting_decision_ids"])
        )
        assert all(
            value in evidence_ids for value in ids(row["supporting_evidence_ids"])
        )
        if row["claimable"] == "false":
            assert row["claimability_reason"] not in {"", "."}

    anchor_counts = Counter(row["block_id"] for row in anchors)
    for row in blocks:
        assert row["provider"] == "MCScanX"
        assert row["context_a_id"] in context_ids
        assert row["context_b_id"] in context_ids
        context_a = context_by_id[row["context_a_id"]]
        context_b = context_by_id[row["context_b_id"]]
        assert row["homology_group_id"] == context_a["homology_group_id"]
        assert row["homology_group_id"] == context_b["homology_group_id"]
        assert (row["genome_a"], row["contig_a"]) == (
            context_a["genome_id"], context_a["contig"]
        )
        assert (row["genome_b"], row["contig_b"]) == (
            context_b["genome_id"], context_b["contig"]
        )
        assert row["wgd_node"] == context_a["wgd_node"]
        assert row["wgd_node"] == context_b["wgd_node"]
        assert int(row["anchor_count"]) == anchor_counts[row["block_id"]]
        assert row["homology_group_id"] in homology_group_ids
        block_anchors = [item for item in anchors if item["block_id"] == row["block_id"]]
        assert sorted(int(item["anchor_rank"]) for item in block_anchors) == list(
            range(1, len(block_anchors) + 1)
        )
        assert sorted(int(item["provider_anchor_rank"]) for item in block_anchors) == list(
            range(len(block_anchors))
        )
    assert all(row["block_id"] in block_ids for row in anchors)

    for row in contexts:
        assert row["homology_group_id"] != "."
        assert int(row["start"]) < int(row["end"])
        assert row["status"] in {"PASS", "METADATA_AMBIGUOUS"}
    for row in te_contexts:
        assert row["context_id"] in context_ids
        assert row["homology_group_id"] in homology_group_ids
        context = context_by_id[row["context_id"]]
        assert row["homology_group_id"] == context["homology_group_id"]
        assert row["genome_id"] == context["genome_id"]
        assert row["assignment"] in {"BRACKETED", "BLOCK_INTERIOR", "AMBIGUOUS"}
        assert 0.0 < float(row["overlap_fraction"]) <= 1.0

    with Path(f"{prefix}.run.json").open(encoding="utf-8") as handle:
        run = json.load(handle)
    assert run["schema_version"] == SCHEMA
    expected_counts = {
        "evidence_observations": len(evidence),
        "candidates": len(candidates),
        "decisions": len(decisions),
        "edges": len(edges),
        "loci": len(loci),
        "synteny_blocks": len(blocks),
        "synteny_anchors": len(anchors),
        "copy_contexts": len(contexts),
        "te_context_assignments": len(te_contexts),
    }
    for field, expected in expected_counts.items():
        assert run["counts"][field] == expected, (field, expected)
    assert run["counts"]["native_evidence_groups"] == len(
        run["alignment_evidence"]
    )
    assert run["counts"]["synteny_evidence_groups"] == len(
        run["synteny_evidence"]
    )
    assert alignment_group_ids <= {
        row["evidence_group_id"] for row in run["alignment_evidence"]
    }
    assert synteny_group_ids == {
        row["evidence_group_id"] for row in run["synteny_evidence"]
    }
    assert run["performance"]["fasta_reopens_after_index"] == 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: validate_schema.py OUTPUT_PREFIX")
    main(sys.argv[1])
