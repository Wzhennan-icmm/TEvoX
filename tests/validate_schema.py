#!/usr/bin/env python3
"""Validate TEvoX schema 1.2 provenance, inference, and missing semantics."""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path


SCHEMA = "1.2.0"


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


def normalized(row: dict[str, str], fields: list[str]) -> list[float]:
    values = [float(row[field]) for field in fields]
    assert all(0.0 <= value <= 1.0 for value in values), (row, fields)
    assert abs(sum(values) - 1.0) <= 2e-7, (row, fields, sum(values))
    return values


def prediction_members(value: str) -> set[str]:
    members = set(ids(value))
    assert members
    return members


def main(prefix_text: str) -> None:
    prefix = Path(prefix_text)
    evidence = rows(Path(f"{prefix}.evidence.tsv"))
    observation_scores = rows(Path(f"{prefix}.observation_scores.tsv"))
    candidates = rows(Path(f"{prefix}.candidates.tsv"))
    candidate_features = rows(Path(f"{prefix}.candidate_features.tsv"))
    candidate_contexts = rows(Path(f"{prefix}.candidate_contexts.tsv"))
    decisions = rows(Path(f"{prefix}.decisions.tsv"))
    edges = rows(Path(f"{prefix}.edges.tsv"))
    relations = rows(Path(f"{prefix}.relations.tsv"))
    solver = rows(Path(f"{prefix}.solver.tsv"))
    loci = rows(Path(f"{prefix}.loci.tsv"))
    instances = rows(Path(f"{prefix}.instances.tsv"))
    blocks = rows(Path(f"{prefix}.synteny.blocks.tsv"))
    anchors = rows(Path(f"{prefix}.synteny.anchors.tsv"))
    contexts = rows(Path(f"{prefix}.contexts.tsv"))
    te_contexts = rows(Path(f"{prefix}.te_contexts.tsv"))

    tables = {
        "evidence": evidence,
        "observation_scores": observation_scores,
        "candidates": candidates,
        "candidate_features": candidate_features,
        "candidate_contexts": candidate_contexts,
        "decisions": decisions,
        "edges": edges,
        "relations": relations,
        "solver": solver,
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
    assert unique(observation_scores, "evidence_id") == evidence_ids
    candidate_ids = unique(candidates, "candidate_id")
    assert unique(candidate_features, "candidate_id") == candidate_ids
    decision_ids = unique(decisions, "decision_id")
    locus_ids = unique(loci, "locus_id")
    edge_ids = unique(edges, "edge_id")
    relation_ids = unique(relations, "relation_id")
    solver_ids = unique(solver, "solver_component_id")
    block_ids = unique(blocks, "block_id")
    anchor_ids = unique(anchors, "anchor_id")
    context_ids = unique(contexts, "context_id")
    te_context_ids = unique(te_contexts, "te_context_id")
    del anchor_ids, te_context_ids, relation_ids
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
        assert int(row["nearby_candidate_count"]) >= int(
            row["graph_candidate_count"]
        )

    observation_by_id = {
        row["evidence_id"]: row for row in observation_scores
    }
    technical_fields = [
        "score_callable", "score_gap", "score_ambiguous", "score_uncallable"
    ]
    biological_fields = [
        "score_present", "score_empty", "score_structural_alternative",
        "score_biological_unknown",
    ]
    annotation_fields = [
        "score_annotation_matched", "score_annotation_missing",
        "score_family_conflict", "score_not_applicable",
        "score_annotation_unknown",
    ]
    for row in observation_scores:
        assert row["model_id"] == "BUILTIN_UNCALIBRATED_V1"
        assert row["calibration_status"] == "UNCALIBRATED"
        assert row["out_of_domain"] in {"true", "false"}
        normalized(row, technical_fields)
        normalized(row, biological_fields)
        normalized(row, annotation_fields)
        assert prediction_members(row["technical_prediction_set"]) <= {
            "CALLABLE", "GAP", "AMBIGUOUS", "UNCALLABLE"
        }
        assert prediction_members(row["biological_prediction_set"]) <= {
            "PRESENT", "EMPTY", "STRUCTURAL_ALTERNATIVE", "UNKNOWN"
        }
        assert prediction_members(row["annotation_prediction_set"]) <= {
            "MATCHED", "MISSING", "FAMILY_CONFLICT", "NOT_APPLICABLE",
            "UNKNOWN",
        }
        for field in (
            "technical_entropy", "biological_entropy", "annotation_entropy"
        ):
            assert 0.0 <= float(row[field]) <= 1.0

    feature_by_candidate = {
        row["candidate_id"]: row for row in candidate_features
    }
    for row in candidate_features:
        assert row["evidence_id"] in evidence_ids
        assert row["model_id"] == "BUILTIN_UNCALIBRATED_V1"
        assert row["calibration_status"] == "UNCALIBRATED"
        assert row["observed_feature_mask"].startswith("0x")
        assert row["out_of_domain"] in {"true", "false"}
        assert row["eligible"] in {"true", "false"}
        score = float(row["membership_score"])
        logit = float(row["membership_logit"])
        assert 0.0 < score < 1.0
        assert abs(score - 1.0 / (1.0 + math.exp(-logit))) <= 1e-8
        assert 0.0 <= float(row["membership_entropy"]) <= 1.0
        assert prediction_members(row["membership_prediction_set"]) <= {
            "SAME_LOCUS", "DIFFERENT_LOCUS"
        }
        missing = set(ids(row["missing_features"]))
        if "flank_min" in missing:
            assert row["flank_min"] == "."
        if "mapq" in missing:
            assert row["mapq_normalized"] == "."
        if "n_fraction" in missing:
            assert row["target_n_fraction"] == "."
        if "identity" in missing:
            assert row["local_identity"] == "."
            assert row["aggregate_identity"] == "."

    assert unique(candidate_contexts, "candidate_id") == candidate_ids
    context_row_by_candidate = {
        row["candidate_id"]: row for row in candidate_contexts
    }
    ranks: dict[str, list[int]] = defaultdict(list)
    for row in candidates:
        assert row["evidence_id"] in evidence_ids
        assert row["decision_id"] in decision_ids
        feature = feature_by_candidate[row["candidate_id"]]
        assert feature["evidence_id"] == row["evidence_id"]
        assert feature["eligible"] == row["eligible"]
        assert row["candidate_rank"].isdigit()
        assert row["graph_retained"] in {"true", "false"}
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
        assert row["model_id"] == "BUILTIN_UNCALIBRATED_V1"
        assert row["calibration_status"] == "UNCALIBRATED"
        assert row["out_of_domain"] in {"true", "false"}
        assert row["matching_selected"] in {"true", "false"}
        assert row["selected"] in {"true", "false"}
        assert row["context_relation"] in {
            "SUPPORTED", "CONFLICT", "AMBIGUOUS", "UNKNOWN"
        }
        if row["context_relation"] == "SUPPORTED":
            assert row["shared_homology_group_id"] in homology_group_ids
        else:
            assert row["shared_homology_group_id"] == "."
        score = float(row["membership_score"])
        logit = float(row["membership_logit"])
        assert 0.0 < score < 1.0
        assert abs(score - 1.0 / (1.0 + math.exp(-logit))) <= 1e-8
        assert 0.0 <= float(row["membership_entropy"]) <= 1.0
        method = row["matching_method"]
        assert method in {
            "NOT_APPLICABLE", "OPTIMAL_HUNGARIAN", "DETERMINISTIC_GREEDY"
        }
        if method == "NOT_APPLICABLE":
            assert row["matching_group_id"] == "."
        else:
            assert row["matching_group_id"].startswith("MAT")
        if row["selected"] == "true":
            assert row["matching_selected"] == "true"
            assert row["solver_component_id"] in solver_ids
            assert row["selection_reason"] in {
                "SELECTED_EXACT", "SELECTED_HEURISTIC"
            }
        elif row["solver_component_id"] != ".":
            assert row["solver_component_id"] in solver_ids

    edge_by_id = {row["edge_id"]: row for row in edges}
    relation_score_fields = [
        "score_ortholog", "score_wgd_homeolog", "score_allelic",
        "score_tandem_paralog", "score_segmental_paralog",
        "score_transposed_paralog", "score_unknown",
    ]
    relation_names = [
        "ORTHOLOG", "WGD_HOMEOLOG", "ALLELIC", "TANDEM_PARALOG",
        "SEGMENTAL_PARALOG", "TRANSPOSED_PARALOG", "UNKNOWN",
    ]
    for row in relations:
        assert row["model_id"] == "BUILTIN_UNCALIBRATED_V1"
        assert row["calibration_status"] == "UNCALIBRATED"
        assert row["out_of_domain"] in {"true", "false"}
        values = normalized(row, relation_score_fields)
        maximum = max(range(len(values)), key=values.__getitem__)
        assert row["predicted_relation"] == relation_names[maximum]
        assert prediction_members(row["prediction_set"]) <= set(relation_names)
        assert 0.0 <= float(row["entropy"]) <= 1.0
        if row["direct_edge"] == "true":
            assert row["edge_id"] in edge_ids
            edge = edge_by_id[row["edge_id"]]
            assert (row["te_a"], row["genome_a"], row["te_b"], row["genome_b"]) == (
                edge["te_a"], edge["genome_a"], edge["te_b"], edge["genome_b"]
            )
        else:
            assert row["edge_id"] == "."
        if row["locus_id"] != ".":
            assert row["locus_id"] in locus_ids

    for row in solver:
        assert int(row["node_count"]) >= 1
        assert int(row["edge_count"]) >= 0
        assert row["method"] in {
            "TRIVIAL", "EXACT_ENUMERATION", "DETERMINISTIC_GREEDY"
        }
        assert row["status"] in {"OPTIMAL", "HEURISTIC"}
        objective = float(row["objective"])
        bound = float(row["upper_bound"])
        gap = float(row["relative_gap"])
        assert bound + 1e-8 >= objective
        assert gap >= 0.0
        assert int(row["states_explored"]) >= 0
        if row["status"] == "OPTIMAL":
            assert row["method"] in {"TRIVIAL", "EXACT_ENUMERATION"}
            assert abs(bound - objective) <= 1e-8
            assert abs(gap) <= 1e-8
        else:
            assert row["method"] == "DETERMINISTIC_GREEDY"

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
    assert run["inference_model"] == {
        "model_id": "BUILTIN_UNCALIBRATED_V1",
        "calibration_status": "UNCALIBRATED",
        "score_semantics": "normalized_scores_not_calibrated_probabilities",
    }
    expected_counts = {
        "evidence_observations": len(evidence),
        "observation_score_rows": len(observation_scores),
        "candidates": len(candidates),
        "candidate_feature_rows": len(candidate_features),
        "decisions": len(decisions),
        "edges": len(edges),
        "relations": len(relations),
        "solver_components": len(solver),
        "loci": len(loci),
        "synteny_blocks": len(blocks),
        "synteny_anchors": len(anchors),
        "copy_contexts": len(contexts),
        "te_context_assignments": len(te_contexts),
    }
    for field, expected in expected_counts.items():
        assert run["counts"][field] == expected, (field, expected)
    assert run["counts"]["internal_candidates"] >= len(candidates)
    assert run["counts"]["internal_candidates"] >= run["counts"][
        "graph_candidates"
    ]
    assert run["counts"]["graph_candidates"] >= sum(
        row["graph_retained"] == "true" for row in candidates
    )
    assert run["config"]["max_candidates"] >= 0
    assert run["config"]["max_graph_candidates"] >= 0
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
    assert {
        "observation_scores.tsv", "candidate_features.tsv", "relations.tsv",
        "solver.tsv",
    } <= set(run["outputs"])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: validate_schema.py OUTPUT_PREFIX")
    main(sys.argv[1])
