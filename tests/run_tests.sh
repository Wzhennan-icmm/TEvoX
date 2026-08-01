#!/usr/bin/env bash
set -euo pipefail

repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
bin="$repo/tevox"
data="$repo/tests/data"
work=$(mktemp -d "${TMPDIR:-/tmp}/tevox.XXXXXX")
trap 'rm -rf "$work"' EXIT

any_row() {
    python3 - "$@" <<'PY'
import csv
import sys

path, *pairs = sys.argv[1:]
if len(pairs) % 2:
    raise SystemExit("column/value arguments must be paired")
expected = dict(zip(pairs[::2], pairs[1::2]))
with open(path, newline="", encoding="utf-8") as handle:
    records = csv.DictReader(handle, delimiter="\t")
    if not any(all(row.get(key) == value for key, value in expected.items()) for row in records):
        raise SystemExit(f"no row in {path} matched {expected}")
PY
}

no_row() {
    if any_row "$@" 2>/dev/null; then
        echo "unexpected matching row: $*" >&2
        return 1
    fi
}

run_pair() {
    local prefix=$1 fasta_a=$2 te_a=$3 fasta_b=$4 te_b=$5 paf=$6
    "$bin" pair \
        --genome-a A --fasta-a "$fasta_a" --te-a "$te_a" \
        --genome-b B --fasta-b "$fasta_b" --te-b "$te_b" \
        --paf "$paf" --flank 20 --candidate-window 10 \
        --output "$prefix" >/dev/null
}

"$bin" --version | grep -F '0.4.0-alpha.1 (schema 1.1.0)' >/dev/null

# Non-finite numeric values must not bypass range validation.
for value in nan inf -inf; do
    if "$bin" pair \
        --genome-a A --fasta-a "$data/pair/A.fa" --te-a "$data/pair/A.gff3" \
        --genome-b B --fasta-b "$data/pair/B.fa" --te-b "$data/pair/B.gff3" \
        --paf "$data/pair/A_B.paf" --min-edge "$value" \
        --output "$work/nonfinite" >/dev/null 2>&1; then
        echo "accepted non-finite numeric option: $value" >&2
        exit 1
    fi
done

# Baseline pair: annotated match plus a cs-supported empty site.
run_pair "$work/pair" \
    "$data/pair/A.fa" "$data/pair/A.gff3" \
    "$data/pair/B.fa" "$data/pair/B.gff3" "$data/pair/A_B.paf"
any_row "$work/pair.states.tsv" genome_id B state EMPTY_SITE_CONFIRMED claimable true
any_row "$work/pair.edges.tsv" te_a A_shared te_b B_shared \
    independent_reciprocal false support_count 1 selected true
python3 "$repo/tests/validate_schema.py" "$work/pair"
diff -u "$repo/tests/golden/pair.states.tsv" "$work/pair.states.tsv"

# A single native PAF and its synthetic reverse share one dependency group and
# must not earn independent reciprocal support.
python3 - "$work/pair.evidence.tsv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
assert {row["origin"] for row in rows} == {"NATIVE", "DERIVED_REVERSE"}
groups = {row["evidence_group_id"] for row in rows}
assert len(groups) == 1
PY

# A separately generated native PAF in the opposite direction is independent
# reciprocal evidence and contributes a second deduplicated support group.
"$bin" graph --manifest "$data/v03/reciprocal/manifest.tsv" \
    --alignments "$data/v03/reciprocal/alignments.tsv" --flank 20 \
    --candidate-window 10 --output "$work/reciprocal" >/dev/null
any_row "$work/reciprocal.edges.tsv" te_a A_shared te_b B_shared \
    independent_reciprocal true support_count 2 selected true
python3 "$repo/tests/validate_schema.py" "$work/reciprocal"

# Annotation dropout, known family conflict, target gap and reverse strand.
run_pair "$work/unannotated" \
    "$data/pair/A.fa" "$data/pair/A.gff3" \
    "$data/pair/B.fa" "$data/unannotated.gff3" "$data/pair/A_B.paf"
any_row "$work/unannotated.states.tsv" genome_id B state PRESENT_UNANNOTATED claimable true

run_pair "$work/conflict" \
    "$data/conflict/A.fa" "$data/conflict/A.gff3" \
    "$data/conflict/B.fa" "$data/conflict/B.gff3" "$data/conflict/A_B.paf"
any_row "$work/conflict.states.tsv" genome_id B \
    state FAMILY_OR_BOUNDARY_DISCORDANCE annotation_state FAMILY_CONFLICT
no_row "$work/conflict.edges.tsv" selected true

run_pair "$work/gap" \
    "$data/conflict/A.fa" "$data/conflict/A.gff3" \
    "$data/gap/B.fa" "$data/conflict/B.gff3" "$data/conflict/A_B.paf"
any_row "$work/gap.states.tsv" genome_id B state ASSEMBLY_GAP technical_state GAP

run_pair "$work/reverse" \
    "$data/conflict/A.fa" "$data/conflict/A.gff3" \
    "$data/conflict/B.fa" "$data/reverse/B.gff3" "$data/reverse/A_B.paf"
any_row "$work/reverse.edges.tsv" te_a A_conflict te_b B_reverse \
    independent_reciprocal false selected true

# Local identity must come from cs/=X, never from whole-record PAF identity.
run_pair "$work/local_cs" \
    "$data/conflict/A.fa" "$data/conflict/A.gff3" \
    "$data/conflict/B.fa" "$data/conflict/B.gff3" "$data/v03/local_cs.paf"
any_row "$work/local_cs.evidence.tsv" source_te_id A_conflict \
    local_identity 0.983333 identity_method CS

run_pair "$work/missing_identity" \
    "$data/conflict/A.fa" "$data/conflict/A.gff3" \
    "$data/conflict/B.fa" "$data/conflict/B.gff3" \
    "$data/v03/missing_identity.paf"
any_row "$work/missing_identity.evidence.tsv" source_te_id A_conflict \
    local_identity . identity_method MISSING

# PAF MAPQ 255 is missing, and neither it nor a low native MAPQ can be rescued
# by the generated reverse view.
run_pair "$work/mapq255" \
    "$data/conflict/A.fa" "$data/conflict/A.gff3" \
    "$data/conflict/B.fa" "$data/conflict/B.gff3" "$data/v03/mapq255.paf"
any_row "$work/mapq255.evidence.tsv" origin NATIVE mapq . \
    mapq_status MISSING_255 technical_state UNCALLABLE
any_row "$work/mapq255.evidence.tsv" origin DERIVED_REVERSE mapq . \
    mapq_status MISSING_255 technical_state UNCALLABLE
no_row "$work/mapq255.edges.tsv" selected true

run_pair "$work/low_mapq" \
    "$data/conflict/A.fa" "$data/conflict/A.gff3" \
    "$data/conflict/B.fa" "$data/conflict/B.gff3" "$data/v03/low_mapq.paf"
any_row "$work/low_mapq.evidence.tsv" origin NATIVE \
    decision_code MAPQ_BELOW_THRESHOLD technical_state UNCALLABLE
any_row "$work/low_mapq.evidence.tsv" origin DERIVED_REVERSE \
    decision_code MAPQ_BELOW_THRESHOLD technical_state UNCALLABLE
no_row "$work/low_mapq.edges.tsv" selected true

# A missing flank at a contig edge is NA/uncallable, not perfect support.
run_pair "$work/edge" \
    "$data/v03/edge/A.fa" "$data/v03/edge/A.gff3" \
    "$data/v03/edge/B.fa" "$data/v03/edge/B.gff3" \
    "$data/v03/edge/A_B.paf"
any_row "$work/edge.evidence.tsv" left_flank . \
    left_flank_status CONTIG_EDGE decision_code SOURCE_CONTIG_EDGE \
    technical_state UNCALLABLE claimable false
no_row "$work/edge.edges.tsv" selected true

# Every near-best observation that triggers a source-level ambiguity must be
# linked by the decision external key, rather than retaining only the winner.
run_pair "$work/ambiguous" \
    "$data/v03/ambiguous/A.fa" "$data/v03/ambiguous/A.gff3" \
    "$data/v03/ambiguous/B.fa" "$data/v03/ambiguous/B.gff3" \
    "$data/v03/ambiguous/A_B.paf"
any_row "$work/ambiguous.decisions.tsv" source_te_id A_multi \
    technical_state AMBIGUOUS observation_count 2 near_best_count 2 \
    decision_code MULTIPLE_NEAR_BEST_PROJECTIONS claimable false
python3 - "$work/ambiguous.decisions.tsv" "$work/ambiguous.evidence.tsv" <<'PY'
import csv, sys
decisions = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
decision = next(row for row in decisions if row["source_te_id"] == "A_multi")
linked = set(decision["linked_evidence_ids"].split(","))
evidence = list(csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8"), delimiter="\t"))
observed = {row["evidence_id"] for row in evidence if row["decision_id"] == decision["decision_id"] and row["near_best"] == "true"}
assert linked == observed and len(linked) == 2
assert all(row["decision_ambiguous"] == "true" for row in evidence if row["evidence_id"] in linked)
PY
python3 "$repo/tests/validate_schema.py" "$work/ambiguous"

# A zero-overlap adjacent TE is not a co-ortholog merely because copy quota=2.
"$bin" graph --manifest "$data/multi/manifest.tsv" \
    --alignments "$data/multi/alignments.tsv" --flank 20 \
    --candidate-window 100 --output "$work/multi" >/dev/null
"$bin" graph --manifest "$data/multi/manifest1.tsv" \
    --alignments "$data/multi/alignments.tsv" --flank 20 \
    --candidate-window 100 --output "$work/multi_one" >/dev/null
no_row "$work/multi.instances.tsv" genome_id G2 copy_count 2
no_row "$work/multi_one.instances.tsv" genome_id G2 copy_count 2
any_row "$work/multi.candidates.tsv" target_te_id G2_A2 \
    reciprocal_overlap 0.000000 eligible false \
    decision_code RECIPROCAL_OVERLAP_LT_0.50
test "$(($(wc -l < "$work/multi.loci.tsv") - 1))" -eq 3
python3 "$repo/tests/validate_schema.py" "$work/multi"

# An unknown-family node cannot bridge two incompatible known families.
"$bin" graph --manifest "$data/v03/bridge/manifest.tsv" \
    --alignments "$data/v03/bridge/alignments.tsv" --flank 20 \
    --candidate-window 10 --output "$work/bridge" >/dev/null
any_row "$work/bridge.edges.tsv" selection_reason COMPONENT_FAMILY_CONFLICT \
    selected false
test "$(($(wc -l < "$work/bridge.loci.tsv") - 1))" -eq 2
if grep -E 'FamA.*FamC|FamC.*FamA' "$work/bridge.loci.tsv" >/dev/null; then
    echo "unknown family bridged incompatible known families" >&2
    exit 1
fi
any_row "$work/bridge.instances.tsv" decision_code TARGET_ASSIGNED_TO_DIFFERENT_LOCUS \
    claimable false
python3 "$repo/tests/validate_schema.py" "$work/bridge"
"$bin" graph --manifest "$data/v03/bridge/manifest.sorted.tsv" \
    --alignments "$data/v03/bridge/alignments.sorted.tsv" --flank 20 \
    --candidate-window 10 --output "$work/bridge_sorted" >/dev/null
for suffix in loci instances edges decisions candidates; do
    diff -u "$work/bridge.$suffix.tsv" "$work/bridge_sorted.$suffix.tsv"
done

# Annotation and PAF record order do not change inference or stable IDs.
run_pair "$work/order_a" \
    "$data/v03/order/A.fa" "$data/v03/order/A.gff3" \
    "$data/v03/order/B.fa" "$data/v03/order/B.gff3" \
    "$data/v03/order/A_B.paf"
run_pair "$work/order_b" \
    "$data/v03/order/A.fa" "$data/v03/order/A.reverse.gff3" \
    "$data/v03/order/B.fa" "$data/v03/order/B.reverse.gff3" \
    "$data/v03/order/A_B.reverse.paf"
for suffix in loci instances edges decisions candidates; do
    diff -u "$work/order_a.$suffix.tsv" "$work/order_b.$suffix.tsv"
done

# MUMmer4 NUCMER delta is normalized into the same evidence chain. Aggregate
# delta identity is not mislabelled as local identity or MAPQ, and a generated
# reverse view remains dependent on the same evidence group.
"$bin" graph --manifest "$data/v04/delta/manifest.tsv" \
    --alignments "$data/v04/delta/alignments.tsv" --flank 20 \
    --candidate-window 10 --output "$work/delta" >/dev/null
any_row "$work/delta.evidence.tsv" provider MUMMER_DELTA origin NATIVE \
    alignment_identity 1.000000 alignment_identity_method DELTA_ERROR_COUNT \
    local_identity . identity_method MISSING mapq . mapq_status NOT_PROVIDED \
    mapping_confidence UNIQUE_ALIGNMENT
any_row "$work/delta.evidence.tsv" provider MUMMER_DELTA \
    origin DERIVED_REVERSE dependency DERIVED_SAME_GROUP
no_row "$work/delta.states.tsv" state EMPTY_SITE_CONFIRMED
python3 "$repo/tests/validate_schema.py" "$work/delta"

"$bin" pair \
    --genome-a A --fasta-a "$data/conflict/A.fa" --te-a "$data/conflict/A.gff3" \
    --genome-b B --fasta-b "$data/conflict/B.fa" --te-b "$data/reverse/B.gff3" \
    --alignment "$data/v04/delta/full_reverse.delta" \
    --alignment-format mummer-delta --flank 20 --candidate-window 10 \
    --output "$work/delta_reverse" >/dev/null
any_row "$work/delta_reverse.edges.tsv" te_a A_conflict te_b B_reverse \
    independent_reciprocal false selected true
python3 "$repo/tests/validate_schema.py" "$work/delta_reverse"

# Canonical file-header paths take precedence over basename fallback. The
# swapped canonical direction must be rejected even when contig IDs and lengths
# are identical in the two genomes.
python3 - "$data/v04/delta/full.delta" "$data/conflict/B.fa" \
    "$data/conflict/A.fa" "$work/delta_canonical.delta" \
    "$work/delta_canonical_swapped.delta" <<'PY'
from pathlib import Path
import sys

template, target, query, correct, swapped = map(Path, sys.argv[1:])
lines = template.read_text(encoding="utf-8").splitlines()
lines[0] = f"{target.resolve()} {query.resolve()}"
correct.write_text("\n".join(lines) + "\n", encoding="utf-8")
lines[0] = f"{query.resolve()} {target.resolve()}"
swapped.write_text("\n".join(lines) + "\n", encoding="utf-8")
PY
"$bin" pair \
    --genome-a A --fasta-a "$data/conflict/A.fa" --te-a "$data/conflict/A.gff3" \
    --genome-b B --fasta-b "$data/conflict/B.fa" --te-b "$data/reverse/B.gff3" \
    --alignment "$work/delta_canonical.delta" --alignment-format delta \
    --output "$work/delta_canonical" >/dev/null
if "$bin" pair \
    --genome-a A --fasta-a "$data/conflict/A.fa" --te-a "$data/conflict/A.gff3" \
    --genome-b B --fasta-b "$data/conflict/B.fa" --te-b "$data/reverse/B.gff3" \
    --alignment "$work/delta_canonical_swapped.delta" --alignment-format delta \
    --output "$work/delta_canonical_swapped" >/dev/null 2>&1; then
    echo "accepted swapped canonical NUCMER reference/query paths" >&2
    exit 1
fi

# If neither header token resolves, basename fallback is legal only when the
# expected target/query FASTA basenames are distinct.
same_name="tevox_delta_same_${BASHPID}.fa"
mkdir -p "$work/same_query" "$work/same_target"
ln -s "$data/conflict/A.fa" "$work/same_query/$same_name"
ln -s "$data/conflict/B.fa" "$work/same_target/$same_name"
python3 - "$data/v04/delta/full.delta" "$same_name" \
    "$work/delta_ambiguous_header.delta" <<'PY'
from pathlib import Path
import sys

template, basename, output = sys.argv[1:]
lines = Path(template).read_text(encoding="utf-8").splitlines()
lines[0] = f"{basename} {basename}"
Path(output).write_text("\n".join(lines) + "\n", encoding="utf-8")
PY
if "$bin" pair \
    --genome-a A --fasta-a "$work/same_query/$same_name" \
    --te-a "$data/conflict/A.gff3" \
    --genome-b B --fasta-b "$work/same_target/$same_name" \
    --te-b "$data/reverse/B.gff3" \
    --alignment "$work/delta_ambiguous_header.delta" --alignment-format delta \
    --output "$work/delta_ambiguous_header" >/dev/null 2>&1; then
    echo "accepted ambiguous unresolved NUCMER header basenames" >&2
    exit 1
fi

for fixture in invalid_promer invalid_unterminated invalid_span invalid_errors \
    invalid_direction; do
    if "$bin" pair \
        --genome-a A --fasta-a "$data/conflict/A.fa" --te-a "$data/conflict/A.gff3" \
        --genome-b B --fasta-b "$data/conflict/B.fa" --te-b "$data/conflict/B.gff3" \
        --alignment "$data/v04/delta/$fixture.delta" --alignment-format delta \
        --output "$work/$fixture" >/dev/null 2>&1; then
        echo "accepted malformed MUMmer delta fixture: $fixture" >&2
        exit 1
    fi
done

# MCScanX is a synteny prior only. It builds stable WGD copy contexts but, in
# the absence of base alignment evidence, cannot produce an empty-site call.
"$bin" graph --manifest "$data/v04/synteny/manifest.tsv" \
    --synteny "$data/v04/synteny/sources.tsv" \
    --output "$work/synteny" >/dev/null
no_row "$work/synteny.states.tsv" state EMPTY_SITE_CONFIRMED
test "$(wc -l < "$work/synteny.evidence.tsv")" -eq 1
python3 "$repo/tests/validate_schema.py" "$work/synteny"
python3 - "$work/synteny.contexts.tsv" "$work/synteny.te_contexts.tsv" <<'PY'
import csv, sys
contexts = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
assignments = list(csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8"), delimiter="\t"))
assert len(contexts) == 3
assert len({row["homology_group_id"] for row in contexts}) == 1
assert sorted(row["syntenic_copy_id"] for row in contexts if row["genome_id"] == "B") == ["copy001", "copy002"]
assert len(assignments) == 4
assert {row["assignment"] for row in assignments} == {"BRACKETED"}
PY

# Gene, block, anchor, and side order cannot change biological stable IDs.
"$bin" graph --manifest "$data/v04/synteny/manifest.tsv" \
    --synteny "$data/v04/synteny/sources.shuffled.tsv" \
    --output "$work/synteny_shuffled" >/dev/null
"$bin" graph --manifest "$data/v04/synteny/manifest.tsv" \
    --synteny "$data/v04/synteny/sources.swapped.tsv" \
    --output "$work/synteny_swapped" >/dev/null
for suffix in contexts te_contexts synteny.anchors; do
    diff -u "$work/synteny.$suffix.tsv" "$work/synteny_shuffled.$suffix.tsv"
    diff -u "$work/synteny.$suffix.tsv" "$work/synteny_swapped.$suffix.tsv"
done

for sources in sources.malformed.tsv sources.bad_prefix.tsv; do
    if "$bin" graph --manifest "$data/v04/synteny/manifest.tsv" \
        --synteny "$data/v04/synteny/$sources" \
        --output "$work/synteny_bad" >/dev/null 2>&1; then
        echo "accepted malformed MCScanX input: $sources" >&2
        exit 1
    fi
done

# Two explicit WGD copy slots may join one ancestral locus even when the
# legacy per-genome fallback quota is one. Same-context copies remain blocked.
"$bin" graph --manifest "$data/v04/synteny/manifest1.tsv" \
    --alignments "$data/v04/synteny/alignments.tsv" \
    --synteny "$data/v04/synteny/sources.tsv" --flank 20 \
    --candidate-window 20 --output "$work/combined" >/dev/null
any_row "$work/combined.instances.tsv" genome_id B copy_count 2
any_row "$work/combined.candidates.tsv" context_relation SUPPORTED \
    shared_homology_group_id HMGfd5f79567dd2c40d context_compatible true
test "$(awk -F '\t' 'NR > 1 && $13 == "true" {n++} END {print n+0}' \
    "$work/combined.edges.tsv")" -eq 2
python3 "$repo/tests/validate_schema.py" "$work/combined"

# Two otherwise supported edges that would place two A annotations from the
# same strong copy context into one locus are rejected during graph merging.
"$bin" graph \
    --manifest "$data/v04/synteny/manifest.copy_context_conflict.tsv" \
    --alignments "$data/v04/synteny/alignments.copy_context_conflict.tsv" \
    --synteny "$data/v04/synteny/sources.tsv" --flank 20 \
    --candidate-window 20 --output "$work/copy_context_conflict" >/dev/null
any_row "$work/copy_context_conflict.candidate_contexts.tsv" \
    source_te_id A_TE target_te_id B1_TE context_relation SUPPORTED \
    context_compatible true
any_row "$work/copy_context_conflict.edges.tsv" te_a A_TE te_b B1_TE \
    selected false selection_reason COPY_CONTEXT_CONFLICT
any_row "$work/copy_context_conflict.edges.tsv" te_a A_TE_same_context \
    te_b B1_TE selected true selection_reason SELECTED
python3 "$repo/tests/validate_schema.py" "$work/copy_context_conflict"

# A context-free bridge may provide individually valid DNA edges, but it must
# not transitively merge TEs assigned to two disconnected HMGs.
"$bin" graph --manifest "$data/v04/synteny/manifest.hmg_conflict.tsv" \
    --alignments "$data/v04/synteny/alignments.hmg_conflict.tsv" \
    --synteny "$data/v04/synteny/sources.hmg_conflict.tsv" --flank 20 \
    --candidate-window 20 --output "$work/hmg_conflict" >/dev/null
any_row "$work/hmg_conflict.candidate_contexts.tsv" \
    source_te_id A_HMG1_TE target_te_id B_CONTEXT_FREE_TE \
    context_relation UNKNOWN context_compatible true
any_row "$work/hmg_conflict.candidate_contexts.tsv" \
    source_te_id A_HMG2_TE target_te_id B_CONTEXT_FREE_TE \
    context_relation UNKNOWN context_compatible true
any_row "$work/hmg_conflict.edges.tsv" te_a A_HMG1_TE \
    te_b B_CONTEXT_FREE_TE selected true selection_reason SELECTED
any_row "$work/hmg_conflict.edges.tsv" te_a A_HMG2_TE \
    te_b B_CONTEXT_FREE_TE selected false \
    selection_reason HOMOLOGY_GROUP_CONFLICT
python3 "$repo/tests/validate_schema.py" "$work/hmg_conflict"

# Conflicting context metadata is explicit ambiguity, never strong support.
"$bin" graph --manifest "$data/v04/synteny/manifest1.tsv" \
    --alignments "$data/v04/synteny/alignments.tsv" \
    --synteny "$data/v04/synteny/sources.metadata_conflict.tsv" --flank 20 \
    --candidate-window 20 --output "$work/context_ambiguous" >/dev/null
any_row "$work/context_ambiguous.contexts.tsv" genome_id A \
    status METADATA_AMBIGUOUS subgenome_id . haplotype_id hap1
any_row "$work/context_ambiguous.candidates.tsv" context_relation AMBIGUOUS \
    context_compatible false decision_code SYNTENY_CONTEXT_AMBIGUOUS
no_row "$work/context_ambiguous.candidates.tsv" context_relation SUPPORTED
python3 "$repo/tests/validate_schema.py" "$work/context_ambiguous"

# Overlapping evidence assigned to different WGD layers must remain separate.
"$bin" graph --manifest "$data/v04/synteny/manifest.tsv" \
    --synteny "$data/v04/synteny/sources.two_wgd.tsv" \
    --output "$work/two_wgd" >/dev/null
python3 - "$work/two_wgd.contexts.tsv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
assert len(rows) == 6
assert len({row["context_id"] for row in rows}) == 6
assert {row["wgd_node"] for row in rows} == {"WGD1", "WGD2"}
assert len({row["homology_group_id"] for row in rows}) == 2
PY
python3 "$repo/tests/validate_schema.py" "$work/two_wgd"

# Bounded top-K and unlimited candidates use the same complete winner order.
for limit in 1 0; do
    "$bin" pair \
        --genome-a A --fasta-a "$data/conflict/A.fa" --te-a "$data/conflict/A.gff3" \
        --genome-b B --fasta-b "$data/conflict/B.fa" --te-b "$data/v04/topk/B.gff3" \
        --paf "$data/conflict/A_B.paf" --flank 20 --candidate-window 10 \
        --max-candidates "$limit" --output "$work/topk_$limit" >/dev/null
    python3 "$repo/tests/validate_schema.py" "$work/topk_$limit"
done
for suffix in decisions edges loci instances states; do
    diff -u "$work/topk_1.$suffix.tsv" "$work/topk_0.$suffix.tsv"
done
python3 - "$work/topk_1.evidence.tsv" "$work/topk_1.candidates.tsv" \
    "$work/topk_0.evidence.tsv" "$work/topk_0.candidates.tsv" <<'PY'
import csv, sys
def read(path):
    return list(csv.DictReader(open(path, newline="", encoding="utf-8"), delimiter="\t"))
ev1, ca1, ev0, ca0 = map(read, sys.argv[1:])
a1 = next(row for row in ev1 if row["source_te_id"] == "A_conflict")
a0 = next(row for row in ev0 if row["source_te_id"] == "A_conflict")
assert a1["nearby_candidate_count"] == "2" and a1["retained_candidate_count"] == "1"
assert a0["nearby_candidate_count"] == "2" and a0["retained_candidate_count"] == "2"
winner1 = next(row for row in ca1 if row["evidence_id"] == a1["evidence_id"] and row["selected"] == "true")
winner0 = next(row for row in ca0 if row["evidence_id"] == a0["evidence_id"] and row["selected"] == "true")
assert winner1["target_te_id"] == winner0["target_te_id"] == a1["selected_target_te_id"]
PY

# Declared PAF query/target direction is validated against FASTA lengths.
if "$bin" pair \
    --genome-a B --fasta-a "$data/pair/B.fa" --te-a "$data/pair/B.gff3" \
    --genome-b A --fasta-b "$data/pair/A.fa" --te-b "$data/pair/A.gff3" \
    --paf "$data/pair/A_B.paf" --output "$work/wrong_direction" \
    >/dev/null 2>&1; then
    echo "accepted a PAF with reversed declared query/target genomes" >&2
    exit 1
fi

if "$bin" pair \
    --genome-a A --fasta-a "$data/conflict/A.fa" --te-a "$data/conflict/A.gff3" \
    --genome-b B --fasta-b "$data/conflict/B.fa" --te-b "$data/conflict/B.gff3" \
    --paf "$data/v03/incompatible_cs.paf" \
    --output "$work/incompatible_cs" >/dev/null 2>&1; then
    echo "accepted incompatible cg:Z and cs:Z paths" >&2
    exit 1
fi

if "$bin" pair \
    --genome-a A --fasta-a "$data/conflict/A.fa" --te-a "$data/conflict/A.gff3" \
    --genome-b B --fasta-b "$data/conflict/B.fa" --te-b "$data/conflict/B.gff3" \
    --paf "$data/invalid.paf" --output "$work/invalid" >/dev/null 2>&1; then
    echo "accepted invalid PAF" >&2
    exit 1
fi

python3 "$repo/scripts/tevox_phylo.py" \
    --states "$data/multi/phylo_states.tsv" --tree "$data/multi/tree.nwk" \
    --output "$work/phylo" >/dev/null
grep -Eq $'^TEL000001\t(gain|loss)\t' "$work/phylo.events.tsv"

python3 -m py_compile "$repo/tests/validate_schema.py" \
    "$repo/scripts/tevox_phylo.py"
echo 'All TEvoX tests passed.'
