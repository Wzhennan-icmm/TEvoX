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

"$bin" --version | grep -F '0.3.0-alpha.1 (schema 1.0.0)' >/dev/null

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
