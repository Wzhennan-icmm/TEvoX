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

"$bin" --version | grep -F '0.5.0-alpha.2 (schema 1.2.0)' >/dev/null

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

# Output paths are preflighted against every frozen input and against one
# another. Neither an input PAF nor two hard-linked output names may be
# truncated before the run is rejected.
cp "$data/pair/A_B.paf" "$work/input_collision.evidence.tsv"
if "$bin" pair \
    --genome-a A --fasta-a "$data/pair/A.fa" --te-a "$data/pair/A.gff3" \
    --genome-b B --fasta-b "$data/pair/B.fa" --te-b "$data/pair/B.gff3" \
    --paf "$work/input_collision.evidence.tsv" \
    --output "$work/input_collision" >/dev/null 2>&1; then
    echo "allowed an output to overwrite its PAF input" >&2
    exit 1
fi
cmp "$work/input_collision.evidence.tsv" "$data/pair/A_B.paf"
test ! -e "$work/input_collision.run.json"

cp "$data/pair/A_B.paf" "$work/output_alias.evidence.tsv"
ln "$work/output_alias.evidence.tsv" "$work/output_alias.candidates.tsv"
if "$bin" pair \
    --genome-a A --fasta-a "$data/pair/A.fa" --te-a "$data/pair/A.gff3" \
    --genome-b B --fasta-b "$data/pair/B.fa" --te-b "$data/pair/B.gff3" \
    --paf "$data/pair/A_B.paf" --output "$work/output_alias" \
    >/dev/null 2>&1; then
    echo "allowed two output names to alias one inode" >&2
    exit 1
fi
cmp "$work/output_alias.evidence.tsv" "$data/pair/A_B.paf"
test ! -e "$work/output_alias.run.json"

ln -s "$work/dangling_output_target.tsv" \
    "$work/dangling_alias.evidence.tsv"
ln -s "$work/dangling_output_target.tsv" \
    "$work/dangling_alias.candidates.tsv"
if "$bin" pair \
    --genome-a A --fasta-a "$data/pair/A.fa" --te-a "$data/pair/A.gff3" \
    --genome-b B --fasta-b "$data/pair/B.fa" --te-b "$data/pair/B.gff3" \
    --paf "$data/pair/A_B.paf" --output "$work/dangling_alias" \
    >/dev/null 2>&1; then
    echo "allowed dangling output symlinks" >&2
    exit 1
fi
test ! -e "$work/dangling_output_target.tsv"
test ! -e "$work/dangling_alias.run.json"

# run.json is the atomic last-written completion marker. A mid-bundle write
# failure must invalidate the previous marker, and a later clean retry must
# restore it without leaving a temporary marker behind.
run_pair "$work/transactional" \
    "$data/pair/A.fa" "$data/pair/A.gff3" \
    "$data/pair/B.fa" "$data/pair/B.gff3" "$data/pair/A_B.paf"
test -s "$work/transactional.run.json"
if (
    ulimit -c 0
    trap '' XFSZ
    ulimit -f 1
    "$bin" pair \
        --genome-a A --fasta-a "$data/pair/A.fa" --te-a "$data/pair/A.gff3" \
        --genome-b B --fasta-b "$data/pair/B.fa" --te-b "$data/pair/B.gff3" \
        --paf "$data/pair/A_B.paf" --flank 10 \
        --output "$work/transactional"
) >/dev/null 2>&1; then
    echo "file-size-limited rerun unexpectedly completed" >&2
    exit 1
fi
test ! -e "$work/transactional.run.json"
run_pair "$work/transactional" \
    "$data/pair/A.fa" "$data/pair/A.gff3" \
    "$data/pair/B.fa" "$data/pair/B.gff3" "$data/pair/A_B.paf"
test -s "$work/transactional.run.json"
if compgen -G "$work/transactional.run.json.tmp.*" >/dev/null; then
    echo "atomic run marker left a temporary file" >&2
    exit 1
fi

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
any_row "$work/missing_identity.candidate_features.tsv" local_identity . \
    aggregate_identity 0.995000 model_id BUILTIN_UNCALIBRATED_V1 \
    calibration_status UNCALIBRATED
python3 "$repo/tests/validate_schema.py" "$work/missing_identity"

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

# Two adjacent target TEs can each cover exactly half of one projected source
# locus and therefore both pass reciprocal=0.5. A permissive copy quota must not
# let the shared source bridge those mutually non-overlapping annotations.
"$bin" graph --manifest "$data/v05/nonoverlap/manifest.tsv" \
    --alignments "$data/v05/nonoverlap/alignments.tsv" --flank 20 \
    --candidate-window 40 --output "$work/nonoverlap_bridge" >/dev/null
python3 - "$work/nonoverlap_bridge.evidence.tsv" \
    "$work/nonoverlap_bridge.candidates.tsv" \
    "$work/nonoverlap_bridge.edges.tsv" \
    "$work/nonoverlap_bridge.instances.tsv" <<'PY'
import csv, sys
read = lambda path: list(csv.DictReader(open(path, newline="", encoding="utf-8"), delimiter="\t"))
evidence, candidates, edges, instances = map(read, sys.argv[1:])
source_by_evidence = {row["evidence_id"]: row["source_te_id"] for row in evidence}
halves = [
    row for row in candidates
    if source_by_evidence[row["evidence_id"]] == "A_bridge"
    and row["target_te_id"] in {"B_left", "B_right"}
]
assert len(halves) == 2
assert {row["target_te_id"] for row in halves} == {"B_left", "B_right"}
assert all(row["reciprocal_overlap"] == "0.500000" for row in halves)
assert all(row["eligible"] == row["graph_retained"] == "true" for row in halves)
bridge_edges = [
    row for row in edges
    if "A_bridge" in {row["te_a"], row["te_b"]}
    and ({row["te_a"], row["te_b"]} & {"B_left", "B_right"})
]
assert len(bridge_edges) == 2
assert sum(row["selected"] == "true" for row in bridge_edges) == 1
rejected = next(row for row in bridge_edges if row["selected"] == "false")
assert rejected["selection_reason"] == "GLOBAL_CONSTRAINT_SEPARATED"
assert not any(
    row["genome_id"] == "B" and row["copy_count"] == "2"
    for row in instances
)
PY
test "$(($(wc -l < "$work/nonoverlap_bridge.loci.tsv") - 1))" -eq 2
python3 "$repo/tests/validate_schema.py" "$work/nonoverlap_bridge"

# An unknown-family node cannot bridge two incompatible known families.
"$bin" graph --manifest "$data/v03/bridge/manifest.tsv" \
    --alignments "$data/v03/bridge/alignments.tsv" --flank 20 \
    --candidate-window 10 --output "$work/bridge" >/dev/null
any_row "$work/bridge.edges.tsv" selection_reason GLOBAL_CONSTRAINT_SEPARATED \
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
for suffix in loci instances edges decisions candidates observation_scores \
    candidate_features relations solver; do
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
for suffix in loci instances edges decisions candidates observation_scores \
    candidate_features relations solver; do
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

# MCScanX control and provider files participate in the same collision guard.
# In particular, a gene table must never be replaced by contexts.tsv.
cp "$data/v04/synteny/genes.tsv" \
    "$work/synteny_gene_collision.contexts.tsv"
python3 - "$data/v04/synteny/blocks.collinearity" \
    "$work/synteny_gene_collision.contexts.tsv" \
    "$work/synteny_gene_collision.sources.tsv" <<'PY'
import sys
from pathlib import Path
blocks, genes, output = map(lambda value: Path(value).resolve(), sys.argv[1:])
output.write_text(
    "source_id\tformat\tcollinearity\tgenes\twgd_node\n"
    f"collision\tmcscanx\t{blocks}\t{genes}\tWGD1\n",
    encoding="utf-8",
)
PY
if "$bin" graph --manifest "$data/v04/synteny/manifest.tsv" \
    --synteny "$work/synteny_gene_collision.sources.tsv" \
    --output "$work/synteny_gene_collision" >/dev/null 2>&1; then
    echo "allowed contexts output to overwrite an MCScanX gene input" >&2
    exit 1
fi
cmp "$work/synteny_gene_collision.contexts.tsv" \
    "$data/v04/synteny/genes.tsv"
test ! -e "$work/synteny_gene_collision.run.json"

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
python3 - "$work/combined.edges.tsv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
assert sum(row["selected"] == "true" for row in rows) == 2
assert all(row["matching_method"] == "OPTIMAL_HUNGARIAN" for row in rows)
assert {row["matching_group_id"] for row in rows} == {
    "MAT3c5287a6ca49c840", "MAT6d922f3907b10ebb"
}
PY
any_row "$work/combined.relations.tsv" predicted_relation WGD_HOMEOLOG \
    calibration_status UNCALIBRATED out_of_domain false
python3 "$repo/tests/validate_schema.py" "$work/combined"
python3 - "$work/combined.run.json" <<'PY'
import hashlib, json, sys
run = json.load(open(sys.argv[1], encoding="utf-8"))
assert run["counts"]["synteny_provider_inputs"] == 1
assert len(run["synteny_inputs"]) == 1
item = run["synteny_inputs"][0]
for path_field, hash_field in (
    ("collinearity_path", "collinearity_sha256"),
    ("gene_table_path", "gene_table_sha256"),
):
    with open(item[path_field], "rb") as handle:
        observed = hashlib.sha256(handle.read()).hexdigest()
    assert item[hash_field] == observed
PY

# Coordinate separation is a hard negative only without explicit WGD copy
# evidence. Two non-overlapping TEs on one contig remain co-orthologous when
# unique PASS contexts place them in distinct copies of the same HMG/WGD node.
"$bin" graph \
    --manifest "$data/v04/synteny/manifest.same_contig_wgd.tsv" \
    --alignments "$data/v04/synteny/alignments.same_contig_wgd.tsv" \
    --synteny "$data/v04/synteny/sources.same_contig_wgd.tsv" --flank 20 \
    --candidate-window 20 --output "$work/same_contig_wgd" >/dev/null
any_row "$work/same_contig_wgd.instances.tsv" genome_id B copy_count 2
any_row "$work/same_contig_wgd.relations.tsv" te_a B1_TE te_b B2_TE \
    predicted_relation WGD_HOMEOLOG
python3 - "$work/same_contig_wgd.contexts.tsv" <<'PY'
import csv, sys
rows = [
    row for row in csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t")
    if row["genome_id"] == "B"
]
assert len(rows) == 2
assert {row["contig"] for row in rows} == {"chrB1"}
assert len({row["context_id"] for row in rows}) == 2
assert len({row["homology_group_id"] for row in rows}) == 1
assert {row["wgd_node"] for row in rows} == {"WGD1"}
assert all(row["status"] == "PASS" for row in rows)
PY
python3 "$repo/tests/validate_schema.py" "$work/same_contig_wgd"

# If a deliberately strict membership gate leaves two explicit WGD copies in
# separate loci, their physical proximity still cannot relabel them TANDEM.
"$bin" graph \
    --manifest "$data/v04/synteny/manifest.same_contig_wgd.tsv" \
    --alignments "$data/v04/synteny/alignments.same_contig_wgd.tsv" \
    --synteny "$data/v04/synteny/sources.same_contig_wgd.tsv" --flank 20 \
    --candidate-window 20 --min-membership 1 --tandem-distance 100000 \
    --output "$work/same_contig_wgd_gated" >/dev/null
any_row "$work/same_contig_wgd_gated.relations.tsv" \
    te_a B1_TE te_b B2_TE predicted_relation UNKNOWN
no_row "$work/same_contig_wgd_gated.relations.tsv" \
    te_a B1_TE te_b B2_TE predicted_relation TANDEM_PARALOG
python3 "$repo/tests/validate_schema.py" "$work/same_contig_wgd_gated"

# Different haplotypes of the same subgenome are allelic metadata, not WGD
# copy evidence. They cannot invoke the same-contig WGD merge exception.
"$bin" graph \
    --manifest "$data/v04/synteny/manifest.same_contig_wgd.tsv" \
    --alignments "$data/v04/synteny/alignments.same_contig_wgd.tsv" \
    --synteny "$data/v04/synteny/sources.same_contig_allelic.tsv" --flank 20 \
    --candidate-window 20 --output "$work/same_contig_allelic" >/dev/null
python3 - "$work/same_contig_allelic.edges.tsv" \
    "$work/same_contig_allelic.instances.tsv" \
    "$work/same_contig_allelic.contexts.tsv" <<'PY'
import csv, sys
edges = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
instances = list(csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8"), delimiter="\t"))
contexts = [
    row for row in csv.DictReader(open(sys.argv[3], newline="", encoding="utf-8"), delimiter="\t")
    if row["genome_id"] == "B"
]
assert sum(row["selected"] == "true" for row in edges) == 1
assert sum(row["selection_reason"] == "GLOBAL_CONSTRAINT_SEPARATED" for row in edges) == 1
assert max(int(row["copy_count"]) for row in instances if row["genome_id"] == "B") == 1
assert {row["subgenome_id"] for row in contexts} == {"B1"}
assert {row["haplotype_id"] for row in contexts} == {"hap1", "hap2"}
PY
python3 "$repo/tests/validate_schema.py" "$work/same_contig_allelic"

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
    selected false matching_selected false \
    selection_reason BLOCK_MATCHING_CONFLICT
any_row "$work/copy_context_conflict.edges.tsv" te_a A_TE_same_context \
    te_b B1_TE selected true matching_selected true \
    selection_reason SELECTED_EXACT
python3 "$repo/tests/validate_schema.py" "$work/copy_context_conflict"

# Both A annotations locally prefer B_MATCH_1, but the exact block objective is
# larger when A_MATCH_2 uses its second-ranked candidate B_MATCH_2. This proves
# matching operates over retained alternatives rather than selected winners.
"$bin" graph --manifest "$data/v05/matching/manifest.tsv" \
    --alignments "$data/v05/matching/alignments.tsv" \
    --synteny "$data/v05/matching/sources.tsv" --flank 20 \
    --candidate-window 30 --output "$work/matching_alternative" >/dev/null
python3 - "$work/matching_alternative.evidence.tsv" \
    "$work/matching_alternative.candidates.tsv" \
    "$work/matching_alternative.edges.tsv" <<'PY'
import csv, sys
evidence = {
    row["evidence_id"]: row
    for row in csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t")
}
candidates = list(csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8"), delimiter="\t"))
edges = list(csv.DictReader(open(sys.argv[3], newline="", encoding="utf-8"), delimiter="\t"))
for source in ("A_MATCH_1", "A_MATCH_2"):
    selected = [
        row for row in candidates
        if evidence[row["evidence_id"]]["source_te_id"] == source
        and evidence[row["evidence_id"]]["query_genome_id"] == "A"
        and row["selected"] == "true"
    ]
    assert len(selected) == 1 and selected[0]["target_te_id"] == "B_MATCH_1"
edge = next(row for row in edges if row["te_a"] == "A_MATCH_2" and row["te_b"] == "B_MATCH_2")
assert edge["matching_method"] == "OPTIMAL_HUNGARIAN"
assert edge["matching_selected"] == edge["selected"] == "true"
conflict = next(
    row for row in edges
    if row["te_a"] == "A_MATCH_1" and row["te_b"] == "B_MATCH_CONFLICT"
)
assert conflict["family_compatible"] == "false"
assert conflict["matching_method"] == "NOT_APPLICABLE"
assert conflict["matching_group_id"] == "."
assert conflict["selected"] == "false"
assert conflict["selection_reason"] == "DIRECT_FAMILY_CONFLICT"
assert next(
    row for row in candidates
    if evidence[row["evidence_id"]]["source_te_id"] == "A_MATCH_2"
    and row["target_te_id"] == "B_MATCH_2"
)["candidate_rank"] == "2"
PY
python3 "$repo/tests/validate_schema.py" "$work/matching_alternative"

# The block matcher reports its deterministic fallback instead of claiming an
# exact result when a context side exceeds the configured Hungarian limit.
"$bin" graph \
    --manifest "$data/v04/synteny/manifest.copy_context_conflict.tsv" \
    --alignments "$data/v04/synteny/alignments.copy_context_conflict.tsv" \
    --synteny "$data/v04/synteny/sources.tsv" --flank 20 \
    --candidate-window 20 --exact-match-nodes 1 \
    --output "$work/copy_context_heuristic" >/dev/null
any_row "$work/copy_context_heuristic.edges.tsv" \
    matching_method DETERMINISTIC_GREEDY matching_selected true \
    selected true selection_reason SELECTED_EXACT
python3 "$repo/tests/validate_schema.py" "$work/copy_context_heuristic"

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
    te_b B_CONTEXT_FREE_TE selected false \
    selection_reason GLOBAL_CONSTRAINT_SEPARATED
any_row "$work/hmg_conflict.edges.tsv" te_a A_HMG2_TE \
    te_b B_CONTEXT_FREE_TE selected true selection_reason SELECTED_EXACT
python3 "$repo/tests/validate_schema.py" "$work/hmg_conflict"

# A truth-controlled four-node component makes the greedy failure explicit:
# its strongest first edge blocks two compatible medium edges. Exact search
# recovers the higher global objective and reports a zero optimality gap.
"$bin" graph --manifest "$data/v05/global/manifest.tsv" \
    --alignments "$data/v05/global/alignments.tsv" --flank 20 \
    --candidate-window 10 --max-candidates 0 --max-graph-candidates 0 \
    --output "$work/global_exact" >/dev/null
"$bin" graph --manifest "$data/v05/global/manifest.tsv" \
    --alignments "$data/v05/global/alignments.tsv" --flank 20 \
    --candidate-window 10 --max-candidates 0 --max-graph-candidates 0 \
    --exact-max-edges 0 \
    --output "$work/global_heuristic" >/dev/null
any_row "$work/global_exact.edges.tsv" te_a A_global te_b B_bridge \
    selected true selection_reason SELECTED_EXACT
any_row "$work/global_exact.edges.tsv" te_a B_bridge te_b C_global \
    selected true selection_reason SELECTED_EXACT
any_row "$work/global_exact.solver.tsv" method EXACT_ENUMERATION \
    status OPTIMAL relative_gap 0.00000000
any_row "$work/global_heuristic.edges.tsv" te_a A_greedy te_b B_bridge \
    selected true selection_reason SELECTED_HEURISTIC
any_row "$work/global_heuristic.solver.tsv" \
    method DETERMINISTIC_GREEDY status HEURISTIC
python3 - "$work/global_exact.solver.tsv" \
    "$work/global_heuristic.solver.tsv" <<'PY'
import csv, sys
def objective(path):
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8"), delimiter="\t"))
    return sum(float(row["objective"]) for row in rows)
assert objective(sys.argv[1]) > objective(sys.argv[2])
PY
python3 "$repo/tests/validate_schema.py" "$work/global_exact"
python3 "$repo/tests/validate_schema.py" "$work/global_heuristic"

# A direct edge can fail its own gate while both endpoints are assigned to the
# same final locus through other selected edges. Relation locus membership must
# follow the final partition, not the direct edge's selected flag.
"$bin" graph --manifest "$data/v05/global/manifest.tsv" \
    --alignments "$data/v05/global/alignments.weak_triangle.tsv" --flank 20 \
    --candidate-window 10 --min-edge 80 \
    --output "$work/relation_weak_triangle" >/dev/null
any_row "$work/relation_weak_triangle.edges.tsv" te_a A_global te_b C_global \
    selected false selection_reason BELOW_EDGE_THRESHOLD
python3 - "$work/relation_weak_triangle.relations.tsv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
row = next(
    item for item in rows
    if {item["te_a"], item["te_b"]} == {"A_global", "C_global"}
)
assert row["direct_edge"] == "true"
assert row["locus_id"] != "."
assert row["predicted_relation"] == "ORTHOLOG"
PY
python3 "$repo/tests/validate_schema.py" "$work/relation_weak_triangle"

# Overlapping or nested annotations are not tandem copies. A nearby,
# non-overlapping pair remains eligible for the tandem relation label.
"$bin" graph --manifest "$data/v04/synteny/manifest.nested.tsv" \
    --alignments "$data/v04/synteny/alignments.tsv" --flank 20 \
    --candidate-window 20 --tandem-distance 100 \
    --output "$work/relation_nested" >/dev/null
python3 - "$work/relation_nested.relations.tsv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
by_pair = {frozenset((row["te_a"], row["te_b"])): row for row in rows}
assert frozenset(("A_TE", "A_TE_nested")) not in by_pair
nearby = by_pair[frozenset(("A_TE", "A_TE_same_context"))]
assert nearby["predicted_relation"] == "TANDEM_PARALOG"
PY
python3 "$repo/tests/validate_schema.py" "$work/relation_nested"

# The nested hard-negative applies to locus reconstruction, not only relation
# naming. A permissive legacy quota and a shared B bridge cannot merge two
# distinct overlapping A annotations without explicit WGD copy contexts.
"$bin" pair \
    --genome-a A --fasta-a "$data/v04/synteny/A.fa" \
    --te-a "$data/v04/synteny/A.nested.gff3" --max-copies-a 2 \
    --genome-b B --fasta-b "$data/v04/synteny/B.fa" \
    --te-b "$data/v04/synteny/B.gff3" --max-copies-b 1 \
    --paf "$data/v04/synteny/A_B.paf" --flank 20 --candidate-window 20 \
    --output "$work/nested_bridge" >/dev/null
python3 - "$work/nested_bridge.loci.tsv" \
    "$work/nested_bridge.instances.tsv" \
    "$work/nested_bridge.edges.tsv" <<'PY'
import csv, sys
read = lambda path: list(csv.DictReader(open(path, newline="", encoding="utf-8"), delimiter="\t"))
loci, instances, edges = map(read, sys.argv[1:])
for row in loci:
    members = set(row["members"].split(","))
    assert not {"A:A_TE", "A:A_TE_nested"} <= members
assert not any(
    row["genome_id"] == "A" and int(row["copy_count"]) > 1
    for row in instances
)
assert any(
    row["selection_reason"] == "GLOBAL_CONSTRAINT_SEPARATED"
    and ({row["te_a"], row["te_b"]} & {"A_TE", "A_TE_nested"})
    for row in edges
)
PY
python3 "$repo/tests/validate_schema.py" "$work/nested_bridge"

# The legacy run-local score audit remains an explicitly uncalibrated
# diagnostic. It now reports coverage, AUPRC, equal-mass ECE and deterministic
# group bootstrap intervals without fitting or blessing a model.
python3 - "$work/global_exact.evidence.tsv" \
    "$work/global_exact.candidates.tsv" "$work/global_truth.tsv" <<'PY'
import csv, sys
evidence = {
    row["evidence_id"]: row
    for row in csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t")
}
truth_pairs = {frozenset(("A_global", "B_bridge")), frozenset(("B_bridge", "C_global"))}
with open(sys.argv[3], "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["candidate_id", "label", "group_id"], delimiter="\t")
    writer.writeheader()
    for row in csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8"), delimiter="\t"):
        observation = evidence[row["evidence_id"]]
        pair = frozenset((observation["source_te_id"], row["target_te_id"]))
        writer.writerow({
            "candidate_id": row["candidate_id"],
            "label": int(pair in truth_pairs),
            "group_id": f'{observation["query_genome_id"]}-{observation["target_genome_id"]}',
        })
PY
python3 "$repo/scripts/tevox_score_audit.py" \
    --features "$work/global_exact.candidate_features.tsv" \
    --truth "$work/global_truth.tsv" --bins 4 --bootstrap-replicates 20 \
    --seed 17 \
    --output "$work/global_audit" >/dev/null
python3 - "$work/global_audit.metrics.json" \
    "$work/global_audit.calibration.tsv" <<'PY'
import csv, json, sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
assert report["input_calibration_status"] == "UNCALIBRATED"
assert report["report_status"] == "EVALUATION_ONLY_NOT_A_CALIBRATED_MODEL"
assert len(report["truth_groups"]) >= 2
assert report["overall"]["n"] > 0
assert report["overall"]["auprc"] is not None
assert report["overall"]["ece_equal_mass"] is not None
assert report["coverage"]["truth_feature_coverage"] == 1.0
assert report["group_bootstrap_95_ci"]["replicates_requested"] == 20
rows = list(csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8"), delimiter="\t"))
assert {row["scope"] for row in rows} >= {"overall"}
PY
python3 - "$work/global_exact.candidate_features.tsv" \
    "$work/global_mixed_model.candidate_features.tsv" <<'PY'
import csv, sys
with open(sys.argv[1], newline="", encoding="utf-8") as source:
    reader = csv.DictReader(source, delimiter="\t")
    fields, rows = reader.fieldnames, list(reader)
rows[0]["model_id"] = "OTHER_UNCALIBRATED_MODEL"
with open(sys.argv[2], "w", newline="", encoding="utf-8") as destination:
    writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
PY
if python3 "$repo/scripts/tevox_score_audit.py" \
    --features "$work/global_mixed_model.candidate_features.tsv" \
    --truth "$work/global_truth.tsv" --output "$work/mixed_model_audit" \
    >/dev/null 2>&1; then
    echo "score audit accepted multiple model_id values" >&2
    exit 1
fi

# Publication metrics must be invariant to truth-row order inside tied scores,
# and a truth-supported class that is completely missed contributes F1=0.
python3 - "$repo/scripts" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
import tevox_benchmark as benchmark

tied = [(0.5, index % 2) for index in range(10)]
ece, bins = benchmark.equal_mass_ece(tied, requested_bins=10)
reversed_ece, reversed_bins = benchmark.equal_mass_ece(
    list(reversed(tied)), requested_bins=10
)
assert ece == 0.0 and reversed_ece == 0.0
assert len(bins) == len(reversed_bins) == 1

partly_wrong = benchmark.confusion_metrics([("A", "B"), ("B", "B")])
assert partly_wrong["per_class"]["A"]["f1"] == 0.0
assert abs(partly_wrong["macro_f1"] - (1.0 / 3.0)) < 1e-12
all_wrong = benchmark.confusion_metrics([("A", "B"), ("B", "A")])
assert all_wrong["macro_f1"] == 0.0

truth_singletons = {("G", "a"): "T1", ("G", "b"): "T2"}
zero = benchmark.clustering_metrics(truth_singletons, {}, exhaustive=True)
assert zero["truth_member_prediction_coverage"] == 0.0
assert zero["b_cubed_f1"] == 0.0
assert zero["intersection_only_b_cubed"]["f1"] is None
try:
    benchmark.clustering_metrics(
        truth_singletons, {("G", "a"): "P1", ("G", "x"): "P2"},
        exhaustive=True,
    )
except benchmark.ContractError:
    pass
else:
    raise AssertionError("exhaustive locus truth accepted a prediction-only node")

tie_truth = {
    ("G", "a"): "T1", ("G", "b"): "T1",
    ("G", "c"): "T2", ("G", "d"): "T2",
}
tie_prediction = {
    ("G", "a"): "P1", ("G", "c"): "P1",
    ("G", "b"): "P2", ("G", "d"): "P2",
}
renamed_prediction = {
    key: {"P1": "Z", "P2": "A"}[value]
    for key, value in tie_prediction.items()
}
mapping, summary = benchmark.maximum_overlap_locus_mapping(
    tie_truth, tie_prediction
)
renamed, renamed_summary = benchmark.maximum_overlap_locus_mapping(
    tie_truth, renamed_prediction
)
assert mapping == renamed == {}
assert {row["status"] for row in summary["mapping_rows"]} == {
    "AMBIGUOUS_OR_WEAK"
}
assert {row["status"] for row in renamed_summary["mapping_rows"]} == {
    "AMBIGUOUS_OR_WEAK"
}

valid_coordinates = {"contig": "chr1", "start": "10", "end": "20"}
assert benchmark.validate_predicted_coordinates(valid_coordinates, 1, 2) == (10, 20)
for malformed, copies in (
    ({"contig": "chr1", "start": "10", "end": "."}, 1),
    ({"contig": "chr1", "start": "20", "end": "10"}, 1),
    ({"contig": "chr1", "start": "10,30", "end": "20,40"}, 2),
):
    try:
        benchmark.validate_predicted_coordinates(malformed, copies, 2)
    except benchmark.ContractError:
        pass
    else:
        raise AssertionError("accepted malformed or multi-copy state coordinates")
PY

# The publication-facing benchmark uses dataset semantic keys, never CAN/EVD/
# TEL identifiers. Exact and heuristic runs share one immutable truth bundle;
# the evaluator must expose the known greedy clustering failure.
python3 - "$work/benchmark_runs.tsv" "$work/global_exact" \
    "$work/global_heuristic" <<'PY'
import csv, sys
fields = [
    "benchmark_schema_version", "run_id", "dataset_id", "method_id",
    "condition_id", "replicate_id", "group_id", "prediction_format", "prefix",
]
with open(sys.argv[1], "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for run_id, condition, prefix in (
        ("exact", "exact", sys.argv[2]),
        ("heuristic", "heuristic", sys.argv[3]),
    ):
        writer.writerow({
            "benchmark_schema_version": "1.0.0", "run_id": run_id,
            "dataset_id": "D_GLOBAL", "method_id": "TEvoX",
            "condition_id": condition, "replicate_id": "1",
            "group_id": "controlled_fixture", "prediction_format": "tevox-1.2",
            "prefix": prefix,
        })
PY
python3 "$repo/scripts/tevox_benchmark.py" \
    --datasets "$data/v05/benchmark/datasets.tsv" \
    --runs "$work/benchmark_runs.tsv" --output "$work/benchmark" >/dev/null
python3 - "$work/benchmark.metrics.json" \
    "$work/benchmark.callability_accuracy.tsv" <<'PY'
import csv, json, sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
assert report["report_status"] == "EVALUATION_ONLY_NOT_A_CALIBRATED_MODEL"
assert report["truth_identity"] == "DATASET_SEMANTIC_KEYS_INDEPENDENT_OF_RUN_ID"
assert report["run_count"] == 2
runs = {row["run_id"]: row for row in report["runs"]}
exact = runs["exact"]
heuristic = runs["heuristic"]
assert exact["locus_clustering"]["b_cubed_f1"] == 1.0
assert exact["locus_clustering"]["adjusted_rand_index"] == 1.0
assert heuristic["locus_clustering"]["b_cubed_f1"] < 1.0
assert heuristic["locus_clustering"]["adjusted_rand_index"] < 1.0
assert exact["state_and_breakpoint"]["axes"]["biological_state"]["accuracy"] == 1.0
assert exact["candidate_membership"]["positive_candidate_generation_recall"] == 2 / 3
curves = list(csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8"), delimiter="\t"))
for run_id in {row["run_id"] for row in curves}:
    coverage = [float(row["coverage"]) for row in curves if row["run_id"] == run_id]
    assert coverage == sorted(coverage, reverse=True)
PY

# A header-only prediction is a valid catastrophic method, not a malformed
# benchmark input. Missing claimable=false must remain NO_PREDICTION, and an
# EMPTY call against biological UNKNOWN is unassessed rather than a false hit.
python3 - "$work/global_exact" "$work/benchmark_zero" \
    "$work/benchmark_unknown_empty" "$work/benchmark_instance_unknown.tsv" \
    "$work/benchmark_strict.datasets.tsv" "$work/benchmark_strict.runs.tsv" \
    "$data/v05/benchmark/truth.candidates.tsv" \
    "$data/v05/benchmark/truth.members.tsv" \
    "$data/v05/benchmark/truth.genomes.tsv" <<'PY'
import csv, shutil, sys
from pathlib import Path

source, zero, empty = map(Path, sys.argv[1:4])
instance_truth, datasets, runs = map(Path, sys.argv[4:7])
candidate_truth, locus_truth, genome_truth = map(
    lambda value: Path(value).resolve(), sys.argv[7:10]
)
for destination in (zero, empty):
    for suffix in ("run.json", "evidence.tsv", "candidates.tsv", "candidate_features.tsv"):
        shutil.copyfile(f"{source}.{suffix}", f"{destination}.{suffix}")

with open(f"{source}.instances.tsv", newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    instance_fields = reader.fieldnames
    instance_rows = list(reader)
with open(f"{zero}.instances.tsv", "w", newline="", encoding="utf-8") as handle:
    csv.DictWriter(handle, fieldnames=instance_fields, delimiter="\t").writeheader()

changed = 0
for row in instance_rows:
    if row["genome_id"] == "A" and row["member_ids"] == "A_global":
        row.update({
            "technical_state": "CALLABLE", "biological_state": "EMPTY",
            "annotation_state": "NOT_APPLICABLE",
            "legacy_state": "EMPTY_SITE_CONFIRMED", "claim_type": "EMPTY_SITE",
            "claimable": "true", "copy_count": "0", "member_ids": ".",
            "decision_code": "CONTROLLED_UNKNOWN_EMPTY_FIXTURE",
        })
        changed += 1
assert changed == 1
with open(f"{empty}.instances.tsv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=instance_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(instance_rows)

truth_fields = [
    "benchmark_schema_version", "dataset_id", "truth_record_id",
    "truth_locus_id", "genome_id", "technical_state", "biological_state",
    "annotation_state", "legacy_state", "claimable", "contig", "start", "end",
    "truth_confidence", "validation_method", "validation_source",
    "validation_batch_id",
]
with instance_truth.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=truth_fields, delimiter="\t")
    writer.writeheader()
    writer.writerow({
        "benchmark_schema_version": "1.0.0", "dataset_id": "D_GLOBAL",
        "truth_record_id": "INSTRUTH_UNKNOWN", "truth_locus_id": "TL_GLOBAL",
        "genome_id": "A", "technical_state": ".", "biological_state": "UNKNOWN",
        "annotation_state": ".", "legacy_state": ".", "claimable": "false",
        "contig": ".", "start": ".", "end": ".", "truth_confidence": "HIGH",
        "validation_method": "CONTROLLED_FIXTURE",
        "validation_source": "tests/data/v05/global",
        "validation_batch_id": "SIM_GLOBAL",
    })

dataset_fields = [
    "benchmark_schema_version", "dataset_id", "genome_truth",
    "candidate_truth", "locus_truth", "instance_truth",
    "candidate_truth_scope", "candidate_sampling_design",
    "candidate_inclusion_probability", "locus_truth_scope",
]
with datasets.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=dataset_fields, delimiter="\t")
    writer.writeheader()
    writer.writerow({
        "benchmark_schema_version": "1.0.0", "dataset_id": "D_GLOBAL",
        "genome_truth": genome_truth,
        "candidate_truth": candidate_truth, "locus_truth": locus_truth,
        "instance_truth": instance_truth.resolve(),
        "candidate_truth_scope": "EXHAUSTIVE",
        "candidate_sampling_design": "CENSUS",
        "candidate_inclusion_probability": "1", "locus_truth_scope": "EXHAUSTIVE",
    })
run_fields = [
    "benchmark_schema_version", "run_id", "dataset_id", "method_id",
    "condition_id", "replicate_id", "group_id", "prediction_format", "prefix",
]
with runs.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=run_fields, delimiter="\t")
    writer.writeheader()
    for run_id, prefix in (("zero", zero), ("unknown_empty", empty)):
        writer.writerow({
            "benchmark_schema_version": "1.0.0", "run_id": run_id,
            "dataset_id": "D_GLOBAL", "method_id": "TEvoX",
            "condition_id": run_id, "replicate_id": "1", "group_id": "strict",
            "prediction_format": "tevox-1.2", "prefix": prefix.resolve(),
        })
PY
python3 "$repo/scripts/tevox_benchmark.py" \
    --datasets "$work/benchmark_strict.datasets.tsv" \
    --runs "$work/benchmark_strict.runs.tsv" \
    --output "$work/benchmark_strict" >/dev/null
python3 - "$work/benchmark_strict.metrics.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
runs = {row["run_id"]: row for row in report["runs"]}
zero = runs["zero"]
assert zero["locus_clustering"]["truth_member_prediction_coverage"] == 0.0
assert zero["locus_clustering"]["b_cubed_f1"] == 0.0
claimable = zero["state_and_breakpoint"]["axes"]["claimable"]
assert claimable["accuracy"] == 0.0
assert claimable["confusion"]["false"]["NO_PREDICTION"] == 1
unknown = runs["unknown_empty"]["state_and_breakpoint"]
assert unknown["empty_site_claims_evaluated"] == 0
assert unknown["unassessed_empty_site_claims"] >= 1
assert unknown["empty_claim_assessment_coverage"] == 0.0
assert unknown["assessed_empty_site_false_discovery_rate"] is None
PY

# A leakage-safe split keeps every shared TE/event/locus/batch/clade binding in
# one partition and fold. Moving one connected record must fail while still
# emitting a machine-readable violation report.
python3 "$repo/scripts/tevox_split_audit.py" \
    --truth "$data/v05/benchmark/truth.candidates.tsv" \
    --locus-truth "$data/v05/benchmark/truth.members.tsv" \
    --splits "$data/v05/benchmark/splits.valid.tsv" --strict \
    --output "$work/split_valid.json" >/dev/null
if python3 "$repo/scripts/tevox_split_audit.py" \
    --truth "$data/v05/benchmark/truth.candidates.tsv" \
    --locus-truth "$data/v05/benchmark/truth.members.tsv" \
    --splits "$data/v05/benchmark/splits.leaky.tsv" \
    --output "$work/split_leaky.json" >/dev/null 2>&1; then
    echo "accepted a candidate split with cross-partition leakage" >&2
    exit 1
fi
python3 - "$work/split_valid.json" "$work/split_leaky.json" <<'PY'
import json, sys
valid, leaky = (json.load(open(path, encoding="utf-8")) for path in sys.argv[1:])
assert valid["passed"] is True and not valid["violations"]
assert valid["model_status"] == "NO_MODEL_FITTED"
assert valid["locus_truth_provided"] is True
assert valid["counts"]["locus_truth_records"] == 4
assert valid["counts"]["binding_entity_counts"]["homology_group"] == 1
assert leaky["passed"] is False
assert {row["code"] for row in leaky["violations"]} == {
    "BINDING_ENTITY_CROSSES_SPLIT"
}
PY

# Strict auditing must fail closed when locus-member truth is omitted, because
# endpoint WGD/HMG bindings cannot then be checked independently.
if python3 "$repo/scripts/tevox_split_audit.py" \
    --truth "$data/v05/benchmark/truth.candidates.tsv" \
    --splits "$data/v05/benchmark/splits.valid.tsv" --strict \
    --output "$work/split_missing_locus_truth.json" >/dev/null 2>&1; then
    echo "strict split audit passed without locus-member truth" >&2
    exit 1
fi
python3 - "$work/split_missing_locus_truth.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
assert report["passed"] is False
assert report["locus_truth_provided"] is False
assert {row["code"] for row in report["unverifiable"]} == {
    "LOCUS_TRUTH_NOT_PROVIDED"
}
PY
cp "$data/v05/benchmark/truth.candidates.tsv" "$work/split_collision.tsv"
if python3 "$repo/scripts/tevox_split_audit.py" \
    --truth "$work/split_collision.tsv" \
    --locus-truth "$data/v05/benchmark/truth.members.tsv" \
    --splits "$data/v05/benchmark/splits.valid.tsv" \
    --output "$work/split_collision.tsv" >/dev/null 2>&1; then
    echo "split audit allowed its output to overwrite truth" >&2
    exit 1
fi
cmp "$work/split_collision.tsv" "$data/v05/benchmark/truth.candidates.tsv"

# Model preparation exports every generated alignment view from a deliberately
# untruncated run, but only raw pre-decision features. Assessed and unlabelled
# rows remain explicit. It freezes hashes and a generator contract without
# fitting a model or exporting the built-in score.
(
    cd "$data/v05/global"
    "$bin" graph --manifest manifest.tsv --alignments alignments.tsv \
        --flank 20 --candidate-window 10 --max-candidates 0 \
        --max-graph-candidates 0 --output "$work/relative_cwd" >/dev/null
)
python3 "$repo/tests/validate_schema.py" "$work/relative_cwd"
python3 "$repo/scripts/tevox_export_training.py" \
    --prefix "$work/relative_cwd" \
    --truth "$data/v05/benchmark/truth.candidates.tsv" \
    --dataset-id D_GLOBAL --run-id relative-cwd-export \
    --output "$work/relative_cwd_export" >/dev/null
python3 - "$work/relative_cwd.run.json" \
    "$work/relative_cwd_export.dataset.json" <<'PY'
import json, sys
from pathlib import Path
run = json.load(open(sys.argv[1], encoding="utf-8"))
dataset = json.load(open(sys.argv[2], encoding="utf-8"))
assert run["input_files"]
assert all(Path(item["path"]).is_absolute() for item in run["input_files"])
assert all(Path(item["path"]).is_absolute() for item in dataset["input_files"])
assert dataset["commit_marker"]["status"] == \
    "COMMITTED_IF_DATASET_JSON_PRESENT"
PY

python3 "$repo/scripts/tevox_export_training.py" \
    --prefix "$work/global_exact" \
    --truth "$data/v05/benchmark/truth.candidates.tsv" \
    --dataset-id D_GLOBAL --run-id export_fixture \
    --output "$work/training_export" >/dev/null
python3 - "$work/training_export.training.tsv" \
    "$work/training_export.unmatched_truth.tsv" \
    "$work/training_export.dataset.json" <<'PY'
import csv, hashlib, json, sys
training = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
unmatched = list(csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8"), delimiter="\t"))
manifest = json.load(open(sys.argv[3], encoding="utf-8"))
digest = lambda path: hashlib.sha256(open(path, "rb").read()).hexdigest()
assert len(training) == 6
assert len(unmatched) == 2
assert manifest["export_schema_version"] == "1.1.0"
assert manifest["export_status"] == "EVALUATION_ONLY"
assert manifest["model_fit_status"] == "NO_MODEL_FITTED"
assert manifest["sampling"]["downsampled"] is False
assert manifest["sampling"]["scope"] == "ALL_GENERATED_CANDIDATE_VIEWS"
assert manifest["counts"]["exported_rows"] == len(training)
assert manifest["counts"]["exported_truth_status_counts"] == {
    "ASSESSED_LABEL": len(training)
}
assert len(manifest["candidate_generator_sha256"]) == 64
assert manifest["outputs"]["training_tsv_sha256"] == digest(sys.argv[1])
assert manifest["outputs"]["unmatched_truth_tsv_sha256"] == digest(sys.argv[2])
assert manifest["commit_marker"]["protocol"] == "DATASET_JSON_LAST_V1"
for forbidden in (
    "membership_logit", "membership_score", "eligible", "selected",
    "candidate_rank", "locus_id", "solver_component_id",
):
    assert forbidden not in training[0]
assert {row["inclusion_probability"] for row in training} == {"1"}
assert all(row["source_te_semantic_key"].startswith("TEKsha256:") for row in training)
directed = {
    row["source_genome_id"]: (
        row["source_truth_locus_id"], row["target_truth_locus_id"],
        row["taxon_a"], row["taxon_b"],
    )
    for row in training if row["truth_record_id"] == "CANTRUTH001"
}
assert directed == {
    "A": ("TL_GREEDY", "TL_GLOBAL", "TaxonA", "TaxonB"),
    "B": ("TL_GLOBAL", "TL_GREEDY", "TaxonB", "TaxonA"),
}
PY

# Stale/tampered feature sidecars and run-time input hashes must fail closed.
python3 - "$work/global_exact" "$work/export_tampered_feature" \
    "$work/export_stale_hash" <<'PY'
import csv, json, shutil, sys
from pathlib import Path
source, feature_prefix, hash_prefix = map(Path, sys.argv[1:])
for destination in (feature_prefix, hash_prefix):
    for suffix in ("run.json", "evidence.tsv", "candidates.tsv", "candidate_features.tsv"):
        shutil.copyfile(f"{source}.{suffix}", f"{destination}.{suffix}")
feature_path = Path(f"{feature_prefix}.candidate_features.tsv")
with feature_path.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = reader.fieldnames
    rows = list(reader)
row = next(item for item in rows if item["flank_min"] != ".")
row["flank_min"] = "0.123456" if row["flank_min"] != "0.123456" else "0.654321"
with feature_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
run_path = Path(f"{hash_prefix}.run.json")
run = json.load(run_path.open(encoding="utf-8"))
run["genomes"][0]["fasta_sha256"] = "0" * 64
with run_path.open("w", encoding="utf-8") as handle:
    json.dump(run, handle, sort_keys=True)
    handle.write("\n")
PY
for rejected_prefix in export_tampered_feature export_stale_hash; do
    if python3 "$repo/scripts/tevox_export_training.py" \
        --prefix "$work/$rejected_prefix" \
        --truth "$data/v05/benchmark/truth.candidates.tsv" \
        --dataset-id D_GLOBAL --run-id "$rejected_prefix" \
        --output "$work/${rejected_prefix}_output" >/dev/null 2>&1; then
        echo "accepted stale or tampered training-export input: $rejected_prefix" >&2
        exit 1
    fi
done

# Removing truth for one generated semantic pair must not remove its feature
# rows or silently turn label ascertainment into negative sampling.
python3 - "$data/v05/benchmark/truth.candidates.tsv" \
    "$work/truth_partial.tsv" <<'PY'
import csv, sys
with open(sys.argv[1], newline="", encoding="utf-8") as source:
    reader = csv.DictReader(source, delimiter="\t")
    rows = [row for row in reader if row["truth_record_id"] == "CANTRUTH002"]
    fields = reader.fieldnames
with open(sys.argv[2], "w", newline="", encoding="utf-8") as destination:
    writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
PY
python3 "$repo/scripts/tevox_export_training.py" \
    --prefix "$work/global_exact" --truth "$work/truth_partial.tsv" \
    --dataset-id D_GLOBAL --run-id export_partial_fixture \
    --output "$work/training_export_partial" >/dev/null
python3 - "$work/training_export_partial.training.tsv" \
    "$work/training_export_partial.dataset.json" <<'PY'
import csv, json, sys
training = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
manifest = json.load(open(sys.argv[2], encoding="utf-8"))
assert len(training) == 6
assert {row["truth_status"] for row in training} == {
    "ASSESSED_LABEL", "UNLABELLED"
}
assert manifest["counts"]["candidate_rows_without_truth"] > 0
assert manifest["sampling"]["truth_ascertainment_probability"] == \
    "UNKNOWN_NOT_ESTIMATED"
PY

# PARTIAL candidate truth reports prediction-only semantic pairs as unassessed.
# Declaring the same incomplete table EXHAUSTIVE must fail closed.
python3 - "$work/truth_partial.tsv" "$work/global_exact" \
    "$work/candidate_partial.datasets.tsv" \
    "$work/candidate_exhaustive_bad.datasets.tsv" \
    "$work/candidate_scope.runs.tsv" \
    "$data/v05/benchmark/truth.genomes.tsv" <<'PY'
import csv, sys
from pathlib import Path
truth, prefix, partial, exhaustive, runs, genome_truth = map(Path, sys.argv[1:])
fields = [
    "benchmark_schema_version", "dataset_id", "genome_truth",
    "candidate_truth", "locus_truth", "instance_truth",
    "candidate_truth_scope", "candidate_sampling_design",
    "candidate_inclusion_probability", "locus_truth_scope",
]
for path, scope, design, probability in (
    (partial, "PARTIAL", "CASE_CONTROL", "."),
    (exhaustive, "EXHAUSTIVE", "CENSUS", "1"),
):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "benchmark_schema_version": "1.0.0", "dataset_id": "D_GLOBAL",
            "genome_truth": genome_truth.resolve(),
            "candidate_truth": truth.resolve(), "locus_truth": ".",
            "instance_truth": ".", "candidate_truth_scope": scope,
            "candidate_sampling_design": design,
            "candidate_inclusion_probability": probability,
            "locus_truth_scope": "NOT_EVALUATED",
        })
run_fields = [
    "benchmark_schema_version", "run_id", "dataset_id", "method_id",
    "condition_id", "replicate_id", "group_id", "prediction_format", "prefix",
]
with runs.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=run_fields, delimiter="\t")
    writer.writeheader()
    writer.writerow({
        "benchmark_schema_version": "1.0.0", "run_id": "candidate_scope",
        "dataset_id": "D_GLOBAL", "method_id": "TEvoX",
        "condition_id": "scope", "replicate_id": "1", "group_id": "scope",
        "prediction_format": "tevox-1.2", "prefix": prefix.resolve(),
    })
PY
python3 "$repo/scripts/tevox_benchmark.py" \
    --datasets "$work/candidate_partial.datasets.tsv" \
    --runs "$work/candidate_scope.runs.tsv" \
    --output "$work/candidate_partial" >/dev/null
python3 - "$work/candidate_partial.metrics.json" <<'PY'
import json, sys
metrics = json.load(open(sys.argv[1], encoding="utf-8"))["runs"][0]["candidate_membership"]
assert metrics["truth_scope"] == "PARTIAL"
assert metrics["unassessed_predicted_semantic_pairs"] > 0
assert metrics["population_prevalence_or_calibration_claimable"] is False
PY
if python3 "$repo/scripts/tevox_benchmark.py" \
    --datasets "$work/candidate_exhaustive_bad.datasets.tsv" \
    --runs "$work/candidate_scope.runs.tsv" \
    --output "$work/candidate_exhaustive_bad" >/dev/null 2>&1; then
    echo "accepted incomplete candidate truth declared EXHAUSTIVE" >&2
    exit 1
fi
"$bin" graph --manifest "$data/v05/global/manifest.tsv" \
    --alignments "$data/v05/global/alignments.tsv" --flank 20 \
    --candidate-window 10 --max-candidates 0 --max-graph-candidates 1 \
    --output "$work/global_export_truncated" >/dev/null
if python3 "$repo/scripts/tevox_export_training.py" \
    --prefix "$work/global_export_truncated" \
    --truth "$data/v05/benchmark/truth.candidates.tsv" \
    --dataset-id D_GLOBAL --run-id rejected_export \
    --output "$work/training_export_rejected" >/dev/null 2>&1; then
    echo "exported training data from a graph-truncated candidate run" >&2
    exit 1
fi

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
assert a1["graph_candidate_count"] == a0["graph_candidate_count"] == "2"
winner1 = next(row for row in ca1 if row["evidence_id"] == a1["evidence_id"] and row["selected"] == "true")
winner0 = next(row for row in ca0 if row["evidence_id"] == a0["evidence_id"] and row["selected"] == "true")
assert winner1["target_te_id"] == winner0["target_te_id"] == a1["selected_target_te_id"]
assert all(row["graph_retained"] == "true" for row in ca1 + ca0)
PY

# Reporting top-K cannot alter inference. Internal pruning is a separate,
# explicitly recorded control and marks candidates excluded from the graph.
"$bin" pair \
    --genome-a A --fasta-a "$data/conflict/A.fa" --te-a "$data/conflict/A.gff3" \
    --genome-b B --fasta-b "$data/conflict/B.fa" --te-b "$data/v04/topk/B.gff3" \
    --paf "$data/conflict/A_B.paf" --flank 20 --candidate-window 10 \
    --max-candidates 0 --max-graph-candidates 1 \
    --output "$work/topk_graph_1" >/dev/null
python3 "$repo/tests/validate_schema.py" "$work/topk_graph_1"
python3 - "$work/topk_graph_1.evidence.tsv" \
    "$work/topk_graph_1.candidates.tsv" "$work/topk_graph_1.run.json" <<'PY'
import csv, json, sys
evidence = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8"), delimiter="\t"))
candidates = list(csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8"), delimiter="\t"))
run = json.load(open(sys.argv[3], encoding="utf-8"))
row = next(item for item in evidence if item["source_te_id"] == "A_conflict")
rows = [item for item in candidates if item["evidence_id"] == row["evidence_id"]]
assert row["retained_candidate_count"] == "2"
assert row["graph_candidate_count"] == "1"
assert [item["graph_retained"] for item in rows] == ["true", "false"]
assert run["config"]["max_candidates"] == 0
assert run["config"]["max_graph_candidates"] == 1
PY

if "$bin" pair \
    --genome-a A --fasta-a "$data/conflict/A.fa" --te-a "$data/conflict/A.gff3" \
    --genome-b B --fasta-b "$data/conflict/B.fa" --te-b "$data/conflict/B.gff3" \
    --paf "$data/conflict/A_B.paf" --min-membership 0.49 \
    --output "$work/low_membership" >/dev/null 2>&1; then
    echo "accepted a negative-log-odds membership threshold" >&2
    exit 1
fi

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
    "$repo/scripts/tevox_phylo.py" "$repo/scripts/tevox_score_audit.py" \
    "$repo/scripts/tevox_benchmark.py" \
    "$repo/scripts/tevox_export_training.py" \
    "$repo/scripts/tevox_split_audit.py"
echo 'All TEvoX tests passed.'
