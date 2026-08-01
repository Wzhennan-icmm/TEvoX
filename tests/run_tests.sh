#!/usr/bin/env bash
set -euo pipefail
repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd);bin="$repo/tevox";d="$repo/tests/data";w=$(mktemp -d "${TMPDIR:-/tmp}/tevox.XXXXXX");trap 'rm -rf "$w"' EXIT
has(){ local p=${2//\\t/$'\t'};grep -Eq "$p" "$1"||{ echo "missing $2" >&2;exit 1;};}
"$bin" --version|grep -F '0.2.0-alpha.1' >/dev/null
"$bin" pair --genome-a A --fasta-a "$d/pair/A.fa" --te-a "$d/pair/A.gff3" --genome-b B --fasta-b "$d/pair/B.fa" --te-b "$d/pair/B.gff3" --paf "$d/pair/A_B.paf" --flank 20 --candidate-window 10 -o "$w/pair"
has "$w/pair.states.tsv" '^TEL[0-9]+\tB\tEMPTY_SITE_CONFIRMED\t';has "$w/pair.edges.tsv" '^A_shared\tA\tB_shared\tB\t.*\ttrue\ttrue\t.*\ttrue$';diff -u "$repo/tests/golden/pair.states.tsv" "$w/pair.states.tsv"
"$bin" pair --genome-a A --fasta-a "$d/pair/A.fa" --te-a "$d/pair/A.gff3" --genome-b B --fasta-b "$d/pair/B.fa" --te-b "$d/unannotated.gff3" --paf "$d/pair/A_B.paf" --flank 20 --candidate-window 10 -o "$w/unannotated";has "$w/unannotated.states.tsv" '^TEL[0-9]+\tB\tPRESENT_UNANNOTATED\t'
"$bin" pair --genome-a A --fasta-a "$d/conflict/A.fa" --te-a "$d/conflict/A.gff3" --genome-b B --fasta-b "$d/conflict/B.fa" --te-b "$d/conflict/B.gff3" --paf "$d/conflict/A_B.paf" --flank 20 --candidate-window 10 -o "$w/conflict"
has "$w/conflict.states.tsv" '^TEL[0-9]+\tB\tFAMILY_OR_BOUNDARY_DISCORDANCE\t';! grep -Eq $'\ttrue$' "$w/conflict.edges.tsv"
"$bin" pair --genome-a A --fasta-a "$d/conflict/A.fa" --te-a "$d/conflict/A.gff3" --genome-b B --fasta-b "$d/gap/B.fa" --te-b "$d/conflict/B.gff3" --paf "$d/conflict/A_B.paf" --flank 20 --candidate-window 10 -o "$w/gap";has "$w/gap.states.tsv" '^TEL[0-9]+\tB\tASSEMBLY_GAP\t'
"$bin" pair --genome-a A --fasta-a "$d/conflict/A.fa" --te-a "$d/conflict/A.gff3" --genome-b B --fasta-b "$d/conflict/B.fa" --te-b "$d/reverse/B.gff3" --paf "$d/reverse/A_B.paf" --flank 20 --candidate-window 10 -o "$w/reverse";has "$w/reverse.edges.tsv" '^A_conflict\tA\tB_reverse\tB\t.*\ttrue\ttrue\t.*\ttrue$'
"$bin" graph --manifest "$d/multi/manifest.tsv" --alignments "$d/multi/alignments.tsv" --flank 20 --candidate-window 100 -o "$w/multi";has "$w/multi.states.tsv" '^TEL[0-9]+\tG2\tPRESENT_ANNOTATED\t.*\t2\t'
"$bin" graph --manifest "$d/multi/manifest1.tsv" --alignments "$d/multi/alignments.tsv" --flank 20 --candidate-window 100 -o "$w/one";! grep -Eq $'^TEL[0-9]+\tG2\tPRESENT_ANNOTATED\t.*\t2\t' "$w/one.states.tsv"
python3 "$repo/scripts/tevox_phylo.py" --states "$d/multi/phylo_states.tsv" --tree "$d/multi/tree.nwk" -o "$w/phylo";has "$w/phylo.events.tsv" '^TEL000001\t(gain|loss)\t'
if "$bin" pair --genome-a A --fasta-a "$d/conflict/A.fa" --te-a "$d/conflict/A.gff3" --genome-b B --fasta-b "$d/conflict/B.fa" --te-b "$d/conflict/B.gff3" --paf "$d/invalid.paf" -o "$w/bad" >/dev/null 2>&1;then echo 'accepted invalid PAF' >&2;exit 1;fi
echo 'All TEvoX tests passed.'
