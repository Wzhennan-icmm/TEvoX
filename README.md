# TEvoX

TEvoX reconstructs orthologous transposable-element (TE) loci across assembled
genomes while keeping biological absence, annotation discordance, assembly gaps
and ambiguous mappings separate.

The `0.2.0-alpha` rewrite replaces the original rule that treated every TE
outside a synteny block as genome-specific. It projects TE boundaries and
paired flanks through base-level PAF CIGAR alignments, scores annotation-aware
candidates and selects a globally consistent multi-genome locus graph.

> Alpha status: the locus model and output schema are tested and usable, but
> the 0–100 quality score is an interpretable heuristic, not a calibrated
> probability. Biological benchmarks remain necessary before manuscript claims.

## Features

- strict FASTA, GFF3/BED, PAF length and `cg:Z` CIGAR validation;
- paired-flank, local-identity, mapping-quality and target-gap evidence;
- joint two- or multi-genome locus reconstruction;
- per-genome WGD/co-ortholog copy quotas;
- eight uncertainty-aware locus states;
- deterministic edge, locus, state and summary TSV outputs;
- exploratory Fitch-parsimony gain/loss inference on Newick trees.

## Build and test

```bash
make
make check
make asan
```

The C core has no runtime library dependencies. Python 3 is needed only for
phylogenetic inference and its regression test.

## Pairwise use

```bash
minimap2 -cx asm5 --cs=long A.fa B.fa > A_B.paf
./tevox pair \
  --genome-a A --fasta-a A.fa --te-a A.te.gff3 \
  --genome-b B --fasta-b B.fa --te-b B.te.gff3 \
  --paf A_B.paf --output results/A_B
```

Only one PAF direction is required; TEvoX derives the reverse CIGAR traversal.

## Multi-genome use

`genomes.tsv`:

```text
genome_id	fasta	te_annotation	max_locus_copies
reference	reference.fa	reference.te.gff3	1
sample_a	sample_a.fa	sample_a.te.gff3	1
polyploid_b	polyploid_b.fa	polyploid_b.te.gff3	2
```

`alignments.tsv`:

```text
query_id	target_id	paf
reference	sample_a	reference_sample_a.paf
reference	polyploid_b	reference_polyploid_b.paf
sample_a	polyploid_b	sample_a_polyploid_b.paf
```

```bash
./tevox graph --manifest genomes.tsv --alignments alignments.tsv \
  --output results/cohort
```

Relative paths resolve from their table's directory. Pairwise alignments need
not form a complete clique, but disconnected genomes cannot share components.

## State vocabulary

| State | Interpretation |
|---|---|
| `PRESENT_ANNOTATED` | compatible annotation assigned to the locus |
| `PRESENT_UNANNOTATED` | aligned TE-length sequence lacks annotation |
| `EMPTY_SITE_CONFIRMED` | paired flanks join across a query TE insertion |
| `STRUCTURAL_ALTERNATIVE` | non-insertion structural difference |
| `FAMILY_OR_BOUNDARY_DISCORDANCE` | annotation family/locus conflict |
| `ASSEMBLY_GAP` | excessive target `N`/gap sequence |
| `PROJECTION_AMBIGUOUS` | competing mappings/candidates disagree |
| `UNCALLABLE` | mapping or flank evidence is insufficient |

Coordinates are zero-based, half-open. For prefix `cohort`, outputs are
`cohort.edges.tsv`, `cohort.loci.tsv`, `cohort.states.tsv` and
`cohort.summary.tsv`.

## Phylogenetic candidates

```bash
python3 scripts/tevox_phylo.py --states results/cohort.states.tsv \
  --tree species_tree.nwk --output results/cohort_phylo
```

Only annotated/unannotated presence is encoded as present and only a confirmed
empty site as absent; all technical states are missing data. Event placements
are exploratory parsimony candidates, not proof of ancestral state.

See [design](docs/DESIGN.md), [input contracts](docs/INPUTS.md),
[validation](docs/VALIDATION.md) and [roadmap](docs/ROADMAP.md). Licensed under
Apache-2.0.
