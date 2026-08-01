# TEvoX

TEvoX reconstructs transposable-element (TE) loci across assembled genomes
while keeping technical callability, biological state, annotation status and
evidence provenance separate.

Version `0.5.0-alpha.1` adds an auditable inference and optimization layer on
top of the v0.4 evidence backends and copy contexts:

- explicit missing-aware candidate and three-axis observation feature tables;
- a versioned `BUILTIN_UNCALIBRATED_V1` membership score;
- exact Hungarian one-to-one matching within MCScanX copy-context pairs;
- exact constrained component optimization for small graphs and a labelled
  deterministic fallback for larger graphs;
- `ORTHOLOG`, `WGD_HOMEOLOG`, `ALLELIC`, tandem/segmental/transposed paralog
  and `UNKNOWN` relation score columns;
- prediction sets, entropy, out-of-domain flags, solver objective/bound/gap;
- evidence schema `1.2.0` with stable edge, relation, matching and solver keys.

> Alpha status: the state rules and schema have regression coverage, but the
> built-in inference scores are explicitly uncalibrated and are not posterior
> probabilities. The discrete state/claimability safety gates remain
> authoritative. Publication-level calibration, simulation and biological
> benchmarks remain future work.

## Build and test

```bash
make
make check
make asan
```

The C core has no runtime library dependencies. Python 3 is used by the schema
validator and exploratory phylogeny helper.

## Pairwise DNA alignment

TEvoX pair mode defines genome A as alignment **query** and genome B as
alignment **target/reference**.

### minimap2 / PAF

minimap2 takes the target first, so use:

```bash
minimap2 -cx asm5 --cs=long B.fa A.fa > A_query_B_target.paf

./tevox pair \
  --genome-a A --fasta-a A.fa --te-a A.te.gff3 \
  --genome-b B --fasta-b B.fa --te-b B.te.gff3 \
  --paf A_query_B_target.paf --output results/A_B
```

`cg:Z` is mandatory and `cs:Z` is recommended. Local identity comes only from
`cs:Z` or an exact `=`/`X` CIGAR; ordinary `M` makes `local_identity=.`. PAF
MAPQ 255 is missing, not a high-confidence value.

### MUMmer4 / NUCMER delta

`nucmer` likewise takes the reference/target first:

```bash
nucmer -p A_query_B_target B.fa A.fa

./tevox pair \
  --genome-a A --fasta-a A.fa --te-a A.te.gff3 \
  --genome-b B --fasta-b B.fa --te-b B.te.gff3 \
  --alignment A_query_B_target.delta --alignment-format delta \
  --output results/A_B_delta
```

TEvoX accepts strict `NUCMER` delta only; `PROMER` is rejected. Delta error
counts produce `alignment_identity` with method `DELTA_ERROR_COUNT`, but delta
does not identify individual substitutions. Consequently `local_identity` and
MAPQ are missing. An unambiguous passing delta alignment can establish
`UNIQUE_ALIGNMENT` for annotated presence/discordance, but cannot by itself
produce `EMPTY_SITE_CONFIRMED`, unannotated sequence presence or a structural
alternative claim.

The delta header must list the target/reference FASTA first and query FASTA
second. TEvoX compares canonical paths when they resolve; otherwise it permits
basename matching only when the two expected FASTA basenames are distinct.
Unresolved equal basenames are rejected because their direction cannot be
verified. Using absolute FASTA paths in the `nucmer` command gives the strongest
validation.

One native alignment direction is sufficient. TEvoX creates a dependent
reverse traversal with the same `evidence_group_id`; only separately generated
opposite native records can earn independent reciprocal support.

## Multi-genome alignment and MCScanX context

`genomes.tsv`:

```text
genome_id	fasta	te_annotation	max_locus_copies
A	A.fa	A.te.gff3	1
B	B.fa	B.te.gff3	1
```

`alignments.tsv` supports the v0.4 four-column form (the v0.3 three-column PAF
form remains accepted):

```text
query_id	target_id	format	path
A	B	paf	A_query_B_target.paf
A	B	delta	A_query_B_target.delta
```

An MCScanX source uses a normalized gene table plus native `.collinearity`:

```text
source_id	format	collinearity	genes	wgd_node
wgd_layer_1	mcscanx	cohort.collinearity	genes.tsv	WGD1
```

```text
genome_id	gene_id	contig	start	end	strand	subgenome_id	haplotype_id
A	A_gene1	chr1	100	200	+	A	hap1
B	B_gene1	chr1	110	210	+	B1	hap1
```

Gene coordinates are explicitly 0-based half-open. MCScanX emits unqualified
gene IDs, so v0.4 requires gene IDs to be globally unique; prefix them by
genome before running MCScanX.

```bash
./tevox graph --manifest genomes.tsv \
  --alignments alignments.tsv \
  --synteny synteny.sources.tsv \
  --output results/cohort
```

MCScanX supplies only a locus/copy-context prior. It never supplies base-level
presence or absence and MCScanX-only runs cannot emit
`EMPTY_SITE_CONFIRMED`. A candidate relation is strong only when both TEs are
bracketed by consecutive anchors in `PASS` contexts. Metadata conflicts become
`AMBIGUOUS`; block-interior and missing contexts remain neutral.

`max_locus_copies` remains a fallback for TEs without a strong context. Strong
contexts replace it with one slot per stable `context_id`; different WGD copy
contexts may coexist, while two members of the same context cannot be merged.
`copy001` labels are display order only and must not be used as join keys.

Current context consolidation merges near-duplicate block sides with at least
80% reciprocal span overlap and 50% anchor Jaccard. It does not yet stitch all
adjacent fragmented MCScanX blocks, so fragmented inputs can over-segment copy
contexts. Inspect `contexts.tsv` before biological interpretation.

## v0.5 inference

Candidate scores preserve missingness, then block-level matching prevents
multiple TEs from competing for one explicit copy slot. The global optimizer
maximizes positive membership log odds subject to family, copy-context, HMG and
fallback-quota constraints. Components with at most 18 active edges are solved
exactly by default; larger components report a heuristic result and a loose
but valid upper-bound gap.

Useful controls are:

```text
--min-membership 0.50
--prediction-mass 0.90
--max-candidates 64
--max-graph-candidates 64
--exact-max-edges 18
--exact-match-nodes 256
--tandem-distance 10000
```

See the exact coefficient, matching, optimization and non-claim contract in
[v0.5 inference](docs/V05_INFERENCE.md).

Given an independent truth TSV with `candidate_id`, binary `label` and a
species-pair/clade `group_id`, the audit helper reports Brier score, log loss,
ECE, AUROC, calibration bins and group-wise held-out evaluation:

```bash
python3 scripts/tevox_score_audit.py \
  --features cohort.candidate_features.tsv \
  --truth independent_truth.tsv --output cohort.audit
```

The audit deliberately retains `EVALUATION_ONLY_NOT_A_CALIBRATED_MODEL`; it
does not fit or bless a calibrated model.

## Evidence schema and outputs

For prefix `cohort`, schema `1.2.0` writes:

| File | Grain | Purpose |
|---|---|---|
| `cohort.evidence.tsv` | directed alignment observation | PAF/delta provenance, measured and missing features |
| `cohort.observation_scores.tsv` | directed observation | three-axis normalized scores, prediction sets, entropy and OOD |
| `cohort.candidates.tsv` | observation × target TE | reported top-K candidates, graph-retention flag, rank, eligibility and context relation |
| `cohort.candidate_features.tsv` | candidate | observed mask, missing features and membership score |
| `cohort.candidate_contexts.tsv` | candidate | source/target context foreign keys |
| `cohort.decisions.tsv` | source TE × target genome | aggregation and every near-best ambiguity trigger |
| `cohort.edges.tsv` | TE pair | membership, block matching, solver key and selection reason |
| `cohort.relations.tsv` | TE pair | relation score vector, prediction set and locus/edge keys |
| `cohort.solver.tsv` | candidate component | method, status, objective, upper bound and gap |
| `cohort.loci.tsv` | locus | selected graph component |
| `cohort.instances.tsv` | locus × genome | three-axis state, claimability and evidence keys |
| `cohort.synteny.blocks.tsv` | MCScanX block | provider record, block/context and homology-group keys |
| `cohort.synteny.anchors.tsv` | anchor pair | normalized and provider ranks plus E-value |
| `cohort.contexts.tsv` | copy context | genome/subgenome/haplotype/copy/WGD tuple |
| `cohort.te_contexts.tsv` | TE × context | bracketed/interior/ambiguous assignment |
| `cohort.states.tsv` | locus × genome | backward-compatible eight-state view |
| `cohort.summary.tsv` | genome × state | counts |
| `cohort.run.json` | run | versions, parameters, providers, counts and performance counters |

All coordinates are 0-based half-open. A literal `.` means unobserved or not
applicable and never numeric zero. See [schema](docs/SCHEMA.md),
[input contracts](docs/INPUTS.md), [design](docs/DESIGN.md),
[v0.5 inference](docs/V05_INFERENCE.md),
[validation](docs/VALIDATION.md) and [roadmap](docs/ROADMAP.md).

## State and phylogeny

The internal axes are:

- technical: `CALLABLE`, `GAP`, `AMBIGUOUS`, `UNCALLABLE`;
- biological: `PRESENT`, `EMPTY`, `STRUCTURAL_ALTERNATIVE`, `UNKNOWN`;
- annotation: `MATCHED`, `MISSING`, `FAMILY_CONFLICT`, `NOT_APPLICABLE`,
  `UNKNOWN`.

`claimable` is evaluated separately. Missing evidence cannot be converted into
absence. The exploratory phylogeny helper treats technical/nonclaimable states
as missing:

```bash
python3 scripts/tevox_phylo.py --states results/cohort.states.tsv \
  --tree species_tree.nwk --output results/cohort_phylo
```

Licensed under Apache-2.0.
