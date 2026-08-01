# TEvoX 0.4 inference design

This document freezes `0.4.0-alpha.1` and schema `1.1.0`. It describes an
evidence-aware deterministic model, not the calibrated probabilistic model
planned for v0.5.

## Evidence layers

TEvoX keeps two independent provider classes:

1. base alignment evidence from PAF or MUMmer4 NUCMER delta;
2. gene-synteny context from MCScanX.

MCScanX narrows and labels locus context but cannot observe TE sequence or an
empty insertion site. Base alignment controls projection and biological state.
The two layers meet only at candidate context relation and graph constraints.

One native alignment record defines an evidence group. Its automatic reverse
traversal has the same group and is dependent. MCScanX blocks receive their own
provider groups and stable block/anchor IDs but never become alignment
observations.

## Alignment normalization

PAF local identity is computed from `cs` or exact `=X`, separately from
`matches/block_len`. MAPQ 255 is missing.

NUCMER delta is converted from 1-based closed reference/query coordinates to a
0-based half-open CIGAR-like path. Positive deltas are reference-only `D` and
negative deltas query-only `I`. Error counts define aggregate identity, while
local identity and MAPQ remain missing. A delta mapping becomes
`UNIQUE_ALIGNMENT` only after source-level comparison finds one unambiguous
passing projection.

Both providers generate the same projection features: paired flank coverage,
TE-body aligned fraction, query insertion fraction, projected coordinates,
target gap fraction and nearby annotation candidates. Provider-unobserved
features remain missing.

## Candidate generation and top-K

Alignments and target TEs are indexed by genome pair/contig and interval.
Candidate-window arithmetic is clamped to contig bounds. Every exact-window
candidate contributes to eligibility, best/second score and ambiguity before
optional output truncation.

The same complete order selects the winner, retains bounded top-K and assigns
rank. This prevents a score-tied winner from being discarded. `max_candidates=0`
appends all candidates and sorts once; bounded mode grows memory incrementally.

A candidate requires minimum reciprocal interval overlap. Zero-overlap adjacent
TEs cannot become loci through family or copy quota. A known MCScanX conflict
marks the projection ambiguous; missing context remains neutral.

## MCScanX contexts

The adapter validates alignment IDs/ranks, anchor count, contigs, orientation,
genomes and gene-table coordinates. Side order and anchors are normalized for
stable IDs.

Near-duplicate sides merge only when they share genome/contig and compatible
WGD/subgenome/haplotype metadata, overlap each other by at least 80% in both
directions and have anchor Jaccard ≥0.50. Opposing sides of the same block and
different WGD layers cannot merge.

Each resulting context records:

```text
(genome, contig interval, subgenome, haplotype,
 syntenic display copy, WGD node, anchor set)
```

Blocks connect contexts into a homology group. Within an HMG, contexts are
canonically ordered per genome to assign display labels `copy001...`; the
stable join key is `context_id`.

A TE is strong only when exactly one `PASS` context contains consecutive left
and right block anchors that fully bracket it. Multiple brackets or metadata
conflict are ambiguous; simple block overlap is `BLOCK_INTERIOR` and neutral.

The alpha deliberately does not infer adjacent-fragment stitching for every
MCScanX output. This avoids silently joining distinct local duplications, but
may over-segment fragmented collinearity blocks. A validated stitching model is
future work.

## State and claimability

Observations retain three axes:

```text
Technical:   CALLABLE | GAP | AMBIGUOUS | UNCALLABLE
Biological:  PRESENT | EMPTY | STRUCTURAL_ALTERNATIVE | UNKNOWN
Annotation:  MATCHED | MISSING | FAMILY_CONFLICT | NOT_APPLICABLE | UNKNOWN
```

Claimability is evaluated after these axes. Technical gap/ambiguity blocks all
claims. Annotated presence can use an established PAF or unique delta mapping.
Unannotated presence, empty site and structural alternative additionally need
observed local identity, which delta does not supply. Therefore aggregate delta
identity cannot leak into an absence claim.

The eight-state compatibility field is a projection of axes plus claimability;
for example `BIO=EMPTY, claimable=false` becomes `UNCALLABLE`, never
`EMPTY_SITE_CONFIRMED`.

## Global locus graph

Only callable, unambiguous decisions with an eligible selected annotation can
support an edge. Evidence support is sorted and aggregated once per TE pair;
native reciprocal directions must use distinct groups.

Edges are processed by deterministic `score, node A, node B` order. A merge
must:

1. pass the edge threshold;
2. have no direct or component-wide known-family conflict;
3. avoid two members in one strong `context_id`;
4. avoid mixing incompatible known HMGs;
5. satisfy the legacy per-genome fallback quota only for nodes without a strong
   context.

Thus two WGD contexts in one HMG can remain co-orthologous copies, while two TEs
bracketed by the same context are not treated as WGD copies. This is a discrete
constraint model; ORTHOLOG/WGD_HOMEOLOG posterior classification and ILP
clustering are deferred to v0.5.

## Performance model

- FASTA gap runs are indexed once; inference records zero post-index reopens.
- PAF/delta and TE interval envelopes use sorted prefix-max indexes.
- MCScanX gene IDs use sorted lookup; block sides use chromosome interval
  sweep rather than all-pairs comparison.
- TE-context assignments are stored in per-node contiguous ranges.
- raw TE edges and evidence support are sorted/reduced; reciprocal detection is
  linear in support records.
- graph edge ordering uses `qsort`.

Component family/quota checks and TE-to-context construction still have
scaling work before chromosome-scale benchmark claims. Runtime and memory must
be measured on publication datasets rather than inferred from unit tests.

## Phylogeny

The helper maps claimable presence to 1, confirmed empty to 0 and every
technical/nonclaimable state to missing before Fitch parsimony. It remains an
exploratory event-candidate module; stochastic Dollo inference is planned only
after locus accuracy is calibrated.
