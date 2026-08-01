# TEvoX

TEvoX reconstructs orthologous transposable-element (TE) loci across assembled
genomes while keeping biological state, technical callability and annotation
status separate.

Version `0.3.0-alpha.1` is a correctness and evidence-contract release. It fixes
PAF direction documentation, computes local identity only from `cs:Z` or exact
`=X` CIGAR operations, treats PAF MAPQ 255 as missing, prevents synthetic
reverse views from masquerading as independent reciprocal evidence, rejects
zero-overlap adjacent TE candidates, and makes contig-edge flanks missing rather
than perfect. Every result is traceable through schema `1.0.0` from alignment
record to directed observation, candidate, decision and locus instance.

> Alpha status: the inference rules and schema have regression coverage, but the
> quality score is an uncalibrated ranking score. Large-scale simulation and
> biological benchmarks are still required before manuscript-level accuracy
> claims.

## Build and test

```bash
make
make check
make asan
```

The C core has no runtime library dependencies. Python 3 is used for schema
validation tests and exploratory phylogenetic inference.

## Pairwise use and PAF direction

PAF columns 1–4 describe the **query** and columns 6–9 describe the **target**.
TEvoX pair mode requires genome A to be the PAF query and genome B to be the PAF
target. Because minimap2 takes the target/reference first, the correct order is:

```bash
minimap2 -cx asm5 --cs=long B.fa A.fa > A_query_B_target.paf

./tevox pair \
  --genome-a A --fasta-a A.fa --te-a A.te.gff3 \
  --genome-b B --fasta-b B.fa --te-b B.te.gff3 \
  --paf A_query_B_target.paf --output results/A_B
```

`cg:Z` is mandatory. `cs:Z` is recommended. If `cs:Z` is absent, local identity
is available only when `cg:Z` uses exact `=`/`X` operations throughout; an
ordinary `M` CIGAR produces `local_identity=.` rather than substituting the
whole-PAF identity. MAPQ 255 is likewise emitted as missing and cannot pass the
MAPQ gate.

Only one PAF direction is needed. TEvoX creates a reverse traversal with the
same `evidence_group_id`. Listing a separately generated opposite-direction PAF
is allowed; a reciprocal bonus is possible only when the two supporting views
belong to distinct evidence groups.

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
reference	sample_a	reference_query_sample_a_target.paf
reference	polyploid_b	reference_query_polyploid_b_target.paf
sample_a	polyploid_b	sample_a_query_polyploid_b_target.paf
```

```bash
./tevox graph --manifest genomes.tsv --alignments alignments.tsv \
  --output results/cohort
```

Relative paths resolve from their table's directory. Genome and TE records are
sorted internally so manifest, annotation and PAF record order do not change
the inferred loci or stable evidence IDs.

`max_locus_copies` is retained only as a legacy graph-capacity guard. It is not
a WGD model and does not establish homeology. In particular, adjacent TEs with
zero reciprocal overlap remain separate regardless of this value. Explicit
subgenome/copy contexts are planned for v0.4.

## Evidence schema and outputs

For prefix `cohort`, TEvoX writes:

| File | Grain | Purpose |
|---|---|---|
| `cohort.evidence.tsv` | directed alignment observation | raw/derived provenance and measured features |
| `cohort.candidates.tsv` | observation × target TE | every nearby candidate, eligibility and score |
| `cohort.decisions.tsv` | source TE × target genome | aggregation of all observations, including every near-best ambiguity trigger |
| `cohort.edges.tsv` | TE pair | globally evaluated graph support and rejection reason |
| `cohort.loci.tsv` | locus | selected graph components |
| `cohort.instances.tsv` | locus × genome | three-axis state, claimability and supporting external keys |
| `cohort.states.tsv` | locus × genome | backward-compatible legacy state matrix plus v0.3 fields |
| `cohort.summary.tsv` | genome × legacy state | counts |
| `cohort.run.json` | run | version, schema, parameters, inputs and row counts |

Schema `1.0.0` uses zero-based half-open coordinates and `.` for an unobserved
value. Missing is never encoded as zero. See [schema](docs/SCHEMA.md).

## State and claim model

Each observation and instance carries three axes:

- technical: `CALLABLE`, `GAP`, `AMBIGUOUS`, `UNCALLABLE`;
- biological: `PRESENT`, `EMPTY`, `STRUCTURAL_ALTERNATIVE`, `UNKNOWN`;
- annotation: `MATCHED`, `MISSING`, `FAMILY_CONFLICT`, `NOT_APPLICABLE`,
  `UNKNOWN`.

`claimable` is a separate gate. For example, an insertion-like geometry with
missing local identity may have biological state `EMPTY`, but it is not emitted
as legacy `EMPTY_SITE_CONFIRMED` and cannot support an absence claim. A
confirmed empty site requires full paired flanks, observed passing MAPQ,
callable target sequence, exact local identity and at least 70% of the source TE
represented as a query insertion.

The compatibility state vocabulary remains:

| State | Interpretation |
|---|---|
| `PRESENT_ANNOTATED` | compatible annotation assigned to the locus |
| `PRESENT_UNANNOTATED` | claimable aligned TE sequence lacks annotation |
| `EMPTY_SITE_CONFIRMED` | claimable paired-flank empty site |
| `STRUCTURAL_ALTERNATIVE` | claimable non-insertion structural difference |
| `FAMILY_OR_BOUNDARY_DISCORDANCE` | callable annotation conflict |
| `ASSEMBLY_GAP` | excessive target `N`/gap sequence |
| `PROJECTION_AMBIGUOUS` | candidates, projections or global loci disagree |
| `UNCALLABLE` | required evidence is absent or below threshold |

## Phylogenetic candidates

```bash
python3 scripts/tevox_phylo.py --states results/cohort.states.tsv \
  --tree species_tree.nwk --output results/cohort_phylo
```

Only annotated/unannotated presence is encoded as present and only a confirmed
empty site as absent; technical states are missing data. Event placements are
exploratory parsimony candidates, not proof of ancestral state.

See [design](docs/DESIGN.md), [input contracts](docs/INPUTS.md),
[validation](docs/VALIDATION.md) and [roadmap](docs/ROADMAP.md). Licensed under
Apache-2.0.
