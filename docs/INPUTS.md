# Input contracts

All validation failures are fatal. TEvoX does not silently clip coordinates or
guess alignment direction. Relative paths resolve from the table that contains
them.

## FASTA and TE annotation

- FASTA contig IDs must be unique. Sequence lines must have fixed width except
  the last line of each contig. Contig lengths are authoritative.
- GFF3 TE coordinates are converted from 1-based closed to 0-based half-open.
  `ID`/`Name` identifies a TE; `family`/`Family`/`classification` supplies its
  family.
- BED uses 0-based half-open `chrom start end id score strand family type`;
  only the first three columns are mandatory.
- TE IDs must be unique within a genome. Records are canonically sorted.
- FASTA `N`, `n`, `-` and `.` runs are indexed once as technical gap sequence;
  inference does not repeatedly reopen FASTA files.

## Genome manifest

```text
genome_id	fasta	te_annotation	max_locus_copies
A	A.fa	A.te.gff3	1
B	B.fa	B.te.gff3	1
```

The positive `max_locus_copies` field is a fallback only for TEs without a
strong copy-context assignment. It is not WGD evidence.

## Alignment table

Preferred schema:

```text
query_id	target_id	format	path
A	B	paf	A_query_B_target.paf
A	B	delta	A_query_B_target.delta
```

Formats are case-insensitive: `paf`, `delta`, `mummer-delta` and
`mummer_delta`. The legacy three-column `query_id target_id paf` table remains
accepted and implies PAF.

Each native record is normalized into one `evidence_group_id`; its automatic
reverse traversal is dependent evidence in the same group. Listing an opposite
direction is appropriate only when it was generated independently.

## PAF

PAF columns 1–4 are query and columns 6–9 target. Pair mode requires A=query
and B=target. Since minimap2 takes target first:

```bash
minimap2 -cx asm5 --cs=long B.fa A.fa > A_query_B_target.paf
```

Contract:

- 12 core fields and `cg:Z` are mandatory;
- `M`, `=`, `X`, `I` and `D` are supported;
- CIGAR consumption must equal declared query/target spans;
- `cs:Z` supports `:`, `=`, `*`, `+` and `-`; splice `~` is rejected;
- `cg` and `cs` must describe the same match/indel path;
- local identity comes from `cs` or a completely exact `=X` CIGAR;
- ordinary `M` gives `local_identity=.`; whole-record identity is retained
  separately as `alignment_identity` with method `PAF_CORE`;
- MAPQ 255 is emitted as `mapq=.` and `mapq_status=MISSING_255`;
- contig IDs, lengths, spans and declared direction must match FASTA.

## MUMmer4 delta

For A=query and B=target/reference:

```bash
nucmer -p A_query_B_target B.fa A.fa
```

The adapter reads native `.delta`, not `show-coords`:

- the file type must be exactly `NUCMER`; `PROMER` is rejected;
- reference is TEvoX target and query is TEvoX query. The delta file-header
  paths are checked in that order against the loaded FASTA paths;
- when both a header token and its expected FASTA path resolve, their canonical
  paths must be identical. If either token cannot be resolved, TEvoX accepts
  basename fallback only when the expected target/query basenames are distinct
  and match the header respectively. Equal unresolved basenames are rejected as
  directionally ambiguous;
- 1-based closed delta coordinates are converted to 0-based half-open;
- positive deltas consume reference-only bases (`D`), negative deltas consume
  query-only bases (`I`);
- reverse-query paths retain `I`/`D` meaning while operation order is reversed;
- path consumption, terminator zero, error counts, FASTA names/lengths and
  coordinate bounds are validated;
- aggregate identity is `(alignment columns - error count) / alignment
  columns` and is labelled `DELTA_ERROR_COUNT`;
- delta supplies neither local substitution placement nor MAPQ, so
  `local_identity=.` and `mapq_status=NOT_PROVIDED`.

`--min-delta-identity` gates the aggregate delta alignment. Even a passing
delta cannot establish an empty-site or other sequence claim that requires
local identity.

## MCScanX synteny sources

Sources table:

```text
source_id	format	collinearity	genes	wgd_node
wgd_layer_1	mcscanx	cohort.collinearity	genes.tsv	WGD1
```

Normalized gene table:

```text
genome_id	gene_id	contig	start	end	strand	subgenome_id	haplotype_id
A	A_gene1	chr1	100	200	+	A	hap1
B	B_gene1	chr1	110	210	+	B1	hap1
```

This gene table is a TEvoX adapter format, not MCScanX native GFF. Coordinates
must be 0-based half-open; convert any 1-based source before import. Use `.` for
unknown subgenome/haplotype/WGD metadata.

MCScanX `.collinearity` requirements:

- header form `## Alignment ID: score=S e_value=E N=N contigA&contigB plus|minus`;
- at least two anchors;
- anchor prefix `ID-rank:` must match the current alignment;
- ranks must uniquely cover `0..N-1`;
- the third anchor field is an E-value, not a score;
- each side must remain in one declared genome/contig and obey orientation;
- alignment IDs must be unique within a source;
- because MCScanX anchor rows contain bare IDs, gene IDs must be globally
  unique across all genomes (genome prefixes are recommended).

The adapter canonicalizes A/B side order for stable biological IDs. Different
WGD nodes and incompatible known subgenome/haplotype values are not merged.
Only near-duplicate overlapping block sides are consolidated in this alpha;
general adjacent-fragment stitching is not yet implemented.

## Phylogeny input

Newick tree tip names must equal state-matrix genome IDs. Quoted labels and
Newick comments are not supported in this alpha.
