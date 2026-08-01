# Input contracts

All validation failures are fatal; TEvoX does not silently clip inconsistent
coordinates.

## FASTA and annotation

- **FASTA:** plain text, unique contig IDs and fixed-width sequence lines except
  the final line. Contig lengths are authoritative. FASTA is indexed in place.
- **GFF3:** paths containing `.gff` are parsed as GFF3. Supply a TE-only file.
  `ID`/`Name` identifies a TE and `family`/`Family`/`classification` supplies
  family. Coordinates are converted from 1-based closed to 0-based half-open.
- **BED:** `chrom start end id score strand family type`; the first three
  0-based half-open columns are required.
- TE IDs must be unique within a genome. Annotation records are sorted by
  contig, start, end and ID before analysis.

## PAF

PAF columns 1–4 are query and columns 6–9 are target. The IDs declared in pair
mode or `alignments.tsv` must agree with that direction. Contig names and PAF
lengths are checked against the corresponding FASTA.

The recommended pair-mode command for A=query and B=target is:

```bash
minimap2 -cx asm5 --cs=long B.fa A.fa > A_query_B_target.paf
```

Requirements and semantics:

- the 12 core fields and `cg:Z` are mandatory;
- supported CIGAR operations are `M`, `=`, `X`, `I`, `D`;
- CIGAR query/target consumption must equal the PAF spans;
- `cs:Z` short and long match syntax is supported for `:`, `=`, `*`, `+`, `-`;
- `cs:Z` and `cg:Z` must describe the same match/indel path, not merely consume
  the same total spans;
- splice (`~`) `cs` records are rejected because v0.3 accepts assembly
  alignments, not spliced-read alignments;
- local identity comes from `cs:Z`, or from `cg:Z` only when it contains no
  ambiguous `M` operations;
- without either source, local identity is `.` and is not replaced by
  `matches/block_len`;
- MAPQ 255 follows the PAF specification and is treated as missing;
- an automatically constructed reverse traversal keeps the original
  `evidence_group_id` and is marked `DERIVED_REVERSE`;
- exact duplicate native records are deduplicated before inference.

## Manifests

- **Manifest:** `genome_id fasta te_annotation max_locus_copies`; the final
  positive integer defaults to 1. Genome rows are sorted by ID.
- **Alignment table:** `query_id target_id paf`. Both directions may be listed
  only when they were generated as separate alignment evidence; synthetic
  reverse views do not count as independent reciprocity.
- Relative paths resolve from the directory containing their table.

`max_locus_copies` is a legacy component-capacity constraint, not evidence of
WGD/homeology. v0.4 will replace it with explicit copy contexts derived from
synteny and subgenome metadata.

Newick tree tip names must equal state-matrix genome IDs. Quoted labels and
Newick comments are not supported in this alpha.
