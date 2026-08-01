# Input contracts

- **FASTA:** plain, unique contig IDs, fixed-width lines except the last line.
  FASTA is indexed in place rather than loaded in memory.
- **GFF3:** paths containing `.gff` are parsed as GFF3. Supply a TE-only file;
  `ID`/`Name` identifies a TE and `family`/`Family`/`classification` supplies
  family. Coordinates are converted from 1-based closed.
- **BED:** `chrom start end id score strand family type`; the first three
  zero-based half-open columns are required.
- **PAF:** 12 core fields plus mandatory `cg:Z`. Supported operations are
  `M`, `=`, `X`, `I`, `D`; contig lengths and CIGAR consumption must match.
- **Manifest:** `genome_id fasta te_annotation max_locus_copies`; the last
  positive integer defaults to 1.
- **Alignment table:** `query_id target_id paf`. Do not list both directions.

Recommended alignment: `minimap2 -cx asm5 --cs=long query.fa target.fa`.
Newick tree tip names must equal state-matrix genome IDs. Quoted labels and
Newick comments are not supported in this alpha.
