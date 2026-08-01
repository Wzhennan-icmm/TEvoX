# Changelog

## 0.4.0-alpha.1 — 2026-08-01

- added strict MUMmer4 NUCMER `.delta` import with file-header direction
  verification, while keeping aggregate alignment identity, local identity and
  mapping confidence as separate evidence fields;
- added strict MCScanX collinearity import and normalized gene tables as an
  independent synteny prior that cannot by itself produce base-level state
  calls or confirmed empty sites;
- introduced stable homology groups and copy contexts across genome,
  subgenome, haplotype, syntenic copy and WGD-node metadata;
- constrained graph merges by copy context and made metadata conflicts
  explicit rather than silently choosing one annotation;
- added sparse FASTA gap indexes, alignment/TE interval indexes, deterministic
  bounded top-K candidates and deduplicated edge-support aggregation;
- introduced schema 1.1.0 sidecars for synteny blocks, anchors, copy contexts
  and TE-to-context assignments, with provider provenance and performance
  counters;
- expanded regression tests for delta parsing, MCScanX order invariance, WGD
  copy slots, context ambiguity, top-K equivalence and index boundaries.

## 0.3.0-alpha.1 — 2026-08-01

- corrected minimap2/PAF direction documentation and added direction checks;
- added exact local identity parsing from `cs:Z` and `=X`, with explicit NA for
  ordinary `M` CIGARs;
- treated MAPQ 255 as missing and grouped native/derived reverse evidence so a
  synthetic reverse cannot earn reciprocal support;
- required 50% reciprocal TE overlap, fixed contig-edge flanks, prevented
  adjacent-TE merges and unknown-family bridging;
- made graph inference invariant to manifest, annotation and PAF record order;
- introduced schema 1.0.0 with evidence, candidate, decision, instance and run
  provenance tables plus three-axis state and claimability gates;
- expanded regression and schema-integrity tests for all v0.3 correctness bugs.

## 0.2.0-alpha.1 — 2026-08-01

- replaced “outside synteny equals unique TE” with locus-level CIGAR evidence;
- added eight uncertainty-aware states and strict input validation;
- added multi-genome graph, component copy quotas and quality diagnostics;
- added Newick/Fitch event candidates, golden tests, sanitizer and CI;
- removed the tracked binary, legacy CLI and obsolete fixtures.
