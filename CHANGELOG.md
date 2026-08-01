# Changelog

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
