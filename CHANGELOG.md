# Changelog

## Unreleased integration — 2026-10-06

- reject invalid GFF coordinates before converting a one-based start, so
  `INT64_MIN` cannot overflow; preserve diagnostics and completion-marker gates;
- normalize native negative-strand PAF cg/cs operations to query-forward
  traversal, fixing asymmetric indel projection;
- report a missing local identity method when a local window contains no
  aligned bases, even if the full alignment has an exact CIGAR; inference
  still treats these unsupported windows as uncallable;
- index per-locus output rows and long CIGAR paths; differential tests compare
  all 17 TSVs with linear traversal on both strands;
- stream TSV output through checked gzip child processes and teach Python
  readers to read compressed tables without materializing them;
- support isolated builds, read-only sources, Python 3.8+, staged helper
  installation and native Intel/ARM/Universal macOS packages;
- retain the Hungarian exact threshold: **both sides** must contain at most
  `--exact-match-nodes` nodes, not merely one side.
- add version-bound published-genome preparation, resume-integrity checks,
  lossless result archiving, external GIAB/PCR state evaluation and controlled
  copy-context/scale studies; preserve uncalibrated model status.

## 0.5.0-alpha.2 — 2026-08-02

- added a run-independent semantic truth contract and multi-run evaluator for
  candidate generation/scores, B-cubed/ARI/pairwise locus reconstruction,
  three-axis states, empty sites, breakpoints and callability–accuracy;
- froze dataset genome/taxon truth and pre-parse FASTA, annotation, alignment,
  manifest and MCScanX gene/collinearity SHA-256 values; inputs are rechecked
  before analysis/output, recorded with canonical absolute paths, and
  benchmark runs fail on assembly drift;
- made candidate truth scope/sampling explicit, counted partial-truth
  prediction-only pairs, accepted zero-prediction methods as all-missing, and
  prevented missing nodes/classes and tied score bins from inflating metrics;
- gated state scoring on reciprocal unique-majority locus mappings, retained
  missing calls as `NO_PREDICTION`, excluded biological `UNKNOWN` from empty
  FDR, and fail closed on unsupported multi-copy state coordinates;
- added safe untruncated raw-feature export with candidate-generator,
  assembly, annotation and input checksums, while excluding downstream model,
  selection, solver and locus fields; export schema 1.1 rebuilds raw features,
  preserves unlabelled candidates and commits hashed TSVs via a JSON marker;
- added partition/fold leakage auditing that binds shared TEs, reciprocal
  pairs, ancestral events, loci, homology groups, validation batches and
  clades, with strict candidate/locus-truth cross-checks and external taxon
  isolation by default;
- expanded the run-local uncalibrated score audit with AUPRC, equal-mass ECE,
  diagnostic calibration intercept/slope, truth coverage and deterministic
  group-bootstrap confidence intervals;
- corrected direct relation rows to use final component membership even when
  the direct edge was not selected, and stopped classifying nested/overlapping
  annotations as tandem copies;
- blocked distinct same-genome/same-contig TEs—including adjacent, overlapping
  and nested annotations—from entering one locus through a third-node bridge,
  except under explicit distinct strong contexts in one known HMG/WGD node and
  known distinct subgenomes; allelic haplotypes cannot invoke the exception,
  and gated WGD pairs cannot be relabelled tandem by proximity;
- rejected input/output and output/output path, symlink and hard-link aliases,
  checked every output close, and made `run.json` an atomic last-written marker
  so failed reruns cannot retain a stale completion record;
- documented the reserved `.tvm` frozen-model contract without fitting,
  loading or shipping a calibrated model; schema 1.2.0 and
  `BUILTIN_UNCALIBRATED_V1` remain unchanged.

## 0.5.0-alpha.1 — 2026-08-01

- added missing-aware candidate features and normalized three-axis observation
  scores under the explicit `BUILTIN_UNCALIBRATED_V1` contract;
- added exact Hungarian one-to-one matching within strong MCScanX copy-context
  pairs, including alternative-candidate assignment, hard-gate prefiltering
  and a deterministic, visibly labelled large-block fallback;
- replaced the greedy locus forest with component-wise constrained partition
  optimization: exhaustive optimal search for bounded components and a
  deterministic heuristic with objective, upper bound and relative gap for
  larger components;
- added normalized relationship scores for ortholog, WGD homeolog, allelic,
  tandem, segmental, transposed and unknown classes, including prediction sets,
  entropy and out-of-domain flags;
- introduced schema 1.2.0 sidecars for observation scores, candidate features,
  relations and solver diagnostics plus stable edge/matching/relation/solver
  identities;
- separated report candidate truncation from graph candidate pruning so output
  size controls cannot silently alter inferred edges or loci;
- added a truth-controlled graph where the exact solver recovers the global
  optimum missed by greedy selection, plus matching fallback, order-invariance
  and schema-integrity regression tests;
- retained v0.4 state/claimability gates and prohibited probability/posterior
  wording until independent holdout calibration exists.

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
