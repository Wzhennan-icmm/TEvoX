# Roadmap

## P0/P1 foundation — 0.2.0-alpha.1

Introduced CIGAR projection, paired-flank evidence, strict input validation,
multi-genome manifests, a quota-constrained locus graph, compatibility states,
golden tests, CI and exploratory Fitch events.

## Correctness and evidence contract — 0.3.0-alpha.1

- corrected minimap2/PAF direction and exact local-identity semantics;
- represented MAPQ 255, ordinary `M`, contig-edge flanks and missing evidence
  explicitly;
- prevented synthetic reciprocity, adjacent-TE merging and family bridges;
- added three-axis state, claimability and schema `1.0.0` foreign keys.

## Synteny, delta and copy contexts — 0.4.0-alpha.1

- strict MUMmer4 NUCMER `.delta` import with provider-specific identity/MAPQ;
- strict MCScanX block/anchor import from a normalized gene table;
- stable block, anchor, HMG, context and TE-assignment identities;
- explicit genome/subgenome/haplotype/syntenic-copy/WGD context constraints;
- MCScanX-only safety: no base-level or confirmed empty-site claim;
- alignment, TE and gap-run indexes; deterministic bounded top-K candidates;
- sorted/reduced edge support and linear reciprocal aggregation;
- schema `1.1.0`, context sidecars, provider provenance and performance counters.

Remaining v0.4 hardening before a beta tag:

- validated stitching of adjacent fragmented MCScanX blocks;
- a side-qualified gene-ID adapter that does not require global prefixes;
- indexed TE-to-context construction and incremental DSU metadata at very
  large graph scale;
- measured runtime/RAM scaling on plant chromosomes;
- optional GENESPACE and AnchorWave adapters after their semantics are frozen.

## Constrained score inference — 0.5.0-alpha.1

- missing-aware candidate features and three-axis observation score tables;
- explicit `BUILTIN_UNCALIBRATED_V1` model/calibration labelling;
- copy-context-pair maximum-weight Hungarian matching with a labelled
  deterministic size fallback;
- constrained component partitioning: exhaustive optimum for bounded small
  components, deterministic heuristic plus upper-bound gap otherwise;
- normalized scores for `ORTHOLOG`, `WGD_HOMEOLOG`, `ALLELIC`,
  `TANDEM_PARALOG`, `SEGMENTAL_PARALOG`, `TRANSPOSED_PARALOG` and `UNKNOWN`;
- prediction sets, entropy, out-of-domain flags and stable edge/matching/
  relation/solver keys in schema `1.2.0`;
- a truth-controlled fixture proving a case where exact global optimization
  improves on the greedy fallback;
- an independent-truth score audit for Brier/log-loss/ECE/AUROC and group-wise
  evaluation that cannot mark the built-in model calibrated.
- separate deterministic limits for reported and graph-inference candidates,
  preventing report truncation from changing reconstructed loci.

## Evaluation and model preparation — 0.5.0-alpha.2

- dataset-semantic truth tables plus multi-run method/condition/replicate
  manifests, with no generated-ID truth joins and frozen assembly/annotation
  hashes;
- candidate coverage and score diagnostics, B-cubed/ARI/pairwise locus
  metrics, exact sparse truth/prediction locus matching, three-axis state
  metrics, empty-site FDR/recall, breakpoint error and callability–accuracy;
- untruncated raw-feature export with candidate-generator/input/assembly
  hashes, no downstream decision features and no negative downsampling;
- partition, nested-fold and external-taxon leakage auditing across shared TEs,
  ancestral events, loci, homology groups, validation batches and clades,
  with strict candidate/locus-truth endpoint validation;
- AUPRC, equal-mass ECE, calibration intercept/slope, coverage and optional
  group-bootstrap intervals in the legacy run-local score audit;
- fixed relation locus assignment for non-selected direct edges whose members
  meet through another path, and separated nested versus tandem relations;
- blocked distinct same-contig TE bridge merges, including nested annotations,
  while retaining an explicit same-HMG/WGD-context exception only for known
  distinct subgenomes, never allelic haplotypes;
- froze and rechecked every raw/control input, rejected path/inode aliases and
  made `run.json` an atomic last-written completion marker;
- reserved `.tvm` contract documented without shipping or fitting a calibrated
  model.

Remaining before v0.5 beta:

- train and freeze a calibrated model on independent truth with species-pair
  or clade holdout; until then every score remains explicitly uncalibrated;
- correct graph top-K ordering so context-aware membership is computed before
  pruning, and make matching plus component partitioning a joint or iterative
  optimization rather than an irreversible prefilter;
- replace exhaustive small-component search with a production ILP/correlation-
  clustering backend if benchmarks justify the dependency;
- add independent duplication evidence before making confident segmental or
  transposed-paralog calls;
- add semantic relation truth and WGD/copy-context per-class metrics; alpha.2
  evaluates locus membership but makes no relation-accuracy claim;
- replace locus×genome instance truth with copy-level semantic member/context
  truth before evaluating multi-copy state and breakpoint accuracy;
- add event/locus-macro estimates and validation-batch/event cluster bootstrap
  intervals so large multi-genome loci do not dominate pairwise metrics;
- extend provider provenance beyond format/path hashes to external tool
  version, command/preset and preprocessing fingerprints before training a
  portable model;
- replace repeated locus×genome full-table output scans with locus/member and
  decision indexes before claiming million-TE scalability;
- benchmark Hungarian fallback thresholds and component constraint checks on
  chromosome-scale plant graphs;
- add DNA-only/synteny-only/combined ablation reports and score calibration
  audit plots.

## Publication validation — planned 0.9

- truth-controlled divergence, gap, rearrangement, nesting, tandem and WGD
  simulations;
- a public diploid pangenome and independently validated allopolyploid case;
- comparisons with panREPET, GraffiTE, Minigraph and cross-species projection
  baselines;
- state and locus accuracy, empty-site FDR, calibration, ablations and scaling;
- frozen containers/workflow/data and one-command manuscript reproduction.

Branch-length-aware stochastic Dollo inference follows only after locus calls
are validated.
