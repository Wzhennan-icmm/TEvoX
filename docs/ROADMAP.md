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

## Probabilistic inference — planned 0.5

- missing-aware feature table and calibrated hierarchical state model;
- block-level maximum-weight matching;
- copy-context-constrained correlation clustering/ILP;
- explicit `ORTHOLOG`, `WGD_HOMEOLOG`, `ALLELIC`, `TANDEM_PARALOG`,
  `SEGMENTAL_PARALOG`, `TRANSPOSED_PARALOG` and `UNKNOWN` posteriors;
- prediction sets, entropy, out-of-domain flags and solver status/gap;
- species-pair holdout and calibration evaluation.

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
