# Roadmap

## P0/P1 foundation — 0.2.0-alpha.1

Introduced CIGAR projection, paired-flank evidence, strict input validation,
multi-genome manifests, a quota-constrained locus graph, compatibility states,
golden tests, CI and exploratory Fitch events.

## Correctness and evidence contract — 0.3.0-alpha.1

- corrected minimap2/PAF direction documentation and validation;
- exact local identity from `cs`/`=X`, with explicit missing values;
- MAPQ 255 missing semantics and provenance-preserving reverse views;
- ≥50% reciprocal-overlap candidate eligibility;
- contig-edge, adjacent-TE, family-bridge and input-order fixes;
- three-axis state plus a separate claimability gate;
- schema `1.0.0`: evidence, candidates, decisions, edges, loci, instances and
  run provenance with foreign-key validation.

## Synteny and copy contexts — planned 0.4

- MCScanX and GENESPACE adapters as gene-synteny priors;
- MUMmer `.delta` and AnchorWave import adapters;
- explicit `(genome, subgenome, haplotype, syntenic copy, WGD node)` context;
- replacement of `max_locus_copies` by context-specific copy slots;
- distinction among ortholog, WGD homeolog, allelic, tandem, segmental and
  transposed paralog relations;
- interval indexes, streaming PAF and deterministic top-K candidates.

MCScanX-only evidence will never be allowed to produce a confirmed empty-site
call.

## Probabilistic inference — planned 0.5

- missing-aware feature table and calibrated hierarchical state model;
- block-level maximum-weight matching;
- copy-context-constrained correlation clustering/ILP;
- posterior probabilities, prediction sets, entropy, out-of-domain flags and
  solver status/optimality gap.

## Publication validation — planned 0.9

Simulation, public truth data, a diploid pangenome and an allopolyploid case;
independent read/PCR or curated validation; competitive baselines, ablations,
calibration, scaling, frozen containers/workflow and one-command figure
reproduction. Branch-length-aware stochastic Dollo inference follows only after
locus calls are validated.
