# Validation strategy

## Automated v0.5 coverage

`make check` compiles with strong warnings and covers:

- annotated, unannotated, empty-site, family-conflict, assembly-gap and reverse
  PAF cases;
- exact `cs`/`=X` local identity versus aggregate PAF identity;
- missing MAPQ 255, contig-edge flanks and derived reverse provenance;
- rejection of zero-overlap candidates, family bridges and inconsistent
  query/target direction, plus non-finite numeric options;
- rejection of real reciprocal=0.5 adjacent and nested/overlapping TE bridges
  under a permissive copy quota, with a separately validated explicit
  same-HMG/WGD exception and an allelic same-subgenome/haplotype negative
  control;
- near-best projection ambiguity and stable output under PAF/GFF/manifest order;
- MUMmer NUCMER plus/reverse operation normalization, exact half-open
  coordinates, 7/9 indel-aware aggregate identity, shared reverse provenance,
  strict header direction, absent local-identity/MAPQ semantics and
  malformed/PROMER rejection;
- MCScanX plus/minus blocks, checked alignment/rank prefixes, unknown genes and
  malformed input rejection;
- invariance to gene/block/anchor order and swapped MCScanX sides;
- one A context versus two B WGD contexts in one HMG;
- two context copies retained when the fallback quota is one;
- same-context, cross-HMG bridge and metadata-ambiguous constraints;
- separation of overlapping evidence assigned to different WGD nodes;
- MCScanX-only inability to emit a confirmed empty-site call;
- report top-K=1 versus unlimited edge/locus equivalence, plus separately
  recorded graph-candidate pruning under score ties;
- missing-aware feature masks, normalized axis/relation scores, entropy,
  prediction sets, OOD and mandatory `UNCALIBRATED` model labels;
- exact Hungarian selection of a second-ranked local alternative, exclusion
  of hard-gated edges from copy slots, conflict rejection, and deterministic
  oversized-block fallback;
- exact component objective/bound/gap reporting and a controlled graph where
  the exact partition beats the greedy fallback;
- run-local score auditing for coverage, Brier/log-loss, AUROC/AUPRC,
  equal-width/equal-mass ECE, diagnostic intercept/slope and deterministic
  group-bootstrap intervals without upgrading calibration status;
- one semantic truth bundle evaluating exact and heuristic runs, with the
  controlled greedy failure detected by B-cubed and ARI;
- frozen assembly/annotation hashes, exhaustive/partial candidate scope,
  prediction-only accounting and zero-prediction end-to-end B-cubed behavior;
- pre-parse and rechecked hashes for all raw/control inputs, input/output and
  output/output alias rejection (including symlink/hard-link cases), output
  close failures, cross-working-directory export and atomic last-written
  `run.json` recovery;
- sparse exact truth/prediction locus matching with reciprocal-majority gates,
  three-axis state accuracy, explicit `NO_PREDICTION`, UNKNOWN-safe assessed
  empty FDR, breakpoint error and monotone quality-selective callability;
- split leakage acceptance/rejection across TE/event/locus/HMG/batch/clade/fold
  bindings, strict failure without locus truth and explicit no-model-fitted
  reports;
- raw training-feature export from an untruncated run, with every generated
  view retained, reverse endpoint truth oriented correctly, feature/evidence
  reconstruction, run-time input-hash verification, commit-marker output
  hashes, forbidden downstream fields absent and fail-closed tamper/truncation
  checks;
- explicit WGD-homeolog and tandem relation rows plus stable edge, matching,
  relation and solver foreign keys;
- final-locus relation labelling for a weak direct edge connected through a
  third node, and no tandem label for nested/overlapping annotations;
- standalone alignment and TE interval-index boundary/counter tests;
- schema `1.2.0` primary/foreign keys, provider-specific missing semantics,
  score/logit normalization, matching/solver invariants, context chains, row
  counts and zero FASTA reopen after indexing;
- golden compatibility output and exploratory phylogenetic events.

`make asan` rebuilds the program plus the standalone index and delta-importer
tests under AddressSanitizer/UndefinedBehaviorSanitizer, then repeats the
complete suite. Both sanitizers halt on their first error. Leak detection is
disabled in ptrace-restricted containers.

The release gate also runs:

```bash
make clean && make check
make asan
make clean && make CC=clang check
make clean && make CFLAGS='-O0 -g -fanalyzer' check
git diff --check
```

## Performance evidence

`run.json` reports alignment/TE interval queries, records actually examined,
gap queries, post-index FASTA reopens and edge-support records. These counters
verify that indexed code paths are active; they are not substitutes for wall
time and peak-RAM benchmarks.

Before a beta or manuscript claim, add synthetic scaling fixtures and real
chromosome tests covering at least 50k genes, 100k anchors, millions of
alignments and dense TE neighborhoods. Report input size, retained K, wall
time, peak RSS and candidate pruning ratio.

## Evidence still required before publication

Code regression does not establish biological accuracy. Later releases must
add:

- simulations spanning divergence, nested/fragmented TE, annotation errors,
  assembly gaps, rearrangements, tandem duplication and WGD;
- curated TE-PAV truth plus read/PCR or independent long-read validation;
- comparisons with panREPET, GraffiTE, Minigraph and relevant cross-species
  projection methods;
- per-state precision/recall/F1, empty-site FDR, B-cubed/ARI locus clustering,
  breakpoint error and callability–accuracy curves;
- calibrated probabilities trained and evaluated under species-pair/clade
  holdout; the built-in v0.5 score is explicitly not calibrated;
- threshold sensitivity and DNA-only/synteny-only/combined ablations;
- stratification by TE class, divergence, repeat density, nesting, assembly
  quality, ploidy and structural context.
