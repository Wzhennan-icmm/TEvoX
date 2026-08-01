# Validation strategy

## Automated v0.4 coverage

`make check` compiles with strong warnings and covers:

- annotated, unannotated, empty-site, family-conflict, assembly-gap and reverse
  PAF cases;
- exact `cs`/`=X` local identity versus aggregate PAF identity;
- missing MAPQ 255, contig-edge flanks and derived reverse provenance;
- rejection of zero-overlap candidates, family bridges and inconsistent
  query/target direction, plus non-finite numeric options;
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
- top-K=1 versus unlimited winner/decision/locus equivalence under score ties;
- standalone alignment and TE interval-index boundary/counter tests;
- schema `1.1.0` primary/foreign keys, provider-specific missing semantics,
  context chains, row counts and zero FASTA reopen after indexing;
- golden compatibility output and exploratory phylogenetic events.

`make asan` rebuilds the program plus the standalone index and delta-importer
tests under AddressSanitizer/UndefinedBehaviorSanitizer, then repeats the
complete suite. Leak detection is disabled in ptrace-restricted containers.

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
- calibrated probabilities after the heuristic model is replaced;
- threshold sensitivity and DNA-only/synteny-only/combined ablations;
- stratification by TE class, divergence, repeat density, nesting, assembly
  quality, ploidy and structural context.
