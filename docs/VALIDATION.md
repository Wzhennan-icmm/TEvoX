# Validation strategy

`make check` covers shared loci, query-insertion empty sites, annotation dropout, family conflicts,
assembly gaps, forward/reverse alignments, malformed PAF rejection,
three-genome reconstruction, WGD quota 2 versus quota 1, golden output and
phylogenetic event inference. `make asan` repeats the suite with Address- and
UndefinedBehaviorSanitizer; leak detection is disabled for ptrace-restricted
containers.

Before a Bioinformatics submission, P2 must add simulated and curated truth
sets; comparisons with relevant graph/pangenome and TE-presence baselines;
per-state precision/recall and locus-clustering metrics; quality-score
calibration; runtime/RAM scaling; threshold sensitivity; and ablations for
flank, family, reciprocal, quota and assembly-gap terms. Results must be
stratified by TE class, divergence, repeat density, nesting, assembly quality
and WGD status.
