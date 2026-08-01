# Validation strategy

## v0.3 regression coverage

`make check` currently covers:

- annotated presence, `cs`-supported empty sites and unannotated sequence;
- known-family conflict, target assembly gap and reverse-strand projection;
- local `cs` identity versus whole-PAF identity and missing identity for `M`;
- PAF MAPQ 255 and low MAPQ preserved through the derived reverse view;
- contig-edge flanks represented as missing/uncallable;
- ≥50% reciprocal-overlap eligibility and rejection of adjacent zero-overlap
  TEs even when the legacy quota is 2;
- source-level multi-projection ambiguity with every near-best observation
  linked to the decision;
- independent versus same-group reciprocal evidence;
- prevention of unknown-family bridging between incompatible known families;
- invariance to annotation, manifest and PAF record order;
- declared PAF query/target and malformed-record rejection;
- evidence-schema primary/foreign keys, missing semantics and `run.json` row
  counts;
- golden compatibility output and exploratory phylogenetic event inference.

`make asan` repeats the suite with AddressSanitizer and
UndefinedBehaviorSanitizer. Leak detection is disabled for ptrace-restricted
containers; an unrestricted CI environment should add LeakSanitizer.

## Evidence still required before publication

Code-level regression tests do not establish biological accuracy. Before a
Bioinformatics submission, later releases must add:

- simulated truth sets spanning divergence, nested/fragmented TE, annotation
  errors, contig edges, assembly gaps, structural rearrangements, tandem
  duplication and WGD;
- curated real TE-PAV truth plus read/PCR validation;
- comparisons with panREPET, GraffiTE, Minigraph and relevant cross-species
  projection baselines;
- per-state precision/recall/F1, empty-site false-discovery rate, B-cubed/ARI
  locus clustering, breakpoint error and callability–accuracy curves;
- probability calibration after the heuristic model is replaced;
- runtime/RAM scaling and threshold sensitivity;
- ablations for DNA alignment, synteny context, family information,
  independent reciprocity and assembly-quality evidence.

Results must be stratified by TE class, divergence, repeat density, nesting,
assembly quality, ploidy and structural context.
