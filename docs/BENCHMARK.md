# Benchmark and model-preparation contract

TEvoX `0.5.0-alpha.2` introduces a run-independent benchmark contract. It is
an evaluation and dataset-preparation layer, not evidence that the built-in
score is calibrated or accurate on biological data.

The central invariant is:

> Truth belongs to a dataset and biological semantic keys. A `run_id` belongs
> only to predictions and reports.

Truth tables must not contain `candidate_id`, `evidence_id`, `decision_id`,
`locus_id` or any other TEvoX-generated join key. This lets one frozen truth
bundle evaluate multiple methods, ablations, seeds and solver configurations
without duplicating or drifting the truth.

## Multi-run manifests

`datasets.tsv` has these columns:

```text
benchmark_schema_version  dataset_id  genome_truth
candidate_truth  locus_truth  instance_truth
candidate_truth_scope  candidate_sampling_design
candidate_inclusion_probability  locus_truth_scope
```

Paths are resolved relative to `datasets.tsv`. `genome_truth` is mandatory and
freezes, for every `genome_id`, its taxon plus assembly and annotation SHA-256.
Every evaluated run must contain exactly those genomes and the same run-time
hashes. Optional task-truth paths are `.`.

`candidate_truth_scope` is `EXHAUSTIVE`, `PARTIAL` or `NOT_EVALUATED`.
Exhaustive candidate truth requires `candidate_sampling_design=CENSUS` and
`candidate_inclusion_probability=1`; any prediction-only semantic pair is then
a contract error. Partial truth reports such pairs as unassessed and labels
all pairwise score metrics as conditional on the assessed sampling design.
`CASE_CONTROL`, `STRATIFIED` and unknown sampling do not support population
prevalence or calibration claims.

`locus_truth_scope` is `EXHAUSTIVE`, `PARTIAL` or `NOT_EVALUATED`.
Exhaustive truth includes every assessed TE, including biological singleton
loci, so clustering precision cannot be inflated by omitting difficult nodes.

`runs.tsv` has:

```text
benchmark_schema_version  run_id  dataset_id  method_id  condition_id
replicate_id  group_id  prediction_format  prefix
```

`run_id` is globally unique within the manifest. `prefix` is resolved relative
to the run manifest. Alpha.2 accepts `prediction_format=tevox-1.2`; external
methods must first be converted into the same normalized prediction bundle.

## Semantic truth tables

All truth rows contain `benchmark_schema_version=1.0.0` and the matching
`dataset_id`.

### Genome truth

Required fields are:

```text
genome_id  assembly_sha256  annotation_sha256  taxon
```

This dataset-level table prevents identical TE IDs from silently joining truth
across a different assembly or annotation release. `run.json` independently
records the run-time SHA-256 of every FASTA, TE annotation, alignment and
synteny-provider input, including control manifests and MCScanX gene tables.
Hashes are frozen before parsing and rechecked before analysis/output; paths
are stored canonically and remain resolvable from a different working
directory.

### Candidate truth

Required fields are:

```text
truth_record_id
source_genome_id  source_te_id  target_genome_id  target_te_id
label
ancestral_event_id
source_truth_locus_id  target_truth_locus_id
truth_confidence
validation_method  validation_source  validation_batch_id
curation_blinded
taxon_a  taxon_b  clade_holdout_id
```

`label` is `SAME_LOCUS`, `DIFFERENT_LOCUS` or `UNKNOWN`; unknown rows are
tracked but excluded from binary score metrics. `truth_confidence` is `HIGH`,
`MEDIUM` or `LOW`. When locus truth is supplied, the evaluator independently
checks that every candidate label and both locus foreign keys agree with truth
membership.

`SAME_LOCUS` requires a known shared `ancestral_event_id`. For a
`DIFFERENT_LOCUS` pair there is no single shared ancestral event, so `.` is a
valid not-applicable value; its two TE and locus keys still bind both endpoint
events during split auditing.

The semantic key is the unordered pair of `(genome_id, te_id)` nodes. TEvoX
IDs are used only to join prediction files internally. Multiple alignment
views of one semantic pair are reduced by a declared deterministic maximum
fixed-score rule and their multiplicity is reported. The evaluator refuses
candidate benchmarking when `max_candidates != 0`, because a reported top-K
is not a complete candidate universe.

Two score summaries are kept separate:

- generated-only metrics estimate fixed-score behavior conditional on the
  candidate having been generated;
- end-to-end metrics assign zero to ungenerated semantic candidates and
  therefore also measure candidate-generation failure.

Neither summary is a calibration fit.

The separate legacy `tevox-score-audit` can fit diagnostic logistic
intercept/slope summaries. Those are goodness-of-fit diagnostics only; the
tool never writes or blesses a deployable calibrated model.

Candidate summaries are micro-averages over assessed semantic pairs. Pairs
sharing an ancestral event or validation batch are dependent, and alpha.2 does
not yet report event-macro estimates or cluster-bootstrap confidence intervals.
The JSON names AUPRC's definition as tied-threshold step average precision.

### Locus-member truth

Required fields are:

```text
truth_record_id  truth_locus_id  genome_id  te_id
ancestral_event_id  homology_group_id  validation_batch_id
```

The semantic node key is `(dataset_id, genome_id, te_id)`, and a node can occur
in only one truth locus. B-cubed precision/recall/F1, pairwise metrics and ARI
are calculated on the declared node universe. Missing predicted nodes become
zero contributions to end-to-end B-cubed; an intersection-only B-cubed summary
is reported separately. Pairwise metrics and ARI still impute missing members
as unique prediction singletons and therefore must be interpreted together
with member coverage. Under exhaustive truth, prediction-only nodes are a
contract error rather than invented truth singletons; under partial truth they
are counted as unassessed.

### Instance truth

Required fields are:

```text
truth_record_id  truth_locus_id  genome_id
technical_state  biological_state  annotation_state  legacy_state
claimable  contig  start  end
truth_confidence  validation_method  validation_source  validation_batch_id
```

A literal `.` means that an axis or coordinate is not assessed. Otherwise the
state vocabularies equal schema 1.2.0. Coordinates are 0-based half-open and
must be either all present (`contig/start/end`) or all missing.

Alpha.2 state/breakpoint evaluation fails closed when a predicted instance has
`copy_count>1`: the current truth row is keyed by locus × genome and cannot
represent separate WGD-copy states or coordinates. Locus clustering still
supports multi-copy membership. A later copy-level truth schema will add
semantic member/copy-context keys.

Predicted and truth loci are never joined by `TEL` ID. The evaluator builds a
sparse node-overlap graph and applies exact Hungarian maximum-overlap matching
within each connected overlap component, followed by a reciprocal unique-best
and strict-majority purity/completeness gate. Tied or weak mappings are retained
as diagnostics but not used for state scoring. Unmatched truth loci are
retained as `NO_PREDICTION`, rather than disappearing from recall denominators.

## Metrics and outputs

Run:

```bash
tevox-benchmark \
  --datasets benchmark/datasets.tsv \
  --runs benchmark/runs.tsv \
  --output results/benchmark
```

It writes:

- `.metrics.json`: candidate coverage, Brier score, log loss, AUROC, AUPRC,
  equal-mass ECE, B-cubed, ARI, pairwise clustering, per-axis confusion and
  class metrics, empty-site FDR/recall, and breakpoint summaries;
- `.callability_accuracy.tsv`: quality-selective call rate, biological-state
  accuracy and selective risk;
- `.breakpoints.tsv`: per-instance boundary errors and contig agreement.

Truth and normalized prediction inputs are SHA-256 hashed in the report;
dataset genome hashes are checked against the run-time input hashes frozen in
`run.json`. Reports carry `run_id`, method, condition, replicate and group
metadata and are permanently labelled:

```text
EVALUATION_ONLY_NOT_A_CALIBRATED_MODEL
```

The callability curve accepts a row only when it is claimable, its biological
state is not `UNKNOWN`, and quality meets the threshold. Quality and normalized
scores remain uncalibrated ranking quantities.

Assessed empty-site FDR uses only calls with
`legacy_state=EMPTY_SITE_CONFIRMED`,
`biological_state=EMPTY` and `claimable=true`. Predicted empty calls outside
assessed truth are reported separately as `unassessed_empty_site_claims`; they
are never silently discarded and are not folded into an FDR denominator whose
truth status is unknown. The report includes empty-claim assessment coverage,
and biological truth state `UNKNOWN` is always unassessed for this estimand.

## Safe model preparation

`tevox-export-training` requires a deliberately untruncated inference run:

```bash
tevox graph ... --max-candidates 0 --max-graph-candidates 0 --output cohort
mkdir -p model_data
tevox-export-training \
  --prefix cohort --truth truth.candidates.tsv \
  --dataset-id D1 --run-id export-001 --output model_data/cohort
```

The exporter writes **every generated candidate view**, including candidates
without matching truth, with raw pre-decision features and explicit
missingness. `truth_status` distinguishes `ASSESSED_LABEL`,
`ASSESSED_UNKNOWN` and `UNLABELLED`; unlabelled rows use `.` in truth fields
and must not be treated as negatives. It does not export current logits/scores,
eligibility, selected edges, loci, solver results or relation calls as training
features. No candidate downsampling is performed and
`inclusion_probability=1`, while truth ascertainment remains explicitly
unknown. For derived-reverse candidate views, endpoint-specific truth fields
(source/target locus and taxon) are swapped into the exported row's direction.
Its manifest freezes the candidate-generator contract, input hashes and genome
assembly/annotation hashes. Export schema 1.1 rebuilds every raw feature from
the evidence/candidate tables, verifies it against the feature sidecar, and
uses `dataset.json` as a last-written commit marker containing both TSV hashes.

Before any split is used, run the leakage audit:

```bash
tevox-split-audit \
  --truth truth.candidates.tsv --locus-truth truth.members.tsv \
  --splits splits.tsv \
  --strict --output model_data/split_audit.json
```

The split table assigns every `truth_record_id` to a development, calibration
or external-test partition plus outer and inner folds. The validator binds the
same TE, reciprocal pair, ancestral event, truth locus, homology group,
validation batch and clade in one partition/fold, and by default prohibits
taxa shared between the external test and model-development data. Strict mode
requires locus-member truth, verifies both candidate endpoint locus foreign
keys, and fails when a locus homology-group binding is unknown. Omitting
`--locus-truth` is permitted only for a non-strict exploratory report, which is
explicitly marked unverifiable for endpoint event/HMG/batch leakage.

## Non-claims

The repository does not ship a trained calibrated model in alpha.2. Until
independent truth, pre-registered leakage-safe splits and external-clade
evaluation exist, TEvoX must not:

- call `membership_score` a probability or posterior;
- refit built-in coefficients on regression fixtures;
- tune a threshold and report performance on the same data;
- transfer candidate score behavior to state, locus or relation accuracy;
- claim WGD-homeolog/allelic/tandem/segmental/transposed relation accuracy:
  alpha.2 has no semantic relation truth or per-class relation benchmark;
- claim copy-level state or breakpoint accuracy for `copy_count>1`;
- claim an empty-site FDR, out-of-domain guarantee or method superiority on
  the basis of the bundled software fixtures.

Fixtures under `tests/data/` establish implementation correctness only. They
are not independent biological evidence.
