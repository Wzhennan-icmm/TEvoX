# Frozen candidate-model format (reserved)

TEvoX alpha.2 reserves the extension `.tvm` for a future, dependency-free,
strictly validated candidate-membership model. No calibrated `.tvm` is shipped
or loaded by the inference engine in this release.

The first supported format will be a canonical UTF-8 TSV stream with this
logical structure:

```text
TEVOX_MODEL  1
META  model_id  ...
META  task  CANDIDATE_SAME_LOCUS
META  feature_schema_version  ...
META  candidate_generator_contract_sha256  ...
META  training_dataset_sha256  ...
META  split_manifest_sha256  ...
NUM   feature_name  center  scale  coefficient  missing_coefficient
CAT   feature_name  level  coefficient
CAL   PLATT  intercept  slope
DOMAIN  ...
END   sha256_of_all_preceding_canonical_bytes
```

A future loader must reject unknown, duplicate or missing features; non-finite
parameters; incompatible evidence/feature schemas; candidate-generator
contract mismatches; malformed canonical ordering; and checksum failures.
Out-of-domain rows or contract mismatches must produce a missing calibrated
field plus an explicit reason, never silently fall back while retaining a
`CALIBRATED` label.

The eventual model will add a separate
`candidate_membership_probability` field. Existing `membership_score` remains
the versioned `BUILTIN_UNCALIBRATED_V1` score so its meaning does not change.
Calibration of candidate membership will not imply calibrated observation
states, selected-edge probabilities, reconstructed-locus confidence or
relationship probabilities; each task requires its own truth and validation.
