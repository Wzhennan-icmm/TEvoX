# Evidence schema 1.2.0

Schema `1.2.0` is the public contract for TEvoX `0.5.x`. All coordinates are
0-based half-open. A literal `.` means unobserved or not applicable; it never
means numeric zero.

## Stable keys

| Prefix | Entity | Biological key |
|---|---|---|
| `EVG` | provider evidence group | normalized native PAF/delta record or MCScanX block |
| `EVD` | directed observation | alignment group × source TE × target genome |
| `CAN` | target candidate | observation × target TE |
| `DEC` | decision | source TE × target genome |
| `SBL` | synteny block | normalized source/WGD/orientation/anchor pairs |
| `SYN` | synteny anchor | block × normalized gene pair |
| `CTX` | copy context | genome/contig/span/anchor set/subgenome/haplotype/WGD |
| `HMG` | context homology group | connected context IDs |
| `TEC` | TE-context assignment | TE × context |
| `EDG` | reduced candidate edge | canonical TE pair |
| `MAT` | block matching group | ordered stable copy-context pair |
| `SOL` | solver component | candidate-connected stable TE set |
| `REL` | relationship score row | canonical TE pair |
| `TEL` | reconstructed locus | selected global component |
| `INS` | locus instance | locus × genome |

IDs do not depend on physical row order or file path. Provider record and line
remain separate provenance fields. `copy001` is a display label, not a stable
key; use `context_id` for joins.

## Alignment observations

### `evidence.tsv`

One row per directed alignment observation. The provider section distinguishes:

- `provider`, path, record and physical line;
- `NATIVE/INDEPENDENT` from `DERIVED_REVERSE/DERIVED_SAME_GROUP`;
- query/target genome and alignment coordinates.

Identity and mapping fields are deliberately separate:

| Provider | `alignment_identity_method` | `identity_method` | `mapq_status` |
|---|---|---|---|
| PAF | `PAF_CORE` | `CS`, `EQX` or `MISSING` | `OBSERVED` or `MISSING_255` |
| MUMmer delta | `DELTA_ERROR_COUNT` | `MISSING` | `NOT_PROVIDED` |

`alignment_identity` describes the complete provider record.
`local_identity` describes only the projected TE/flank window and is never
imputed from aggregate identity. Provider-specific error fields are `.` when
not supplied.

`mapping_confidence` is one of `NOT_ESTABLISHED`, `MAPQ_PASS`,
`UNIQUE_ALIGNMENT`, `UNIQUE_COPY_CONTEXT`, `AMBIGUOUS`, or `FAILED`. It is an
explicit claim gate, not a probability.

`nearby_candidate_count` counts every exact-window candidate examined;
`retained_candidate_count` is the deterministic top-K written to candidate
tables; `graph_candidate_count` is the separately bounded top-K admitted to
matching and global inference. Eligibility, the best and second-best score,
and ambiguity are computed over all candidates before either truncation.

### `candidates.tsv` and `candidate_contexts.tsv`

`candidates.tsv` contains retained target annotations ordered by the single
complete comparator:

1. eligible before ineligible;
2. decreasing score;
3. decreasing reciprocal overlap;
4. increasing breakpoint distance;
5. stable candidate ID.

Ranks are contiguous within the reported prefix and exactly one eligible winner
may have `selected=true`. `graph_retained` says whether a reported candidate
also entered inference; changing `--max-candidates` alone cannot change edges
or loci. Reciprocal overlap is

\[
\min\left(\frac{|P\cap T|}{|P|},\frac{|P\cap T|}{|T|}\right).
\]

`candidate_contexts.tsv` is a one-to-one extension keyed by `candidate_id`. It
contains source/target context foreign keys, a shared HMG when present, and
`SUPPORTED`, `CONFLICT`, `AMBIGUOUS`, or `UNKNOWN`.

### `observation_scores.tsv` and `candidate_features.tsv`

These are one-to-one sidecars keyed by `evidence_id` and `candidate_id`.
Every row records `model_id=BUILTIN_UNCALIBRATED_V1` and
`calibration_status=UNCALIBRATED`.

`observation_scores.tsv` contains normalized scores, entropy and a cumulative
score prediction set for each state axis. `candidate_features.tsv` contains an
observed bit mask, named missing fields, raw feature values, membership logit,
normalized score, entropy, prediction set and out-of-domain flag. A missing
feature remains `.` and contributes through an explicit missing term; it is
never imputed as numeric zero.

These normalized values are not calibrated probabilities. The legacy state
and `claimable` fields remain the state-call contract. Exact coefficients and
non-claims are specified in [the v0.5 inference contract](V05_INFERENCE.md).

### `decisions.tsv`

One row per source TE × target genome. `winner_evidence_id` is the deterministic
observation winner. `linked_evidence_ids` lists every near-best observation
used for consistency testing, and its cardinality equals `near_best_count`.

### `edges.tsv`

One row per undirected annotated TE pair with stable `edge_id`. Evidence IDs
are listed once and evidence groups are deduplicated. `support_count` equals
the number of unique alignment groups. `independent_reciprocal=true` requires
native opposite directions from different groups; an automatic reverse never
qualifies.

Membership fields retain the model/calibration label. An otherwise eligible
strong-context edge also has `matching_group_id`, `matching_method` and
`matching_selected`; a pre-matching hard-gate failure remains
`NOT_APPLICABLE` and records its gate-specific selection reason.
Selected graph edges reference a `solver_component_id`; rejected rows preserve
the exact gate or global-separation reason.

### `relations.tsv` and `solver.tsv`

`relations.tsv` contains direct edge relations plus implicit within-locus and
nearby tandem relations. It exposes a seven-class normalized score vector,
prediction set, entropy, OOD flag, and optional edge/locus foreign keys.

`solver.tsv` contains one row for every candidate-connected component,
including trivial singletons. `EXACT_ENUMERATION/OPTIMAL` rows have
`objective == upper_bound` and zero gap. `DETERMINISTIC_GREEDY/HEURISTIC` rows
report the sum of all positive edge weights as a valid, potentially loose upper
bound. Matching and global fallbacks are never labelled optimal.

`OPTIMAL` is conditional on the retained candidate graph and preceding block
matching. Candidate pruning occurs before matching, and component constraints
do not cause a rejected block-matching alternative to be reconsidered when a
previously matched edge is blocked. A zero component gap therefore does not
certify a joint optimum over all original candidates, matchings and partitions.

## Synteny and copy context

### `synteny.blocks.tsv`

One validated MCScanX block with provider record/line, normalized side
coordinates, orientation, score/E-value, WGD node, HMG and two context foreign
keys. MCScanX groups occupy the same `EVG` namespace but cannot appear as an
alignment observation.

### `synteny.anchors.tsv`

One normalized gene pair per block. `anchor_rank` is canonical coordinate order;
`provider_anchor_rank` preserves the checked MCScanX `0..N-1` rank. The third
provider field is stored as `reported_evalue`.

### `contexts.tsv`

One copy context with stable context/HMG IDs, genome interval, subgenome,
haplotype, syntenic display copy, WGD node, anchor count and status. A metadata
conflict is `METADATA_AMBIGUOUS` and cannot act as a strong slot.

### `te_contexts.tsv`

One overlapping TE × context assignment:

- `BRACKETED`: consecutive anchors enclose the complete TE;
- `BLOCK_INTERIOR`: overlap without a valid two-anchor bracket;
- `AMBIGUOUS`: more than one strong context brackets the TE.

Only `BRACKETED` assignments in `PASS` contexts are strong. Block-interior and
missing context are neutral; ambiguous context cannot be converted to support.

## Loci, instances and compatibility

- `loci.tsv` describes selected graph components.
- `instances.tsv` is the primary locus × genome biological output with the
  three state axes, claim type, claimability, decisions and evidence keys.
- `states.tsv` retains the eight-state compatibility matrix.
- `summary.tsv` counts compatibility states by genome.
- `run.json` records versions, parameters, provider groups, row counts,
  inference model/calibration status, performance counters and the complete
  output inventory. Each genome entry freezes `fasta_sha256` and
  `te_annotation_sha256`; every alignment/synteny evidence entry freezes its
  provider `path_sha256`. Its `input_files` inventory also freezes genome,
  alignment and synteny manifests plus MCScanX gene/collinearity inputs, while
  `synteny_inputs` binds both provider tables explicitly. Inputs are hashed
  before parsing and rechecked before analysis and output; drift is fatal.
  Recorded input paths are canonical absolute paths, so later validation does
  not depend on the process working directory.
  Output/input and output/output aliases are rejected before writing, and
  `run.json` is an atomic last-written completion marker.

`run.json` distinguishes the exact-window `candidate_observations`, union
`internal_candidates`, `graph_candidates`, and reported `candidates` counts.
`max_candidates` controls only report size; `max_graph_candidates` controls the
inference graph and therefore can change loci.

Strong copy-context constraints replace the legacy quota for assigned TEs:
one ancestral locus may contain at most one member per `context_id`, while
different contexts in the same HMG can be retained as WGD copies. TEs without
a strong context still use the manifest fallback quota.

The fallback quota cannot merge two distinct TEs on one genome/contig,
including adjacent, overlapping or nested annotations. That conflict is
lifted only for different unique strong contexts in the same known HMG and
explicit WGD node. Explicit WGD pairs that fail locus membership gates remain
`UNKNOWN`; proximity alone cannot relabel them as tandem copies.

## Claimability invariants

| Claim | Required evidence |
|---|---|
| annotated presence/discordance | callable paired flanks, established mapping confidence and eligible overlapping annotation |
| unannotated sequence presence | above plus ≥50% TE-body alignment and local identity ≥0.50 |
| confirmed empty site | callable paired flanks, established mapping, local identity ≥0.50, insertion fraction ≥0.70 and no excessive target gap |
| structural alternative | callable evidence and local identity ≥0.50 |

MCScanX, missing annotation, aggregate delta identity and a shortened projected
span cannot independently produce `EMPTY_SITE_CONFIRMED`.

## Integrity validation

`tests/validate_schema.py PREFIX` validates schema versions, primary and foreign
keys, provider-specific missing semantics, score normalization, logits,
prediction sets, matching/solver state, context chains, ambiguity links,
edge-group deduplication, row counts and the zero-FASTA-reopen invariant. It is
executed by `make check` and under ASan/UBSan by `make asan`.

Schema 1.2.0 prediction IDs are stable within their documented evidence and
input contract, but they are not valid primary keys for independent truth.
Publication benchmarking joins `(dataset_id, genome_id, te_id)` semantic nodes
and truth locus IDs as specified in [BENCHMARK.md](BENCHMARK.md). Generated
`CAN`, `EVD`, `DEC`, `TEL` and `INS` values are used only as foreign keys inside
one prediction bundle.
