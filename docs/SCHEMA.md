# Evidence schema 1.1.0

Schema `1.1.0` is the public contract for TEvoX `0.4.x`. All coordinates are
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
tables. Eligibility, the best and second-best score, and ambiguity are computed
over all candidates before truncation.

### `candidates.tsv` and `candidate_contexts.tsv`

`candidates.tsv` contains retained target annotations ordered by the single
complete comparator:

1. eligible before ineligible;
2. decreasing score;
3. decreasing reciprocal overlap;
4. increasing breakpoint distance;
5. stable candidate ID.

Ranks are contiguous within an observation and exactly one eligible winner may
have `selected=true`. Reciprocal overlap is

\[
\min\left(\frac{|P\cap T|}{|P|},\frac{|P\cap T|}{|T|}\right).
\]

`candidate_contexts.tsv` is a one-to-one extension keyed by `candidate_id`. It
contains source/target context foreign keys, a shared HMG when present, and
`SUPPORTED`, `CONFLICT`, `AMBIGUOUS`, or `UNKNOWN`.

### `decisions.tsv`

One row per source TE × target genome. `winner_evidence_id` is the deterministic
observation winner. `linked_evidence_ids` lists every near-best observation
used for consistency testing, and its cardinality equals `near_best_count`.

### `edges.tsv`

One row per undirected annotated TE pair. Evidence IDs are listed once and
evidence groups are deduplicated. `support_count` equals the number of unique
alignment groups. `independent_reciprocal=true` requires native opposite
directions from different groups; an automatic reverse never qualifies.

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
  performance counters and the complete output inventory.

Strong copy-context constraints replace the legacy quota for assigned TEs:
one ancestral locus may contain at most one member per `context_id`, while
different contexts in the same HMG can be retained as WGD copies. TEs without
a strong context still use the manifest fallback quota.

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
keys, provider-specific missing semantics, context chains, ambiguity links,
edge-group deduplication, row counts and the zero-FASTA-reopen invariant. It is
executed by `make check` and under ASan/UBSan by `make asan`.
