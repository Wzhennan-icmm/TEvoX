# Evidence schema 1.0.0

Schema `1.0.0` is the public contract for TEvoX `0.3.x`. Coordinates are
0-based half-open. A literal `.` means unobserved or not applicable; it never
means numeric zero.

## Identity and dependency keys

| Prefix | Entity | Uniqueness |
|---|---|---|
| `EVG` | native evidence group | one normalized native PAF record plus dependent views |
| `EVD` | directed observation | evidence group × source TE × target genome |
| `CAN` | annotation candidate | observation × target TE |
| `DEC` | source-level decision | source TE × target genome |
| `TEL` | reconstructed locus | selected global component |
| `INS` | locus instance | locus × genome |

IDs are deterministic hashes of normalized biological keys rather than input
row numbers. `provider_record` remains available for provenance but does not
control inference IDs.

## Tables

### `evidence.tsv`

One row per directed alignment observation. It contains provider path/record,
origin and dependency, source/alignment/projected coordinates, observed
features, missing-status fields, selected annotation, three-axis state,
claimability, quality/completeness, rank and decision codes.

`NATIVE` observations have dependency `INDEPENDENT`. An automatically inverted
view has origin `DERIVED_REVERSE`, dependency `DERIVED_SAME_GROUP`, and the same
`evidence_group_id` as its native record.

### `candidates.tsv`

One row for every target annotation within `candidate_window` of a projection.
`reciprocal_overlap` is:

\[
\min\left(\frac{|P\cap T|}{|P|},\frac{|P\cap T|}{|T|}\right).
\]

An ineligible candidate remains in the table with an explicit `decision_code`.
Only one candidate per observation may have `selected=true`.

`local_identity` is exact match identity among aligned query bases in the local
window. Query insertions are excluded from its denominator and reported
separately as `insertion_fraction`; this prevents a true long TE insertion from
artificially lowering the identity of its flanking empty site.

### `decisions.tsv`

One row per source TE × target genome. `winner_evidence_id` identifies the
deterministic ranking winner. `linked_evidence_ids` lists every near-best
observation used to test consistency. Its cardinality must equal
`near_best_count`, including ambiguous decisions.

### `edges.tsv`

One row per undirected annotated TE pair. `evidence_ids` lists directed support;
`evidence_group_ids` is deduplicated. `support_count` equals the number of
unique evidence groups. `independent_reciprocal=true` requires opposite
directions from distinct groups. `selection_reason` records global acceptance
or rejection.

### `loci.tsv` and `instances.tsv`

`loci.tsv` describes selected components. `instances.tsv` is the primary
biological output at locus × genome grain, with the three state axes,
`claim_type`, `claimable`, `claimability_reason`, supporting decision/evidence
keys and an aggregation code. A target annotation assigned to a different
global component cannot silently support presence; it becomes
`TARGET_ASSIGNED_TO_DIFFERENT_LOCUS`.

### Compatibility and provenance files

- `states.tsv` retains the legacy eight-state matrix used by the phylogeny
  helper and appends v0.3 axes, claim fields and external keys.
- `summary.tsv` counts legacy states by genome.
- `run.json` records software/schema versions, coordinate convention,
  parameters, input paths and table row counts.

## Claimability rules

| Proposed claim | Minimum v0.3 gate |
|---|---|
| annotated presence | callable full paired flanks, observed passing MAPQ, eligible ≥50% reciprocal-overlap annotation |
| unannotated presence | annotated-presence technical gates plus ≥50% TE-body alignment and observed local identity ≥0.50 |
| confirmed empty site | callable full paired flanks, observed passing MAPQ, callable target sequence, observed local identity ≥0.50, query-insertion fraction ≥0.70 |
| structural alternative | callable evidence and observed local identity ≥0.50 |
| annotation discordance | callable eligible annotation with known family conflict |

MCScanX/synteny priors, missing annotation alone and a shortened projected span
cannot produce `EMPTY_SITE_CONFIRMED`.

## Integrity validation

`tests/validate_schema.py PREFIX` verifies schema versions, primary/foreign
keys, missing identity and MAPQ semantics, ambiguity links, deduplicated edge
groups and `run.json` row counts. It runs as part of `make check`.
