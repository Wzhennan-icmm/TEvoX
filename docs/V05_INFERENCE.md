# v0.5 inference contract

TEvoX `0.5.0-alpha.2` retains a missing-aware score layer and constrained global
optimizer. The built-in model is deliberately named
`BUILTIN_UNCALIBRATED_V1`: its normalized scores are **not calibrated
probabilities** and must not be reported as posterior probabilities. The
v0.4 three-axis state and claimability rules remain the authoritative safety
gate for presence, empty-site and structural claims.

## Candidate membership score

For a candidate TE pair, the built-in logit is

\[
\begin{aligned}
z={}&-5.5+2.4r+1.4b+0.4t\\
 &+F_{flank}+F_{identity}+F_{mapq}+F_N\\
 &+F_{family}+F_{context},
\end{aligned}
\]

where `r` is reciprocal TE overlap, `b` is the boundary score and `t` is
TE-body aligned fraction. Observed and missing terms are:

| Term | Observed contribution | Missing/unknown contribution |
|---|---:|---:|
| paired flank | `1.2 × min(left,right)` | `-1.2` |
| local identity | `1.0 × identity` | use aggregate term if available |
| aggregate identity fallback | `0.7 × identity - 0.15` | net `-0.05` if all identity is missing |
| MAPQ | `0.5 × min(MAPQ/60,1)` | `0` plus explicit missing flag |
| target N | `0.6 × (1-N)` | `-0.3` |
| family | match `+1.0`; conflict `-2.5` | `-0.15` |
| copy context | supported `+1.2`; conflict `-3.0`; ambiguous `-2.0` | `0` |

The reported membership score is `sigmoid(z)`. The feature table records an
observed mask and named missing fields, so missing measurements are never
silently converted to numeric zero. For a reduced edge, every additional
independent evidence group contributes `0.15 ln(1 + groups - 1)` and genuine
native reciprocal evidence contributes `0.35`. A generated reverse view never
qualifies as independent evidence.

The default graph gate is `membership_score >= 0.50`; values below 0.50 are not
accepted because the global objective uses positive log odds. `--min-membership`
may make this gate stricter, not weaker.

## Three-axis observation scores

`observation_scores.tsv` reports normalized energy scores for the technical,
biological and annotation axes. The v0.4 discrete state receives an anchor
term, then observable flank, identity, mapping, gap, TE-body and annotation
features modify the other energies. These scores expose uncertainty and
prediction sets without changing claimability. The row is marked
`out_of_domain=true` when fewer than three core evidence classes are observed
or mapping context is ambiguous.

## Block-level matching

Candidate reporting and graph inference have separate deterministic limits.
`--max-candidates` controls only candidate/feature/context table rows;
`--max-graph-candidates` (default 64, `0` for unlimited) controls which ranked
candidates can form edges. The latter is recorded because changing it can
change matching and loci; changing the former alone cannot.

An MCScanX-supported edge is assigned to the ordered pair of strong copy
contexts that bracket its two TEs. Edges first pass the legacy edge threshold,
known-family compatibility and positive membership threshold; a pre-gate
failure cannot occupy a matching slot. Within each context pair, TEvoX solves
a maximum-weight one-to-one assignment using the Hungarian algorithm and
positive membership log odds. This prevents two TEs in one copy slot from
competing into the same ancestral locus.

The exact path is used while both sides contain at most
`--exact-match-nodes` nodes (default 256). Larger blocks use a deterministic
weight/ID-ordered greedy fallback. Every edge reports its matching group,
method and selection result; the fallback is never labelled optimal. Edges
without strong context are `NOT_APPLICABLE` and proceed to the global solver
under the legacy quota and other hard constraints.

## Component optimization

For each connected candidate component, TEvoX maximizes

\[
\sum_{(i,j):\;L_i=L_j}\log\frac{s_{ij}}{1-s_{ij}},
\]

where `s` is the uncalibrated membership score and `L` is the reconstructed
locus assignment. A merge is feasible only if it:

1. has no incompatible known TE families;
2. contains at most one member per strong `context_id`;
3. does not mix distinct known homology groups;
4. respects the manifest quota for every member lacking strong context.

Distinct TEs from the same genome and contig—including adjacent, overlapping
and nested annotations—are an additional component hard negative. They can
coexist only when both have different unique strong contexts in the same
nonzero HMG, explicitly share one known WGD node and are not an allelic
haplotype pair. This is a narrow, evidence-backed WGD exception, not a
copy-quota shortcut. When such an explicit WGD pair remains in separate loci
after membership gating, physical proximity cannot produce a tandem label.

Components with at most `--exact-max-edges` active edges (default 18, maximum
24) are exhaustively enumerated. Their status is `OPTIMAL`, upper bound equals
the achieved objective and relative gap is zero. Larger components use a
deterministic greedy initializer. Its upper bound is the sum of all positive
edge weights, so the reported gap is valid but can be loose. The solver table
records method, status, objective, bound, gap and states explored per stable
component ID.

The exact implementation is a constrained partition enumeration, not an
external ILP solver. Calling the heuristic path “optimal”, or calling either
path a calibrated Bayesian posterior, is outside this contract.

`OPTIMAL` is conditional on graph-candidate pruning and the preceding
one-to-one block matching. The current solver does not jointly reconsider an
alternative rejected by Hungarian after a later component constraint blocks a
selected edge.

## Relation scores

`relations.tsv` contains normalized scores over:

```text
ORTHOLOG  WGD_HOMEOLOG  ALLELIC  TANDEM_PARALOG
SEGMENTAL_PARALOG  TRANSPOSED_PARALOG  UNKNOWN
```

Cross-genome members of one reconstructed locus are labelled `ORTHOLOG`.
Same-genome members in distinct contexts of one HMG can be labelled
`WGD_HOMEOLOG` when their WGD metadata agrees, or `ALLELIC` when distinct
haplotypes in a compatible subgenome are explicit. Nearby compatible-family
TEs outside a locus can be labelled `TANDEM_PARALOG`. A strong DNA edge that
conflicts with known synteny can receive a tentative transposed-paralog score,
but it is marked out of domain. TEvoX does not assert segmental or transposed
duplication from family strings alone.

Prediction sets accumulate normalized score mass up to
`--prediction-mass` (default 0.90). They express ambiguity in this fixed score
system, not frequentist coverage or Bayesian credible sets.

## Required calibration path

A future calibrated model must be trained on independent truth with
species-pair or clade-level holdout, then report Brier score, log loss,
calibration error and calibration curves. Until such a model is versioned and
validated, all v0.5 score tables and `run.json` retain
`calibration_status=UNCALIBRATED`.

Alpha.2 provides semantic multi-run evaluation, raw pre-decision feature
export and split leakage auditing for that future work. Their exact contracts
and non-claims are in [the benchmark contract](BENCHMARK.md).
