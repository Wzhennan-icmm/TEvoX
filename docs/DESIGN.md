# TEvoX 0.3 correctness model

This document freezes the algorithm implemented by `0.3.0-alpha.1`. All
intervals are zero-based and half-open. Schema `1.0.0` separates observations,
inference decisions and biological claims.

## Evidence hierarchy

One native PAF record defines an `evidence_group_id`. Its automatic reverse
traversal belongs to the same group and is explicitly dependent. It can provide
the reverse coordinate view but cannot create independent reciprocal support.

For each `source TE → target genome`, TEvoX retains every spanning alignment as
a directed evidence observation. Each observation retains every nearby target
TE as a candidate. Candidates need at least 50% reciprocal interval overlap to
be eligible; proximity, family similarity or a permissive copy quota cannot
make a zero-overlap adjacent TE eligible.

All observations within five quality points of the best are compared. If their
target contigs, coordinates, target annotations or biological interpretations
disagree, the source-level decision is `AMBIGUOUS`. Every triggering
`evidence_id` is preserved in `linked_evidence_ids`; only retaining the winner
would violate schema `1.0.0`.

## Observed features and missing values

The directed observation measures:

- full-length left and right flank coverage;
- local identity among aligned query bases from `cs:Z` or exact `=X` CIGAR
  operations (query insertions are measured separately, not counted as
  mismatches);
- TE-body aligned fraction and query-insertion fraction;
- observed MAPQ, with PAF 255 represented as missing;
- projected coordinates and target ambiguous-base fraction;
- nearby candidate overlap, boundary distance and family relation.

A flank truncated by a source contig boundary is missing, not 1.0. Ordinary
`M` operations do not identify matches versus mismatches, so local identity is
missing rather than estimated from the whole PAF record. Ranking scores are
normalized over observed terms and multiplied by evidence completeness; the
score is not a probability.

## Three-axis state and claimability

Internal state is factored into:

```text
Technical:   CALLABLE | GAP | AMBIGUOUS | UNCALLABLE
Biological:  PRESENT | EMPTY | STRUCTURAL_ALTERNATIVE | UNKNOWN
Annotation:  MATCHED | MISSING | FAMILY_CONFLICT | NOT_APPLICABLE | UNKNOWN
```

`claimable` is evaluated after these axes. Technical gap, ambiguity or
uncallability always blocks a biological claim. Annotated presence can be
claimed from callable paired-flank geometry and an eligible overlapping target
annotation. Unannotated presence, empty-site and structural-alternative claims
additionally require observed local identity of at least 0.50. A confirmed
empty site also requires at least 70% query-insertion coverage of the source TE.

The legacy eight-state field is a projection of these axes and the claim gate.
In particular, `BIO=EMPTY` with `claimable=false` maps to `UNCALLABLE`, not
`EMPTY_SITE_CONFIRMED`.

## Global locus graph

Only a callable, unambiguous decision with an eligible selected annotation
candidate can support an edge. Independent reciprocal support adds five points
only when opposite directions use different evidence groups.

Edges are processed in deterministic score and node order. A component merge
must:

1. exceed `min_edge_score`;
2. have no direct known-family conflict;
3. avoid a known-family conflict across the two complete components, preventing
   an unknown-family node from bridging incompatible families;
4. respect the legacy per-genome component quota.

Every rejected edge records `selection_reason`. TE and PAF records are sorted
canonically so input row order does not change graph selection or locus IDs.

The component quota is not described as WGD-aware: no subgenome, haplotype,
syntenic-copy or WGD-node context exists in v0.3. That model belongs to v0.4.

## Phylogeny

The exploratory helper maps present states to 1, confirmed empty sites to 0,
and every technical/nonclaimable state to `{0,1}` before Fitch parsimony. It
marks events touching ambiguous nodes. Branch-length-aware stochastic Dollo
inference is future work.
