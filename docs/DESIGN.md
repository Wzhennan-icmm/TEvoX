# TEvoX 0.2 locus model

This document freezes the P0/P1 algorithm implemented by `0.2.0-alpha.1`.
All intervals are zero-based and half-open internally. FASTA contig lengths are
authoritative; inconsistent annotation and PAF records are rejected.

## Directed evidence

For each `source TE -> target genome` hypothesis, the CIGAR projects both TE
boundaries and measures aligned left/right flanks, local identity, MAPQ, target
ambiguous-base fraction and nearby annotations. A supplied PAF is inverted;
reverse-strand operations are reversed and insertion/deletion are swapped.

Decision order is: low mapping/flanks → `UNCALLABLE`; high target N →
`ASSEMBLY_GAP`; near ties → `PROJECTION_AMBIGUOUS`; known family conflict →
`FAMILY_OR_BOUNDARY_DISCORDANCE`; compatible annotation →
`PRESENT_ANNOTATED`; ≥70% of the TE in a query insertion →
`EMPTY_SITE_CONFIRMED`; projected span ≥50% TE length →
`PRESENT_UNANNOTATED`; otherwise `STRUCTURAL_ALTERNATIVE`.

## Quality score

```text
100 * (0.35 * paired_flank + 0.25 * identity + 0.15 * min(MAPQ/60,1)
       + 0.15 * (1-target_N_fraction) + 0.10 * uniqueness)
```

Uncallable, assembly-gap, ambiguity and discordance calls receive documented
penalties. This is a ranking score, not a probability; calibration is P2.

## Global locus graph

Candidate edges combine flank support (25), identity (20), reciprocal interval
overlap (20), breakpoint proximity (15), MAPQ (10) and family agreement (10).
Reciprocal evidence adds 5. Edges are processed in deterministic descending
order. A component merge must exceed `min_edge_score`, have no known family
conflict, and respect every genome's `max_locus_copies` quota.

This quota-constrained maximum-spanning-forest heuristic supports WGD
co-orthologs while keeping every accepted/rejected edge auditable. It is not
claimed to be a probabilistic maximum-likelihood orthology solution.

## Phylogeny

The helper maps present states to 1, confirmed empty sites to 0 and all other
states to `{0,1}` before Fitch parsimony. It emits assigned and possible states
and marks events touching ambiguous nodes.
