#include "tevox.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    TvState state;
    TvTechnicalState technical_state;
    TvBiologicalState biological_state;
    TvAnnotationState annotation_state;
    TvClaimType claim_type;
    bool claimable;
    double quality;
    double evidence_completeness;
    int copies;
    const TvDecision *decision;
    char decision_code[96];
    char claimability_reason[96];
} InstanceCall;

static FILE *open_output(const char *prefix, const char *suffix, char **path)
{
    size_t length = strlen(prefix) + strlen(suffix) + 1;
    FILE *stream;

    *path = malloc(length);
    if (*path == NULL) {
        tv_print_error("out of memory");
        return NULL;
    }
    (void)snprintf(*path, length, "%s%s", prefix, suffix);
    stream = fopen(*path, "w");
    if (stream == NULL) {
        tv_print_error("cannot write '%s': %s", *path, strerror(errno));
        free(*path);
        *path = NULL;
    }
    return stream;
}

static void print_bool(FILE *stream, bool value)
{
    fputs(value ? "true" : "false", stream);
}

static void print_double(FILE *stream, double value, int precision)
{
    if (isnan(value)) {
        fputc('.', stream);
    } else {
        (void)fprintf(stream, "%.*f", precision, value);
    }
}

static void print_number(FILE *stream, double value)
{
    if (isnan(value)) {
        fputc('.', stream);
    } else {
        (void)fprintf(stream, "%.17g", value);
    }
}

static void print_text(FILE *stream, const char *value)
{
    fputs(value != NULL && *value != '\0' ? value : TEVOX_UNKNOWN, stream);
}

static void print_evidence_id(FILE *stream, uint64_t id)
{
    (void)fprintf(stream, "EVD%016llx", (unsigned long long)id);
}

static void print_group_id(FILE *stream, uint64_t id)
{
    (void)fprintf(stream, "EVG%016llx", (unsigned long long)id);
}

static void print_decision_id(FILE *stream, uint64_t id)
{
    (void)fprintf(stream, "DEC%016llx", (unsigned long long)id);
}

static void print_candidate_id(FILE *stream, uint64_t id)
{
    (void)fprintf(stream, "CAN%016llx", (unsigned long long)id);
}

static void print_block_id(FILE *stream, uint64_t id)
{
    (void)fprintf(stream, "SBL%016llx", (unsigned long long)id);
}

static void print_anchor_id(FILE *stream, uint64_t id)
{
    (void)fprintf(stream, "SYN%016llx", (unsigned long long)id);
}

static void print_context_id(FILE *stream, uint64_t id)
{
    (void)fprintf(stream, "CTX%016llx", (unsigned long long)id);
}

static void print_homology_group_id(FILE *stream, uint64_t id)
{
    if (id == 0) {
        fputc('.', stream);
    } else {
        (void)fprintf(stream, "HMG%016llx", (unsigned long long)id);
    }
}

static void print_te_context_id(FILE *stream, uint64_t id)
{
    (void)fprintf(stream, "TEC%016llx", (unsigned long long)id);
}

static bool decision_disagrees(const TvRun *run, const TvDecision *left,
                               const TvDecision *right)
{
    const TvProjection *a;
    const TvProjection *b;

    if (left->winner_projection < 0 || right->winner_projection < 0) {
        return left->winner_projection != right->winner_projection;
    }
    if (left->state != right->state) {
        return true;
    }
    a = &run->projections[left->winner_projection];
    b = &run->projections[right->winner_projection];
    if (strcmp(a->contig, b->contig) != 0) {
        return true;
    }
    if (llabs(a->start - b->start) > run->cfg.candidate_window
        || llabs(a->end - b->end) > run->cfg.candidate_window) {
        return true;
    }
    if (a->target_te != b->target_te) {
        return true;
    }
    return false;
}

static InstanceCall instance_call(TvRun *run, int locus, int genome)
{
    InstanceCall call = {
        .state = TV_UNCALLABLE,
        .technical_state = TV_TECH_UNCALLABLE,
        .biological_state = TV_BIO_UNKNOWN,
        .annotation_state = TV_ANN_UNKNOWN,
        .claim_type = TV_CLAIM_NONE
    };
    double best_edge_score = 0.0;

    for (size_t node = 0; node < run->n_nodes; node++) {
        if (run->components[node] == locus
            && run->nodes[node].genome == genome) {
            call.copies++;
        }
    }
    if (call.copies > 0) {
        call.state = TV_PRESENT_ANNOTATED;
        call.technical_state = TV_TECH_CALLABLE;
        call.biological_state = TV_BIO_PRESENT;
        call.annotation_state = TV_ANN_MATCHED;
        call.claim_type = TV_CLAIM_PRESENCE;
        call.claimable = true;
        call.evidence_completeness = 1.0;
        for (size_t edge_index = 0; edge_index < run->n_edges; edge_index++) {
            TvEdge *edge = &run->edges[edge_index];

            if (edge->selected
                && (run->components[edge->a] == locus
                    || run->components[edge->b] == locus)
                && edge->score > best_edge_score) {
                best_edge_score = edge->score;
            }
        }
        call.quality = best_edge_score > 0.0
            ? 70.0 + 0.30 * best_edge_score : 70.0;
        (void)snprintf(call.decision_code, sizeof(call.decision_code),
                       "ANNOTATED_COMPONENT");
        (void)snprintf(call.claimability_reason,
                       sizeof(call.claimability_reason),
                       "DIRECT_ANNOTATION_PRESENT");
        return call;
    }
    for (size_t index = 0; index < run->n_decisions; index++) {
        TvDecision *decision = &run->decisions[index];

        if (decision->target_genome != genome
            || run->components[decision->source_te] != locus) {
            continue;
        }
        if (call.decision == NULL || decision->quality > call.decision->quality
            || (decision->quality == call.decision->quality
                && decision->decision_id < call.decision->decision_id)) {
            call.decision = decision;
        }
    }
    if (call.decision == NULL) {
        (void)snprintf(call.decision_code, sizeof(call.decision_code),
                       "NO_EVIDENCE");
        (void)snprintf(call.claimability_reason,
                       sizeof(call.claimability_reason), "NO_EVIDENCE");
        return call;
    }
    call.state = call.decision->state;
    call.technical_state = call.decision->technical_state;
    call.biological_state = call.decision->biological_state;
    call.annotation_state = call.decision->annotation_state;
    call.claim_type = call.decision->claim_type;
    call.claimable = call.decision->claimable;
    call.quality = call.decision->quality;
    call.evidence_completeness = call.decision->evidence_completeness;
    (void)snprintf(call.decision_code, sizeof(call.decision_code), "%s",
                   call.decision->decision_code);
    (void)snprintf(call.claimability_reason,
                   sizeof(call.claimability_reason), "%s",
                   call.decision->claimability_reason);
    if (call.decision->winner_projection >= 0) {
        const TvProjection *winner =
            &run->projections[call.decision->winner_projection];

        if (call.state == TV_PRESENT_ANNOTATED
            && winner->target_te >= 0
            && run->components[winner->target_te] != locus) {
            call.state = TV_PROJECTION_AMBIGUOUS;
            call.technical_state = TV_TECH_AMBIGUOUS;
            call.biological_state = TV_BIO_UNKNOWN;
            call.annotation_state = TV_ANN_UNKNOWN;
            call.claim_type = TV_CLAIM_NONE;
            call.claimable = false;
            call.quality *= 0.75;
            (void)snprintf(call.decision_code, sizeof(call.decision_code),
                           "TARGET_ASSIGNED_TO_DIFFERENT_LOCUS");
            (void)snprintf(call.claimability_reason,
                           sizeof(call.claimability_reason),
                           "GLOBAL_LOCUS_CONFLICT");
        }
    }
    bool component_decisions_conflict = false;
    for (size_t index = 0; index < run->n_decisions; index++) {
        TvDecision *decision = &run->decisions[index];

        if (decision != call.decision
            && decision->target_genome == genome
            && run->components[decision->source_te] == locus
            && decision->quality
               >= call.decision->quality - run->cfg.near_best_delta
            && decision_disagrees(run, call.decision, decision)) {
            component_decisions_conflict = true;
            break;
        }
    }
    if (component_decisions_conflict) {
        call.state = TV_PROJECTION_AMBIGUOUS;
        call.technical_state = TV_TECH_AMBIGUOUS;
        call.biological_state = TV_BIO_UNKNOWN;
        call.annotation_state = TV_ANN_UNKNOWN;
        call.claim_type = TV_CLAIM_NONE;
        call.claimable = false;
        call.quality *= 0.75;
        (void)snprintf(call.decision_code, sizeof(call.decision_code),
                       "CONFLICTING_COMPONENT_DECISIONS");
        (void)snprintf(call.claimability_reason,
                       sizeof(call.claimability_reason),
                       "TECHNICAL_AMBIGUITY");
    }
    return call;
}

static void print_locus_members(FILE *stream, TvRun *run, int locus, int field)
{
    bool first = true;

    for (size_t index = 0; index < run->n_nodes; index++) {
        TvTE *te;

        if (run->components[index] != locus) {
            continue;
        }
        if (!first) {
            fputc(',', stream);
        }
        first = false;
        te = &run->nodes[index];
        if (field == 0) {
            (void)fprintf(stream, "%s:%s", run->genomes[te->genome].id,
                          te->id);
        } else if (field == 1) {
            fputs(te->family, stream);
        } else {
            (void)fprintf(stream, "%s:%lld-%lld", te->contig,
                          (long long)te->start, (long long)te->end);
        }
    }
    if (first) {
        fputc('.', stream);
    }
}

static void print_present(FILE *stream, TvRun *run, int locus, int genome,
                          int field)
{
    bool first = true;

    for (size_t index = 0; index < run->n_nodes; index++) {
        TvTE *te = &run->nodes[index];

        if (run->components[index] != locus || te->genome != genome) {
            continue;
        }
        if (!first) {
            fputc(',', stream);
        }
        first = false;
        if (field == 0) {
            fputs(te->id, stream);
        } else if (field == 1) {
            fputs(te->contig, stream);
        } else if (field == 2) {
            (void)fprintf(stream, "%lld", (long long)te->start);
        } else {
            (void)fprintf(stream, "%lld", (long long)te->end);
        }
    }
    if (first) {
        fputc('.', stream);
    }
}

static void print_decision_evidence_ids(FILE *stream, const TvRun *run,
                                        const TvDecision *decision)
{
    bool first = true;
    size_t end = decision->projection_start + decision->projection_count;

    for (size_t index = decision->projection_start; index < end; index++) {
        const TvProjection *projection = &run->projections[index];

        if (!projection->near_best) {
            continue;
        }
        if (!first) {
            fputc(',', stream);
        }
        first = false;
        print_evidence_id(stream, projection->evidence_id);
    }
    if (first) {
        fputc('.', stream);
    }
}

static void print_instance_decisions(FILE *stream, TvRun *run, int locus,
                                     int genome, bool evidence)
{
    bool first = true;

    for (size_t index = 0; index < run->n_decisions; index++) {
        TvDecision *decision = &run->decisions[index];

        if (decision->target_genome != genome
            || run->components[decision->source_te] != locus) {
            continue;
        }
        if (!evidence) {
            if (!first) {
                fputc(',', stream);
            }
            first = false;
            print_decision_id(stream, decision->decision_id);
        } else {
            size_t end = decision->projection_start + decision->projection_count;

            for (size_t projection_index = decision->projection_start;
                 projection_index < end; projection_index++) {
                TvProjection *projection = &run->projections[projection_index];

                if (!projection->near_best) {
                    continue;
                }
                if (!first) {
                    fputc(',', stream);
                }
                first = false;
                print_evidence_id(stream, projection->evidence_id);
            }
        }
    }
    if (first) {
        fputc('.', stream);
    }
}

static void print_edge_support(FILE *stream, const TvRun *run,
                               const TvEdge *edge, bool groups)
{
    bool first = true;
    uint64_t previous = 0;
    bool have_previous = false;
    size_t end = edge->support_start + edge->support_record_count;

    for (size_t index = edge->support_start; index < end; index++) {
        const TvEdgeSupport *support = &run->edge_support[index];
        uint64_t value = groups ? support->evidence_group_id
                                : support->evidence_id;

        if (have_previous && value == previous) {
            continue;
        }
        if (!first) {
            fputc(',', stream);
        }
        first = false;
        if (groups) {
            print_group_id(stream, support->evidence_group_id);
        } else {
            print_evidence_id(stream, support->evidence_id);
        }
        previous = value;
        have_previous = true;
    }
    if (first) {
        fputc('.', stream);
    }
}

static int write_evidence(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".evidence.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tevidence_id\tevidence_group_id\tdecision_id\tprovider"
          "\tprovider_path\tprovider_record\tprovider_line\torigin\tdependency"
          "\tquery_genome_id"
          "\ttarget_genome_id\tsource_te_id\tsource_contig\tsource_start\tsource_end"
          "\talignment_query_contig\talignment_query_start\talignment_query_end"
          "\talignment_target_contig\talignment_target_start\talignment_target_end"
          "\tprojection_contig\tprojection_start\tprojection_end"
          "\talignment_identity\talignment_identity_method\talignment_error_count"
          "\tsimilarity_error_count\tnonalpha_count\tlocal_identity"
          "\tidentity_method\tmapq\tmapq_status\tmapping_confidence"
          "\tleft_flank\tright_flank"
          "\tleft_flank_status\tright_flank_status\tte_aligned_fraction"
          "\tinsertion_fraction\ttarget_n_fraction\tnearby_candidate_count"
          "\tretained_candidate_count\teligible_candidate_count"
          "\tselected_target_te_id\ttechnical_state"
          "\tbiological_state\tannotation_state\tlegacy_state\tclaim_type"
          "\tclaimable\tquality\tevidence_completeness\tobservation_rank\tprimary"
          "\tnear_best\tdecision_ambiguous\tdecision_code\tclaimability_reason\n",
          stream);
    for (size_t index = 0; index < run->n_projections; index++) {
        TvProjection *projection = &run->projections[index];
        TvPaf *paf = &run->pafs[projection->paf_index];
        TvTE *source = &run->nodes[projection->source_te];
        TvDecision *decision = &run->decisions[projection->decision_index];

        (void)fprintf(stream, "%s\t", TEVOX_SCHEMA_VERSION);
        print_evidence_id(stream, projection->evidence_id);
        fputc('\t', stream);
        print_group_id(stream, paf->evidence_group_id);
        fputc('\t', stream);
        print_decision_id(stream, decision->decision_id);
        (void)fprintf(stream, "\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s"
                      "\t%lld\t%lld\t%s\t%lld\t%lld\t%s\t%lld\t%lld"
                      "\t%s\t%lld\t%lld\t",
                      tv_alignment_provider_name(paf->provider),
                      paf->source_path, paf->source_record, paf->source_line,
                      tv_evidence_origin_name(paf->origin),
                      paf->origin == TV_EVIDENCE_NATIVE
                          ? "INDEPENDENT" : "DERIVED_SAME_GROUP",
                      run->genomes[paf->query_genome].id,
                      run->genomes[paf->target_genome].id,
                      source->id, source->contig,
                      (long long)source->start, (long long)source->end,
                      paf->qname, (long long)paf->qstart,
                      (long long)paf->qend, paf->tname,
                      (long long)paf->tstart, (long long)paf->tend,
                      projection->contig, (long long)projection->start,
                      (long long)projection->end);
        print_double(stream, paf->aggregate_identity, 6);
        (void)fprintf(stream, "\t%s\t%lld\t",
                      tv_aggregate_identity_method_name(
                          paf->aggregate_identity_method),
                      (long long)paf->error_count);
        if (paf->provider == TV_ALIGNMENT_MUMMER_DELTA) {
            (void)fprintf(stream, "%lld\t%lld\t",
                          (long long)paf->similarity_error_count,
                          (long long)paf->nonalpha_count);
        } else {
            fputs(".\t.\t", stream);
        }
        print_double(stream, projection->identity, 6);
        (void)fprintf(stream, "\t%s\t",
                      tv_identity_method_name(paf->identity_method));
        if (projection->mapq_observed) {
            (void)fprintf(stream, "%d\t%s\t", projection->mapq,
                          tv_mapq_status_name(paf->mapq_status));
        } else {
            (void)fprintf(stream, ".\t%s\t",
                          tv_mapq_status_name(paf->mapq_status));
        }
        (void)fprintf(stream, "%s\t",
                      tv_mapping_confidence_name(
                          projection->mapping_confidence));
        print_double(stream, projection->left_flank, 6);
        fputc('\t', stream);
        print_double(stream, projection->right_flank, 6);
        (void)fprintf(stream, "\t%s\t%s\t%.6f\t%.6f\t",
                      projection->left_flank_observed ? "OBSERVED" : "CONTIG_EDGE",
                      projection->right_flank_observed ? "OBSERVED" : "CONTIG_EDGE",
                      projection->te_aligned_fraction,
                      projection->insertion_fraction);
        print_double(stream, projection->n_fraction, 6);
        (void)fprintf(stream, "\t%d\t%zu\t%d\t",
                      projection->nearby_candidate_count,
                      projection->candidate_count,
                      projection->eligible_candidate_count);
        if (projection->target_te >= 0) {
            fputs(run->nodes[projection->target_te].id, stream);
        } else {
            fputc('.', stream);
        }
        (void)fprintf(stream, "\t%s\t%s\t%s\t%s\t%s\t",
                      tv_technical_state_name(projection->technical_state),
                      tv_biological_state_name(projection->biological_state),
                      tv_annotation_state_name(projection->annotation_state),
                      tv_state_name(projection->state),
                      tv_claim_type_name(projection->claim_type));
        print_bool(stream, projection->claimable);
        (void)fprintf(stream, "\t%.2f\t%.6f\t%d\t", projection->quality,
                      projection->evidence_completeness, projection->rank);
        print_bool(stream, projection->primary);
        fputc('\t', stream);
        print_bool(stream, projection->near_best);
        fputc('\t', stream);
        print_bool(stream, projection->decision_ambiguous);
        (void)fprintf(stream, "\t%s\t%s\n", projection->decision_code,
                      projection->claimability_reason);
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_candidates(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".candidates.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tcandidate_id\tevidence_id\tdecision_id\ttarget_genome_id"
          "\ttarget_te_id\ttarget_contig\ttarget_start\ttarget_end\tscore"
          "\treciprocal_overlap\tboundary_score\tbreakpoint_distance"
          "\tfamily_relation\tcontext_relation\tshared_homology_group_id"
          "\tcontext_compatible\tcandidate_rank\teligible\tselected"
          "\tdecision_code\n", stream);
    for (size_t index = 0; index < run->n_candidates; index++) {
        TvCandidate *candidate = &run->candidates[index];
        TvProjection *projection = &run->projections[candidate->projection_index];
        TvDecision *decision = &run->decisions[projection->decision_index];
        TvTE *target = &run->nodes[candidate->target_te];

        (void)fprintf(stream, "%s\t", TEVOX_SCHEMA_VERSION);
        print_candidate_id(stream, candidate->candidate_id);
        fputc('\t', stream);
        print_evidence_id(stream, projection->evidence_id);
        fputc('\t', stream);
        print_decision_id(stream, decision->decision_id);
        (void)fprintf(stream, "\t%s\t%s\t%s\t%lld\t%lld\t%.2f\t%.6f\t%.6f"
                      "\t%d\t%s\t%s\t",
                      run->genomes[target->genome].id, target->id, target->contig,
                      (long long)target->start, (long long)target->end,
                      candidate->score, candidate->reciprocal_overlap,
                      candidate->boundary_score, candidate->breakpoint_distance,
                      tv_family_relation_name(candidate->family_relation),
                      tv_context_relation_name(candidate->context_relation));
        print_homology_group_id(stream,
                                candidate->shared_homology_group_id);
        fputc('\t', stream);
        print_bool(stream, candidate->context_compatible);
        (void)fprintf(stream, "\t%d\t", candidate->rank);
        print_bool(stream, candidate->eligible);
        fputc('\t', stream);
        print_bool(stream, candidate->selected);
        (void)fprintf(stream, "\t%s\n", candidate->decision_code);
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_synteny_blocks(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".synteny.blocks.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tblock_id\tevidence_group_id\thomology_group_id"
          "\tprovider\tsource_id\tprovider_path\tprovider_record"
          "\tprovider_line"
          "\tgenome_a\tcontig_a\tstart_a\tend_a\tgenome_b\tcontig_b"
          "\tstart_b\tend_b\torientation\tanchor_count\treported_score"
          "\treported_evalue\twgd_node\tcontext_a_id\tcontext_b_id"
          "\tstatus\treason\n", stream);
    for (size_t index = 0; index < run->n_synteny_blocks; index++) {
        const TvSyntenyBlock *block = &run->synteny_blocks[index];

        (void)fprintf(stream, "%s\t", TEVOX_SCHEMA_VERSION);
        print_block_id(stream, block->block_id);
        fputc('\t', stream);
        print_group_id(stream, block->evidence_group_id);
        fputc('\t', stream);
        print_homology_group_id(stream, block->homology_group_id);
        (void)fprintf(stream, "\tMCScanX\t%s\t%s\t%d\t%d\t%s\t%s\t%lld\t%lld"
                      "\t%s\t%s\t%lld\t%lld\t%c\t%d\t",
                      block->source_id, block->source_path,
                      block->source_record, block->source_line,
                      run->genomes[block->genome_a].id, block->contig_a,
                      (long long)block->start_a, (long long)block->end_a,
                      run->genomes[block->genome_b].id, block->contig_b,
                      (long long)block->start_b, (long long)block->end_b,
                      block->orientation, block->anchor_count);
        print_number(stream, block->reported_score);
        fputc('\t', stream);
        print_number(stream, block->reported_evalue);
        fputc('\t', stream);
        print_text(stream, block->wgd_node);
        fputc('\t', stream);
        print_context_id(stream, run->contexts[block->context_a].context_id);
        fputc('\t', stream);
        print_context_id(stream, run->contexts[block->context_b].context_id);
        fputc('\t', stream);
        print_text(stream, block->status);
        fputc('\t', stream);
        print_text(stream, block->reason);
        fputc('\n', stream);
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_synteny_anchors(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".synteny.anchors.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tanchor_id\tblock_id\tanchor_rank"
          "\tprovider_anchor_rank\tgene_a_id\tgenome_a\tgene_b_id"
          "\tgenome_b\treported_evalue\n", stream);
    for (size_t index = 0; index < run->n_synteny_anchors; index++) {
        const TvSyntenyAnchor *anchor = &run->synteny_anchors[index];
        const TvSyntenyBlock *block =
            &run->synteny_blocks[anchor->block_index];
        const TvGene *gene_a = &run->genes[anchor->gene_a];
        const TvGene *gene_b = &run->genes[anchor->gene_b];

        (void)fprintf(stream, "%s\t", TEVOX_SCHEMA_VERSION);
        print_anchor_id(stream, anchor->anchor_id);
        fputc('\t', stream);
        print_block_id(stream, block->block_id);
        (void)fprintf(stream, "\t%d\t%d\t%s\t%s\t%s\t%s\t", anchor->rank,
                      anchor->provider_rank,
                      gene_a->id, run->genomes[gene_a->genome].id,
                      gene_b->id, run->genomes[gene_b->genome].id);
        print_number(stream, anchor->reported_evalue);
        fputc('\n', stream);
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_contexts(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".contexts.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tcontext_id\thomology_group_id\tgenome_id\tcontig"
          "\tstart\tend\tsubgenome_id\thaplotype_id\tsyntenic_copy_id"
          "\twgd_node\tanchor_count\tstatus\n", stream);
    for (size_t index = 0; index < run->n_contexts; index++) {
        const TvCopyContext *context = &run->contexts[index];

        (void)fprintf(stream, "%s\t", TEVOX_SCHEMA_VERSION);
        print_context_id(stream, context->context_id);
        fputc('\t', stream);
        print_homology_group_id(stream, context->homology_group_id);
        (void)fprintf(stream, "\t%s\t%s\t%lld\t%lld\t",
                      run->genomes[context->genome].id, context->contig,
                      (long long)context->start, (long long)context->end);
        print_text(stream, context->subgenome_id);
        fputc('\t', stream);
        print_text(stream, context->haplotype_id);
        fputc('\t', stream);
        print_text(stream, context->syntenic_copy_id);
        fputc('\t', stream);
        print_text(stream, context->wgd_node);
        (void)fprintf(stream, "\t%d\t", context->anchor_count);
        print_text(stream, context->status);
        fputc('\n', stream);
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_te_contexts(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".te_contexts.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tte_context_id\tgenome_id\tte_id\tcontext_id"
          "\thomology_group_id\tassignment\toverlap_fraction"
          "\tleft_anchor_gene_id\tright_anchor_gene_id\n", stream);
    for (size_t index = 0; index < run->n_te_contexts; index++) {
        const TvTEContext *assignment = &run->te_contexts[index];
        const TvTE *te = &run->nodes[assignment->te_node];
        const TvCopyContext *context =
            &run->contexts[assignment->context_index];

        (void)fprintf(stream, "%s\t", TEVOX_SCHEMA_VERSION);
        print_te_context_id(stream, assignment->te_context_id);
        (void)fprintf(stream, "\t%s\t%s\t",
                      run->genomes[te->genome].id, te->id);
        print_context_id(stream, context->context_id);
        fputc('\t', stream);
        print_homology_group_id(stream, context->homology_group_id);
        (void)fprintf(stream, "\t%s\t%.6f\t",
                      tv_context_assignment_name(assignment->assignment),
                      assignment->overlap_fraction);
        if (assignment->left_anchor >= 0) {
            fputs(run->genes[assignment->left_anchor].id, stream);
        } else {
            fputc('.', stream);
        }
        fputc('\t', stream);
        if (assignment->right_anchor >= 0) {
            fputs(run->genes[assignment->right_anchor].id, stream);
        } else {
            fputc('.', stream);
        }
        fputc('\n', stream);
    }
    fclose(stream);
    free(path);
    return 0;
}

static void print_node_contexts(FILE *stream, const TvRun *run, int node)
{
    bool first = true;
    size_t begin = run->te_context_offsets == NULL
        ? 0 : run->te_context_offsets[node];
    size_t end = run->te_context_offsets == NULL
        ? 0 : run->te_context_offsets[node + 1];

    for (size_t index = begin; index < end; index++) {
        const TvTEContext *assignment = &run->te_contexts[index];

        if (!first) {
            fputc(',', stream);
        }
        first = false;
        print_context_id(stream,
                         run->contexts[assignment->context_index].context_id);
    }
    if (first) {
        fputc('.', stream);
    }
}

static int write_candidate_contexts(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".candidate_contexts.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tcandidate_id\tevidence_id\tsource_te_id"
          "\ttarget_te_id\tcontext_relation\tcontext_compatible"
          "\tshared_homology_group_id\tsource_context_ids"
          "\ttarget_context_ids\n", stream);
    for (size_t index = 0; index < run->n_candidates; index++) {
        const TvCandidate *candidate = &run->candidates[index];
        const TvProjection *projection =
            &run->projections[candidate->projection_index];

        (void)fprintf(stream, "%s\t", TEVOX_SCHEMA_VERSION);
        print_candidate_id(stream, candidate->candidate_id);
        fputc('\t', stream);
        print_evidence_id(stream, projection->evidence_id);
        (void)fprintf(stream, "\t%s\t%s\t%s\t",
                      run->nodes[projection->source_te].id,
                      run->nodes[candidate->target_te].id,
                      tv_context_relation_name(candidate->context_relation));
        print_bool(stream, candidate->context_compatible);
        fputc('\t', stream);
        print_homology_group_id(stream,
                                candidate->shared_homology_group_id);
        fputc('\t', stream);
        print_node_contexts(stream, run, projection->source_te);
        fputc('\t', stream);
        print_node_contexts(stream, run, candidate->target_te);
        fputc('\n', stream);
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_decisions(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".decisions.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tdecision_id\tsource_genome_id\tsource_te_id"
          "\ttarget_genome_id\ttechnical_state\tbiological_state\tannotation_state"
          "\tlegacy_state\tclaim_type\tclaimable\tquality\tevidence_completeness"
          "\tobservation_count\tnear_best_count\twinner_evidence_id"
          "\tlinked_evidence_ids\tdecision_code\tclaimability_reason\n", stream);
    for (size_t index = 0; index < run->n_decisions; index++) {
        TvDecision *decision = &run->decisions[index];
        TvTE *source = &run->nodes[decision->source_te];

        (void)fprintf(stream, "%s\t", TEVOX_SCHEMA_VERSION);
        print_decision_id(stream, decision->decision_id);
        (void)fprintf(stream, "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
                      run->genomes[source->genome].id, source->id,
                      run->genomes[decision->target_genome].id,
                      tv_technical_state_name(decision->technical_state),
                      tv_biological_state_name(decision->biological_state),
                      tv_annotation_state_name(decision->annotation_state),
                      tv_state_name(decision->state),
                      tv_claim_type_name(decision->claim_type));
        print_bool(stream, decision->claimable);
        (void)fprintf(stream, "\t%.2f\t%.6f\t%zu\t%d\t", decision->quality,
                      decision->evidence_completeness,
                      decision->projection_count, decision->near_best_count);
        if (decision->winner_projection >= 0) {
            print_evidence_id(
                stream, run->projections[decision->winner_projection].evidence_id);
        } else {
            fputc('.', stream);
        }
        fputc('\t', stream);
        print_decision_evidence_ids(stream, run, decision);
        (void)fprintf(stream, "\t%s\t%s\n", decision->decision_code,
                      decision->claimability_reason);
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_edges(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".edges.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tte_a\tgenome_a\tte_b\tgenome_b\tscore"
          "\tindependent_reciprocal\tfamily_compatible\tbreakpoint_distance"
          "\tsupport_count\tevidence_ids\tevidence_group_ids\tselected"
          "\tselection_reason\n", stream);
    for (size_t index = 0; index < run->n_edges; index++) {
        TvEdge *edge = &run->edges[index];
        TvTE *a = &run->nodes[edge->a];
        TvTE *b = &run->nodes[edge->b];

        (void)fprintf(stream, "%s\t%s\t%s\t%s\t%s\t%.2f\t",
                      TEVOX_SCHEMA_VERSION, a->id,
                      run->genomes[a->genome].id, b->id,
                      run->genomes[b->genome].id, edge->score);
        print_bool(stream, edge->independent_reciprocal);
        fputc('\t', stream);
        print_bool(stream, edge->family_compatible);
        (void)fprintf(stream, "\t%d\t%zu\t", edge->breakpoint_distance,
                      edge->support_count);
        print_edge_support(stream, run, edge, false);
        fputc('\t', stream);
        print_edge_support(stream, run, edge, true);
        fputc('\t', stream);
        print_bool(stream, edge->selected);
        (void)fprintf(stream, "\t%s\n", edge->selection_reason);
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_loci(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".loci.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tlocus_id\tmember_count\tgenome_count\tmembers"
          "\tfamilies\tcoordinates\tquality\n", stream);
    for (int locus = 0; locus < run->n_loci; locus++) {
        int members = 0;
        int genome_count = 0;
        double quality = 0.0;

        for (size_t node = 0; node < run->n_nodes; node++) {
            if (run->components[node] == locus) {
                members++;
            }
        }
        for (size_t genome = 0; genome < run->n_genomes; genome++) {
            InstanceCall call = instance_call(run, locus, (int)genome);

            if (call.copies > 0) {
                genome_count++;
            }
            quality += call.quality;
        }
        (void)fprintf(stream, "%s\tTEL%06d\t%d\t%d\t",
                      TEVOX_SCHEMA_VERSION, locus + 1, members, genome_count);
        print_locus_members(stream, run, locus, 0);
        fputc('\t', stream);
        print_locus_members(stream, run, locus, 1);
        fputc('\t', stream);
        print_locus_members(stream, run, locus, 2);
        (void)fprintf(stream, "\t%.2f\n", quality / (double)run->n_genomes);
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_instances(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".instances.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("schema_version\tinstance_id\tlocus_id\tgenome_id\ttechnical_state"
          "\tbiological_state\tannotation_state\tlegacy_state\tclaim_type"
          "\tclaimable\tclaimability_reason\tquality\tevidence_completeness"
          "\tcopy_count\tmember_ids\tcontig\tstart\tend\tsupporting_decision_ids"
          "\tsupporting_evidence_ids\tdecision_code\n", stream);
    for (int locus = 0; locus < run->n_loci; locus++) {
        for (size_t genome = 0; genome < run->n_genomes; genome++) {
            InstanceCall call = instance_call(run, locus, (int)genome);
            const TvProjection *projection = call.decision != NULL
                && call.decision->winner_projection >= 0
                ? &run->projections[call.decision->winner_projection] : NULL;

            (void)fprintf(stream,
                          "%s\tINS%06d:%s\tTEL%06d\t%s\t%s\t%s\t%s\t%s\t%s\t",
                          TEVOX_SCHEMA_VERSION, locus + 1,
                          run->genomes[genome].id, locus + 1,
                          run->genomes[genome].id,
                          tv_technical_state_name(call.technical_state),
                          tv_biological_state_name(call.biological_state),
                          tv_annotation_state_name(call.annotation_state),
                          tv_state_name(call.state),
                          tv_claim_type_name(call.claim_type));
            print_bool(stream, call.claimable);
            (void)fprintf(stream, "\t%s\t%.2f\t%.6f\t%d\t",
                          call.claimability_reason, call.quality,
                          call.evidence_completeness, call.copies);
            if (call.copies > 0) {
                print_present(stream, run, locus, (int)genome, 0);
                fputc('\t', stream);
                print_present(stream, run, locus, (int)genome, 1);
                fputc('\t', stream);
                print_present(stream, run, locus, (int)genome, 2);
                fputc('\t', stream);
                print_present(stream, run, locus, (int)genome, 3);
            } else if (projection != NULL) {
                (void)fprintf(stream, ".\t%s\t%lld\t%lld", projection->contig,
                              (long long)projection->start,
                              (long long)projection->end);
            } else {
                fputs(".\t.\t.\t.", stream);
            }
            fputc('\t', stream);
            print_instance_decisions(stream, run, locus, (int)genome, false);
            fputc('\t', stream);
            print_instance_decisions(stream, run, locus, (int)genome, true);
            (void)fprintf(stream, "\t%s\n", call.decision_code);
        }
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_states(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".states.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("locus_id\tgenome_id\tstate\tquality\tcopy_count\tmember_ids\tcontig"
          "\tstart\tend\tidentity\tleft_flank\tright_flank\tn_fraction\tmapq"
          "\tevidence\ttechnical_state\tbiological_state\tannotation_state"
          "\tclaim_type\tclaimable\tclaimability_reason\tdecision_id"
          "\tevidence_ids\tschema_version\n", stream);
    for (int locus = 0; locus < run->n_loci; locus++) {
        for (size_t genome = 0; genome < run->n_genomes; genome++) {
            InstanceCall call = instance_call(run, locus, (int)genome);
            const TvProjection *projection = call.decision != NULL
                && call.decision->winner_projection >= 0
                ? &run->projections[call.decision->winner_projection] : NULL;

            (void)fprintf(stream, "TEL%06d\t%s\t%s\t%.2f\t%d\t", locus + 1,
                          run->genomes[genome].id, tv_state_name(call.state),
                          call.quality, call.copies);
            if (call.copies > 0) {
                print_present(stream, run, locus, (int)genome, 0);
                fputc('\t', stream);
                print_present(stream, run, locus, (int)genome, 1);
                fputc('\t', stream);
                print_present(stream, run, locus, (int)genome, 2);
                fputc('\t', stream);
                print_present(stream, run, locus, (int)genome, 3);
                fputs("\t.\t.\t.\t.\t.\tANNOTATED_COMPONENT", stream);
            } else if (projection != NULL) {
                (void)fprintf(stream, ".\t%s\t%lld\t%lld\t",
                              projection->contig, (long long)projection->start,
                              (long long)projection->end);
                print_double(stream, projection->identity, 4);
                fputc('\t', stream);
                print_double(stream, projection->left_flank, 4);
                fputc('\t', stream);
                print_double(stream, projection->right_flank, 4);
                fputc('\t', stream);
                print_double(stream, projection->n_fraction, 4);
                fputc('\t', stream);
                if (projection->mapq_observed) {
                    (void)fprintf(stream, "%d", projection->mapq);
                } else {
                    fputc('.', stream);
                }
                (void)fprintf(stream, "\t%s", call.decision_code);
            } else {
                fputs(".\t.\t.\t.\t.\t.\t.\t.\t.\tNO_EVIDENCE", stream);
            }
            (void)fprintf(stream, "\t%s\t%s\t%s\t%s\t",
                          tv_technical_state_name(call.technical_state),
                          tv_biological_state_name(call.biological_state),
                          tv_annotation_state_name(call.annotation_state),
                          tv_claim_type_name(call.claim_type));
            print_bool(stream, call.claimable);
            (void)fprintf(stream, "\t%s\t", call.claimability_reason);
            if (call.decision != NULL) {
                print_decision_id(stream, call.decision->decision_id);
            } else {
                fputc('.', stream);
            }
            fputc('\t', stream);
            print_instance_decisions(stream, run, locus, (int)genome, true);
            (void)fprintf(stream, "\t%s\n", TEVOX_SCHEMA_VERSION);
        }
    }
    fclose(stream);
    free(path);
    return 0;
}

static int write_summary(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".summary.tsv", &path);

    if (stream == NULL) {
        return -1;
    }
    fputs("genome_id\tstate\tcount\n", stream);
    for (size_t genome = 0; genome < run->n_genomes; genome++) {
        int counts[TV_UNCALLABLE + 1] = {0};

        for (int locus = 0; locus < run->n_loci; locus++) {
            InstanceCall call = instance_call(run, locus, (int)genome);

            counts[call.state]++;
        }
        for (int state = 0; state <= TV_UNCALLABLE; state++) {
            (void)fprintf(stream, "%s\t%s\t%d\n", run->genomes[genome].id,
                          tv_state_name((TvState)state), counts[state]);
        }
    }
    fclose(stream);
    free(path);
    return 0;
}

static void json_string(FILE *stream, const char *text)
{
    fputc('"', stream);
    while (*text != '\0') {
        unsigned char character = (unsigned char)*text++;

        if (character == '"' || character == '\\') {
            fputc('\\', stream);
            fputc(character, stream);
        } else if (character == '\n') {
            fputs("\\n", stream);
        } else if (character == '\r') {
            fputs("\\r", stream);
        } else if (character == '\t') {
            fputs("\\t", stream);
        } else if (character < 0x20) {
            (void)fprintf(stream, "\\u%04x", (unsigned int)character);
        } else {
            fputc(character, stream);
        }
    }
    fputc('"', stream);
}

static int write_run_json(TvRun *run, const char *prefix)
{
    char *path = NULL;
    FILE *stream = open_output(prefix, ".run.json", &path);
    size_t native_groups = 0;
    size_t candidate_observations = 0;

    for (size_t index = 0; index < run->n_pafs; index++) {
        if (run->pafs[index].origin == TV_EVIDENCE_NATIVE) {
            native_groups++;
        }
    }
    for (size_t index = 0; index < run->n_projections; index++) {
        int count = run->projections[index].nearby_candidate_count;

        if (count > 0 && (size_t)count <= SIZE_MAX - candidate_observations) {
            candidate_observations += (size_t)count;
        }
    }

    if (stream == NULL) {
        return -1;
    }
    (void)fprintf(stream,
                  "{\n  \"software\": \"TEvoX\",\n  \"version\": \"%s\",\n"
                  "  \"schema_version\": \"%s\",\n"
                  "  \"coordinate_system\": \"0-based-half-open\",\n"
                  "  \"config\": {\n"
                  "    \"flank\": %d,\n"
                  "    \"candidate_window\": %d,\n"
                  "    \"min_mapq\": %d,\n"
                  "    \"min_flank_fraction\": %.6f,\n"
                  "    \"min_edge_score\": %.6f,\n"
                  "    \"max_n_fraction\": %.6f,\n"
                  "    \"min_reciprocal_overlap\": %.6f,\n"
                  "    \"near_best_delta\": %.6f,\n"
                  "    \"min_delta_identity\": %.6f,\n"
                  "    \"max_candidates\": %d\n"
                  "  },\n  \"counts\": {\n"
                  "    \"genomes\": %zu,\n    \"tes\": %zu,\n"
                  "    \"alignment_views\": %zu,\n"
                  "    \"native_evidence_groups\": %zu,\n"
                  "    \"synteny_evidence_groups\": %zu,\n"
                  "    \"evidence_observations\": %zu,\n"
                  "    \"candidate_observations\": %zu,\n"
                  "    \"candidates\": %zu,\n    \"decisions\": %zu,\n"
                  "    \"genes\": %zu,\n    \"synteny_blocks\": %zu,\n"
                  "    \"synteny_anchors\": %zu,\n"
                  "    \"copy_contexts\": %zu,\n"
                  "    \"te_context_assignments\": %zu,\n"
                  "    \"edges\": %zu,\n    \"loci\": %d\n  },\n"
                  "  \"genomes\": [\n",
                  TEVOX_VERSION, TEVOX_SCHEMA_VERSION, run->cfg.flank,
                  run->cfg.candidate_window, run->cfg.min_mapq,
                  run->cfg.min_flank_fraction, run->cfg.min_edge_score,
                  run->cfg.max_n_fraction, run->cfg.min_reciprocal_overlap,
                  run->cfg.near_best_delta, run->cfg.min_delta_identity,
                  run->cfg.max_candidates, run->n_genomes, run->n_nodes,
                  run->n_pafs, native_groups, run->n_synteny_blocks,
                  run->n_projections, candidate_observations,
                  run->n_candidates, run->n_decisions, run->n_genes,
                  run->n_synteny_blocks, run->n_synteny_anchors,
                  run->n_contexts, run->n_te_contexts, run->n_edges,
                  run->n_loci);
    for (size_t index = 0; index < run->n_genomes; index++) {
        TvGenome *genome = &run->genomes[index];

        fputs("    {\"genome_id\": ", stream);
        json_string(stream, genome->id);
        fputs(", \"fasta\": ", stream);
        json_string(stream, genome->fasta_path);
        fputs(", \"te_annotation\": ", stream);
        json_string(stream, genome->te_path);
        (void)fprintf(stream, ", \"legacy_max_locus_copies\": %d}%s\n",
                      genome->max_locus_copies,
                      index + 1 == run->n_genomes ? "" : ",");
    }
    fputs("  ],\n  \"alignment_evidence\": [\n", stream);
    size_t emitted = 0;
    for (size_t index = 0; index < run->n_pafs; index++) {
        TvPaf *paf = &run->pafs[index];
        char group_id[24];

        if (paf->origin != TV_EVIDENCE_NATIVE) {
            continue;
        }
        (void)snprintf(group_id, sizeof(group_id), "EVG%016llx",
                       (unsigned long long)paf->evidence_group_id);
        fputs("    {\"evidence_group_id\": ", stream);
        json_string(stream, group_id);
        fputs(", \"provider\": ", stream);
        json_string(stream, tv_alignment_provider_name(paf->provider));
        fputs(", \"path\": ", stream);
        json_string(stream, paf->source_path);
        (void)fprintf(stream,
                      ", \"record\": %d, \"line\": %d, \"query_genome_id\": ",
                      paf->source_record, paf->source_line);
        json_string(stream, run->genomes[paf->query_genome].id);
        fputs(", \"target_genome_id\": ", stream);
        json_string(stream, run->genomes[paf->target_genome].id);
        emitted++;
        (void)fprintf(stream, "}%s\n", emitted == native_groups ? "" : ",");
    }
    fputs("  ],\n  \"synteny_evidence\": [\n", stream);
    for (size_t index = 0; index < run->n_synteny_blocks; index++) {
        const TvSyntenyBlock *block = &run->synteny_blocks[index];
        char group_id[24];
        char block_id[24];

        (void)snprintf(group_id, sizeof(group_id), "EVG%016llx",
                       (unsigned long long)block->evidence_group_id);
        (void)snprintf(block_id, sizeof(block_id), "SBL%016llx",
                       (unsigned long long)block->block_id);
        fputs("    {\"evidence_group_id\": ", stream);
        json_string(stream, group_id);
        fputs(", \"block_id\": ", stream);
        json_string(stream, block_id);
        fputs(", \"provider\": \"MCScanX\", \"source_id\": ", stream);
        json_string(stream, block->source_id);
        fputs(", \"path\": ", stream);
        json_string(stream, block->source_path);
        (void)fprintf(stream, ", \"record\": %d, \"line\": %d}%s\n",
                      block->source_record, block->source_line,
                      index + 1 == run->n_synteny_blocks ? "" : ",");
    }
    (void)fprintf(stream,
                  "  ],\n  \"performance\": {\n"
                  "    \"paf_interval_queries\": %llu,\n"
                  "    \"alignment_records_examined\": %llu,\n"
                  "    \"te_interval_queries\": %llu,\n"
                  "    \"te_records_examined\": %llu,\n"
                  "    \"gap_queries\": %llu,\n"
                  "    \"fasta_reopens_after_index\": %llu,\n"
                  "    \"edge_support_records\": %llu\n"
                  "  },\n",
                  (unsigned long long)run->performance.paf_interval_queries,
                  (unsigned long long)run->performance.paf_records_examined,
                  (unsigned long long)run->performance.te_interval_queries,
                  (unsigned long long)run->performance.te_records_examined,
                  (unsigned long long)run->performance.gap_queries,
                  (unsigned long long)run->performance.fasta_reopens_after_index,
                  (unsigned long long)run->performance.edge_support_records);
    fputs("  \"outputs\": [\"evidence.tsv\", \"candidates.tsv\", "
          "\"candidate_contexts.tsv\", \"decisions.tsv\", \"edges.tsv\", "
          "\"loci.tsv\", \"instances.tsv\", \"states.tsv\", "
          "\"summary.tsv\", \"synteny.blocks.tsv\", "
          "\"synteny.anchors.tsv\", \"contexts.tsv\", "
          "\"te_contexts.tsv\"]\n}\n", stream);
    fclose(stream);
    free(path);
    return 0;
}

int tv_write_outputs(TvRun *run, const char *prefix)
{
    if (prefix == NULL || *prefix == '\0') {
        return -1;
    }
    if (write_evidence(run, prefix) != 0
        || write_candidates(run, prefix) != 0
        || write_candidate_contexts(run, prefix) != 0
        || write_synteny_blocks(run, prefix) != 0
        || write_synteny_anchors(run, prefix) != 0
        || write_contexts(run, prefix) != 0
        || write_te_contexts(run, prefix) != 0
        || write_decisions(run, prefix) != 0
        || write_edges(run, prefix) != 0
        || write_loci(run, prefix) != 0
        || write_instances(run, prefix) != 0
        || write_states(run, prefix) != 0
        || write_summary(run, prefix) != 0
        || write_run_json(run, prefix) != 0) {
        return -1;
    }
    return 0;
}
