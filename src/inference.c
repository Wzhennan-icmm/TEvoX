#include "tevox.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

typedef struct {
    int edge_index;
    int left_node;
    int right_node;
    int left_context;
    int right_context;
    uint64_t group_id;
    uint64_t edge_id;
} TvMatchItem;

typedef struct {
    int edge_index;
    int u;
    int v;
    double weight;
    uint64_t edge_id;
} TvLocalEdge;

typedef struct {
    const TvRun *run;
    const int *nodes;
    int node_count;
    const TvLocalEdge *edges;
    int edge_count;
    int *state_rows;
    int *best_parents;
    double best_objective;
    bool have_best;
    uint64_t states_explored;
} TvExactSearch;

typedef struct {
    int root;
    uint64_t id;
    size_t slot;
} TvComponentSeed;

static void *checked_alloc(size_t count, size_t width)
{
    void *result;

    if (count == 0) {
        count = 1;
    }
    if (width != 0 && count > SIZE_MAX / width) {
        tv_print_error("allocation overflow");
        exit(EXIT_FAILURE);
    }
    result = calloc(count, width);
    if (result == NULL) {
        tv_print_error("out of memory");
        exit(EXIT_FAILURE);
    }
    return result;
}

static double clamp_value(double value, double low, double high)
{
    if (value < low) {
        return low;
    }
    if (value > high) {
        return high;
    }
    return value;
}

static double logistic(double logit)
{
    if (logit >= 0.0) {
        double term = exp(-logit);
        return 1.0 / (1.0 + term);
    }
    double term = exp(logit);
    return term / (1.0 + term);
}

static double binary_entropy(double score)
{
    score = clamp_value(score, 1e-12, 1.0 - 1e-12);
    return -(score * log(score) + (1.0 - score) * log(1.0 - score))
        / log(2.0);
}

static double normalize_scores(double *values, size_t count)
{
    double maximum = -DBL_MAX;
    double total = 0.0;
    double entropy = 0.0;

    for (size_t index = 0; index < count; index++) {
        if (values[index] > maximum) {
            maximum = values[index];
        }
    }
    for (size_t index = 0; index < count; index++) {
        values[index] = exp(values[index] - maximum);
        total += values[index];
    }
    if (total <= 0.0 || !isfinite(total)) {
        double uniform = 1.0 / (double)count;
        for (size_t index = 0; index < count; index++) {
            values[index] = uniform;
        }
        return 1.0;
    }
    for (size_t index = 0; index < count; index++) {
        values[index] /= total;
        if (values[index] > 0.0) {
            entropy -= values[index] * log(values[index]);
        }
    }
    return count > 1 ? entropy / log((double)count) : 0.0;
}

static int node_key_compare(const TvRun *run, int left, int right)
{
    const TvTE *a = &run->nodes[left];
    const TvTE *b = &run->nodes[right];
    int comparison = strcmp(run->genomes[a->genome].id,
                            run->genomes[b->genome].id);

    if (comparison != 0) {
        return comparison;
    }
    comparison = strcmp(a->id, b->id);
    if (comparison != 0) {
        return comparison;
    }
    comparison = strcmp(a->contig, b->contig);
    if (comparison != 0) {
        return comparison;
    }
    if (a->start != b->start) {
        return a->start < b->start ? -1 : 1;
    }
    return a->end == b->end ? 0 : (a->end < b->end ? -1 : 1);
}

static void sort_nodes(const TvRun *run, int *nodes, size_t count)
{
    for (size_t index = 1; index < count; index++) {
        int value = nodes[index];
        size_t position = index;

        while (position > 0
               && node_key_compare(run, value, nodes[position - 1]) < 0) {
            nodes[position] = nodes[position - 1];
            position--;
        }
        nodes[position] = value;
    }
}

static uint64_t hash_separator(uint64_t hash)
{
    const unsigned char separator = 0;
    return tv_hash_bytes(hash, &separator, 1);
}

static uint64_t hash_u64(uint64_t hash, uint64_t value)
{
    return tv_hash_bytes(hash, &value, sizeof(value));
}

static uint64_t pair_hash(const TvRun *run, const char *kind, int left,
                          int right)
{
    int a = left;
    int b = right;
    uint64_t hash = UINT64_C(1469598103934665603);

    if (node_key_compare(run, a, b) > 0) {
        int swap = a;
        a = b;
        b = swap;
    }
    hash = tv_hash_text(hash, kind);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, run->genomes[run->nodes[a].genome].id);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, run->nodes[a].id);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, run->genomes[run->nodes[b].genome].id);
    hash = hash_separator(hash);
    return tv_hash_text(hash, run->nodes[b].id);
}

static uint64_t node_set_hash(const TvRun *run, const char *kind,
                              const int *nodes, size_t count)
{
    uint64_t hash = UINT64_C(1469598103934665603);

    hash = tv_hash_text(hash, kind);
    for (size_t index = 0; index < count; index++) {
        const TvTE *te = &run->nodes[nodes[index]];
        hash = hash_separator(hash);
        hash = tv_hash_text(hash, run->genomes[te->genome].id);
        hash = hash_separator(hash);
        hash = tv_hash_text(hash, te->id);
    }
    return hash;
}

static bool family_known_text(const char *family)
{
    return family != NULL && family[0] != '\0' && strcmp(family, ".") != 0
        && strcasecmp(family, "unknown") != 0;
}

static void score_projection(const TvRun *run, TvProjection *projection)
{
    const TvPaf *alignment = &run->pafs[projection->paf_index];
    double flank = projection->left_flank_observed
        && projection->right_flank_observed
        ? fmin(projection->left_flank, projection->right_flank) : 0.0;
    double identity = !isnan(projection->identity)
        ? projection->identity
        : (!isnan(alignment->aggregate_identity)
           ? alignment->aggregate_identity : 0.0);
    double mapq = projection->mapq_observed
        ? clamp_value((double)projection->mapq / 60.0, 0.0, 1.0) : 0.0;
    double non_gap = !isnan(projection->n_fraction)
        ? 1.0 - projection->n_fraction : 0.0;
    int core_observed = 0;

    for (size_t index = 0; index < 4; index++) {
        projection->technical_scores[index] = -0.5;
        projection->biological_scores[index] = -0.5;
    }
    for (size_t index = 0; index < 5; index++) {
        projection->annotation_scores[index] = -0.5;
    }

    projection->technical_scores[TV_TECH_CALLABLE] +=
        1.2 * flank + 0.8 * identity + 0.5 * mapq + 0.5 * non_gap;
    projection->technical_scores[TV_TECH_GAP] +=
        isnan(projection->n_fraction) ? 0.2 : 3.0 * projection->n_fraction;
    projection->technical_scores[TV_TECH_AMBIGUOUS] +=
        projection->eligible_candidate_count > 1 ? 1.5 : 0.0;
    projection->technical_scores[TV_TECH_UNCALLABLE] +=
        (!projection->left_flank_observed
         || !projection->right_flank_observed) ? 2.0 : 0.0;
    if (alignment->provider == TV_ALIGNMENT_PAF
        && !projection->mapq_observed) {
        projection->technical_scores[TV_TECH_UNCALLABLE] += 1.5;
    }
    projection->technical_scores[projection->technical_state] += 2.0;

    projection->biological_scores[TV_BIO_PRESENT] +=
        1.5 * projection->te_aligned_fraction
        + (projection->target_te >= 0 ? 1.5 : 0.0);
    projection->biological_scores[TV_BIO_EMPTY] +=
        2.5 * projection->insertion_fraction;
    if (projection->start >= 0 && projection->end >= projection->start) {
        const TvTE *source = &run->nodes[projection->source_te];
        int64_t source_length = source->end - source->start;
        int64_t target_length = projection->end - projection->start;
        if (source_length > 0 && target_length < source_length / 2) {
            projection->biological_scores[TV_BIO_STRUCTURAL_ALTERNATIVE]
                += 1.5;
        }
    }
    projection->biological_scores[TV_BIO_UNKNOWN] +=
        isnan(projection->identity) ? 0.8 : 0.0;
    projection->biological_scores[projection->biological_state] += 2.0;

    if (projection->target_te >= 0) {
        const TvTE *source = &run->nodes[projection->source_te];
        const TvTE *target = &run->nodes[projection->target_te];
        bool source_known = family_known_text(source->family);
        bool target_known = family_known_text(target->family);

        if (source_known && target_known
            && strcasecmp(source->family, target->family) == 0) {
            projection->annotation_scores[TV_ANN_MATCHED] += 1.5;
        } else if (source_known && target_known) {
            projection->annotation_scores[TV_ANN_FAMILY_CONFLICT] += 1.5;
        } else {
            projection->annotation_scores[TV_ANN_UNKNOWN] += 0.8;
        }
    } else {
        projection->annotation_scores[TV_ANN_MISSING] += 1.0;
    }
    projection->annotation_scores[projection->annotation_state] += 2.0;

    projection->technical_entropy = normalize_scores(
        projection->technical_scores, 4);
    projection->biological_entropy = normalize_scores(
        projection->biological_scores, 4);
    projection->annotation_entropy = normalize_scores(
        projection->annotation_scores, 5);

    if (projection->left_flank_observed
        && projection->right_flank_observed) {
        core_observed++;
    }
    if (!isnan(projection->identity)
        || !isnan(alignment->aggregate_identity)) {
        core_observed++;
    }
    if (projection->mapq_observed
        || alignment->provider == TV_ALIGNMENT_MUMMER_DELTA) {
        core_observed++;
    }
    if (!isnan(projection->n_fraction)) {
        core_observed++;
    }
    projection->inference_out_of_domain = core_observed < 3
        || projection->mapping_confidence == TV_MAPPING_AMBIGUOUS;
}

static void score_candidate(const TvRun *run, TvCandidate *candidate)
{
    const TvProjection *projection =
        &run->projections[candidate->projection_index];
    const TvPaf *alignment = &run->pafs[projection->paf_index];
    double logit = -5.5;
    int observed_core = 0;

    candidate->observed_feature_mask =
        TV_FEATURE_RECIPROCAL_OVERLAP | TV_FEATURE_BOUNDARY
        | TV_FEATURE_TE_ALIGNMENT | TV_FEATURE_INSERTION;
    logit += 2.4 * candidate->reciprocal_overlap;
    logit += 1.4 * candidate->boundary_score;
    logit += 0.4 * projection->te_aligned_fraction;

    if (projection->left_flank_observed
        && projection->right_flank_observed) {
        candidate->observed_feature_mask |= TV_FEATURE_FLANK;
        logit += 1.2 * fmin(projection->left_flank,
                            projection->right_flank);
        observed_core++;
    } else {
        logit -= 1.2;
    }
    if (!isnan(projection->identity)) {
        candidate->observed_feature_mask |= TV_FEATURE_LOCAL_IDENTITY;
        logit += projection->identity;
        observed_core++;
    } else if (!isnan(alignment->aggregate_identity)) {
        candidate->observed_feature_mask |= TV_FEATURE_AGGREGATE_IDENTITY;
        logit += 0.7 * alignment->aggregate_identity;
        logit -= 0.15;
        observed_core++;
    } else {
        logit += 0.5;
        logit -= 0.55;
    }
    if (projection->mapq_observed) {
        candidate->observed_feature_mask |= TV_FEATURE_MAPQ;
        logit += 0.5 * clamp_value((double)projection->mapq / 60.0,
                                   0.0, 1.0);
        observed_core++;
    } else {
        logit += 0.25;
        logit -= 0.25;
    }
    if (!isnan(projection->n_fraction)) {
        candidate->observed_feature_mask |= TV_FEATURE_N_FRACTION;
        logit += 0.6 * (1.0 - projection->n_fraction);
        observed_core++;
    } else {
        logit -= 0.3;
    }

    if (candidate->family_relation == TV_FAMILY_MATCH) {
        candidate->observed_feature_mask |= TV_FEATURE_FAMILY;
        logit += 1.0;
    } else if (candidate->family_relation == TV_FAMILY_CONFLICT) {
        candidate->observed_feature_mask |= TV_FEATURE_FAMILY;
        logit -= 2.5;
    } else {
        logit -= 0.15;
    }
    if (candidate->context_relation == TV_CONTEXT_RELATION_SUPPORTED) {
        candidate->observed_feature_mask |= TV_FEATURE_CONTEXT;
        logit += 1.2;
    } else if (candidate->context_relation == TV_CONTEXT_RELATION_CONFLICT) {
        candidate->observed_feature_mask |= TV_FEATURE_CONTEXT;
        logit -= 3.0;
    } else if (candidate->context_relation
               == TV_CONTEXT_RELATION_AMBIGUOUS) {
        candidate->observed_feature_mask |= TV_FEATURE_CONTEXT;
        logit -= 2.0;
    }

    candidate->membership_logit = clamp_value(logit, -20.0, 20.0);
    candidate->membership_score = logistic(candidate->membership_logit);
    candidate->membership_entropy = binary_entropy(
        candidate->membership_score);
    candidate->inference_out_of_domain =
        projection->inference_out_of_domain || observed_core < 3
        || candidate->context_relation == TV_CONTEXT_RELATION_AMBIGUOUS;
}

void tv_score_inference(TvRun *run)
{
    for (size_t index = 0; index < run->n_projections; index++) {
        score_projection(run, &run->projections[index]);
    }
    for (size_t index = 0; index < run->n_candidates; index++) {
        score_candidate(run, &run->candidates[index]);
    }
}

static int parent_root(int *parents, int node)
{
    if (parents[node] != node) {
        parents[node] = parent_root(parents, parents[node]);
    }
    return parents[node];
}

static int parent_root_const(const int *parents, int node)
{
    while (parents[node] != node) {
        node = parents[node];
    }
    return node;
}

static void parent_union(int *parents, int left, int right)
{
    int a = parent_root(parents, left);
    int b = parent_root(parents, right);

    if (a != b) {
        if (a < b) {
            parents[b] = a;
        } else {
            parents[a] = b;
        }
    }
}

static bool assignment_strong(const TvRun *run,
                              const TvTEContext *assignment)
{
    return assignment->assignment == TV_CONTEXT_BRACKETED
        && assignment->context_index >= 0
        && (size_t)assignment->context_index < run->n_contexts
        && strcmp(run->contexts[assignment->context_index].status,
                  "PASS") == 0;
}

static int strong_context_for_group(const TvRun *run, int node,
                                    uint64_t homology_group)
{
    int result = -1;
    size_t begin = run->te_context_offsets == NULL
        ? 0 : run->te_context_offsets[node];
    size_t end = run->te_context_offsets == NULL
        ? 0 : run->te_context_offsets[node + 1];

    for (size_t index = begin; index < end; index++) {
        const TvTEContext *assignment = &run->te_contexts[index];
        if (!assignment_strong(run, assignment)) {
            continue;
        }
        const TvCopyContext *context =
            &run->contexts[assignment->context_index];
        if (homology_group != 0
            && context->homology_group_id != homology_group) {
            continue;
        }
        if (result >= 0 && result != assignment->context_index) {
            return -2;
        }
        result = assignment->context_index;
    }
    return result;
}

static uint64_t matching_group_hash(uint64_t left, uint64_t right)
{
    uint64_t hash = UINT64_C(1469598103934665603);

    if (left > right) {
        uint64_t swap = left;
        left = right;
        right = swap;
    }
    hash = tv_hash_text(hash, "BLOCK_MATCHING");
    hash = hash_separator(hash);
    hash = hash_u64(hash, left);
    hash = hash_separator(hash);
    return hash_u64(hash, right);
}

static int compare_match_item(const void *left, const void *right)
{
    const TvMatchItem *a = left;
    const TvMatchItem *b = right;

    if (a->group_id != b->group_id) {
        return a->group_id < b->group_id ? -1 : 1;
    }
    if (a->edge_id != b->edge_id) {
        return a->edge_id < b->edge_id ? -1 : 1;
    }
    return a->edge_index == b->edge_index
        ? 0 : (a->edge_index < b->edge_index ? -1 : 1);
}

static void finalize_edge_scores(TvRun *run)
{
    for (size_t index = 0; index < run->n_edges; index++) {
        TvEdge *edge = &run->edges[index];
        double logit = edge->membership_logit;

        edge->edge_id = pair_hash(run, "EDGE", edge->a, edge->b);
        if (edge->support_count > 1) {
            logit += 0.15 * log1p((double)(edge->support_count - 1));
        }
        if (edge->independent_reciprocal) {
            logit += 0.35;
        }
        edge->membership_logit = clamp_value(logit, -20.0, 20.0);
        edge->membership_score = logistic(edge->membership_logit);
        edge->membership_entropy = binary_entropy(edge->membership_score);
        edge->matching_selected = true;
        edge->matching_method = TV_MATCH_NOT_APPLICABLE;
        edge->matching_group_id = 0;
        edge->solver_component = -1;
    }
}

static bool append_unique_node(int *nodes, size_t *count, int node)
{
    for (size_t index = 0; index < *count; index++) {
        if (nodes[index] == node) {
            return false;
        }
    }
    nodes[(*count)++] = node;
    return true;
}

static int node_position(const int *nodes, size_t count, int node)
{
    for (size_t index = 0; index < count; index++) {
        if (nodes[index] == node) {
            return (int)index;
        }
    }
    return -1;
}

static void exact_block_matching(TvRun *run, const TvMatchItem *items,
                                 size_t count)
{
    int *left_nodes = checked_alloc(count, sizeof(*left_nodes));
    int *right_nodes = checked_alloc(count, sizeof(*right_nodes));
    size_t left_count = 0;
    size_t right_count = 0;

    for (size_t index = 0; index < count; index++) {
        (void)append_unique_node(left_nodes, &left_count,
                                 items[index].left_node);
        (void)append_unique_node(right_nodes, &right_count,
                                 items[index].right_node);
    }
    sort_nodes(run, left_nodes, left_count);
    sort_nodes(run, right_nodes, right_count);
    size_t dimension = left_count > right_count ? left_count : right_count;
    double *weights = checked_alloc(dimension * dimension,
                                    sizeof(*weights));
    int *cell_edge = checked_alloc(dimension * dimension,
                                   sizeof(*cell_edge));
    double *u = checked_alloc(dimension + 1, sizeof(*u));
    double *v = checked_alloc(dimension + 1, sizeof(*v));
    int *p = checked_alloc(dimension + 1, sizeof(*p));
    int *way = checked_alloc(dimension + 1, sizeof(*way));
    double *minimum = checked_alloc(dimension + 1, sizeof(*minimum));
    bool *used = checked_alloc(dimension + 1, sizeof(*used));

    for (size_t index = 0; index < dimension * dimension; index++) {
        cell_edge[index] = -1;
    }
    for (size_t index = 0; index < count; index++) {
        int row = node_position(left_nodes, left_count,
                                items[index].left_node);
        int column = node_position(right_nodes, right_count,
                                   items[index].right_node);
        double weight = fmax(0.0,
            run->edges[items[index].edge_index].membership_logit);
        size_t cell = (size_t)row * dimension + (size_t)column;

        if (cell_edge[cell] < 0 || weight > weights[cell]
            || (weight == weights[cell]
                && items[index].edge_id
                   < run->edges[cell_edge[cell]].edge_id)) {
            weights[cell] = weight;
            cell_edge[cell] = items[index].edge_index;
        }
        run->edges[items[index].edge_index].matching_selected = false;
        run->edges[items[index].edge_index].matching_method = TV_MATCH_OPTIMAL;
        run->edges[items[index].edge_index].matching_group_id =
            items[index].group_id;
    }

    for (size_t row = 1; row <= dimension; row++) {
        p[0] = (int)row;
        int column0 = 0;
        for (size_t column = 0; column <= dimension; column++) {
            minimum[column] = DBL_MAX;
            used[column] = false;
        }
        do {
            used[column0] = true;
            int row0 = p[column0];
            double delta = DBL_MAX;
            int column1 = 0;
            for (size_t column = 1; column <= dimension; column++) {
                if (used[column]) {
                    continue;
                }
                double cost = -weights[(size_t)(row0 - 1) * dimension
                                       + column - 1];
                double current = cost - u[row0] - v[column];
                if (current < minimum[column]) {
                    minimum[column] = current;
                    way[column] = column0;
                }
                if (minimum[column] < delta) {
                    delta = minimum[column];
                    column1 = (int)column;
                }
            }
            for (size_t column = 0; column <= dimension; column++) {
                if (used[column]) {
                    u[p[column]] += delta;
                    v[column] -= delta;
                } else {
                    minimum[column] -= delta;
                }
            }
            column0 = column1;
        } while (p[column0] != 0);
        do {
            int column1 = way[column0];
            p[column0] = p[column1];
            column0 = column1;
        } while (column0 != 0);
    }
    for (size_t column = 1; column <= dimension; column++) {
        int row = p[column];
        size_t cell = (size_t)(row - 1) * dimension + column - 1;
        int edge_index = cell_edge[cell];

        if (edge_index >= 0 && weights[cell] > 0.0) {
            run->edges[edge_index].matching_selected = true;
        }
    }
    free(used);
    free(minimum);
    free(way);
    free(p);
    free(v);
    free(u);
    free(cell_edge);
    free(weights);
    free(right_nodes);
    free(left_nodes);
}

static int compare_weighted_edge(const void *left, const void *right)
{
    const TvLocalEdge *a = left;
    const TvLocalEdge *b = right;

    if (a->weight != b->weight) {
        return a->weight > b->weight ? -1 : 1;
    }
    if (a->edge_id != b->edge_id) {
        return a->edge_id < b->edge_id ? -1 : 1;
    }
    return a->edge_index == b->edge_index
        ? 0 : (a->edge_index < b->edge_index ? -1 : 1);
}

static void heuristic_block_matching(TvRun *run, const TvMatchItem *items,
                                     size_t count)
{
    TvLocalEdge *order = checked_alloc(count, sizeof(*order));
    int *used_left = checked_alloc(count, sizeof(*used_left));
    int *used_right = checked_alloc(count, sizeof(*used_right));
    size_t left_count = 0;
    size_t right_count = 0;

    for (size_t index = 0; index < count; index++) {
        TvEdge *edge = &run->edges[items[index].edge_index];
        order[index] = (TvLocalEdge){
            .edge_index = items[index].edge_index,
            .u = items[index].left_node,
            .v = items[index].right_node,
            .weight = edge->membership_logit,
            .edge_id = edge->edge_id
        };
        edge->matching_selected = false;
        edge->matching_method = TV_MATCH_HEURISTIC;
        edge->matching_group_id = items[index].group_id;
    }
    qsort(order, count, sizeof(*order), compare_weighted_edge);
    for (size_t index = 0; index < count; index++) {
        bool left_seen = false;
        bool right_seen = false;
        for (size_t used = 0; used < left_count; used++) {
            left_seen = left_seen || used_left[used] == order[index].u;
        }
        for (size_t used = 0; used < right_count; used++) {
            right_seen = right_seen || used_right[used] == order[index].v;
        }
        if (!left_seen && !right_seen && order[index].weight > 0.0) {
            run->edges[order[index].edge_index].matching_selected = true;
            used_left[left_count++] = order[index].u;
            used_right[right_count++] = order[index].v;
        }
    }
    free(used_right);
    free(used_left);
    free(order);
}

static bool edge_passes_pre_matching_gates(const TvRun *run,
                                           const TvEdge *edge)
{
    return edge->score >= run->cfg.min_edge_score
        && edge->family_compatible
        && edge->membership_score >= run->cfg.min_membership_score
        && edge->membership_logit > 0.0;
}

static void match_edges(TvRun *run)
{
    TvMatchItem *items = checked_alloc(run->n_edges, sizeof(*items));
    size_t count = 0;

    for (size_t index = 0; index < run->n_edges; index++) {
        TvEdge *edge = &run->edges[index];
        uint64_t group = edge->shared_homology_group_id;
        TvContextRelation relation = edge->context_relation;

        if (relation != TV_CONTEXT_RELATION_SUPPORTED || group == 0
            || !edge_passes_pre_matching_gates(run, edge)) {
            continue;
        }
        int context_a = strong_context_for_group(run, edge->a, group);
        int context_b = strong_context_for_group(run, edge->b, group);
        if (context_a < 0 || context_b < 0) {
            edge->matching_selected = false;
            edge->matching_method = TV_MATCH_NOT_APPLICABLE;
            continue;
        }
        const TvCopyContext *a = &run->contexts[context_a];
        const TvCopyContext *b = &run->contexts[context_b];
        bool forward = a->context_id < b->context_id;
        uint64_t match_group = matching_group_hash(a->context_id,
                                                   b->context_id);
        items[count++] = (TvMatchItem){
            .edge_index = (int)index,
            .left_node = forward ? edge->a : edge->b,
            .right_node = forward ? edge->b : edge->a,
            .left_context = forward ? context_a : context_b,
            .right_context = forward ? context_b : context_a,
            .group_id = match_group,
            .edge_id = edge->edge_id
        };
    }
    qsort(items, count, sizeof(*items), compare_match_item);
    size_t start = 0;
    while (start < count) {
        size_t end = start + 1;
        size_t left_count = 0;
        size_t right_count = 0;
        int *left = checked_alloc(count - start, sizeof(*left));
        int *right = checked_alloc(count - start, sizeof(*right));

        while (end < count && items[end].group_id == items[start].group_id) {
            end++;
        }
        for (size_t index = start; index < end; index++) {
            (void)append_unique_node(left, &left_count,
                                     items[index].left_node);
            (void)append_unique_node(right, &right_count,
                                     items[index].right_node);
        }
        if (left_count <= (size_t)run->cfg.exact_matching_max_nodes
            && right_count <= (size_t)run->cfg.exact_matching_max_nodes) {
            exact_block_matching(run, &items[start], end - start);
        } else {
            heuristic_block_matching(run, &items[start], end - start);
        }
        free(right);
        free(left);
        start = end;
    }
    free(items);
}

static bool node_has_strong_context(const TvRun *run, int node)
{
    size_t begin = run->te_context_offsets == NULL
        ? 0 : run->te_context_offsets[node];
    size_t end = run->te_context_offsets == NULL
        ? 0 : run->te_context_offsets[node + 1];

    for (size_t index = begin; index < end; index++) {
        if (assignment_strong(run, &run->te_contexts[index])) {
            return true;
        }
    }
    return false;
}

static bool local_merge_allowed(const TvRun *run, const int *nodes,
                                int node_count, const int *parents,
                                int left_root, int right_root)
{
    for (int left = 0; left < node_count; left++) {
        if (parent_root_const(parents, left) != left_root) {
            continue;
        }
        const TvTE *a = &run->nodes[nodes[left]];
        bool a_known = family_known_text(a->family);
        for (int right = 0; right < node_count; right++) {
            if (parent_root_const(parents, right) != right_root) {
                continue;
            }
            const TvTE *b = &run->nodes[nodes[right]];
            bool b_known = family_known_text(b->family);
            if (a_known && b_known
                && strcasecmp(a->family, b->family) != 0) {
                return false;
            }
        }
    }

    for (int left = 0; left < node_count; left++) {
        if (parent_root_const(parents, left) != left_root) {
            continue;
        }
        size_t left_begin = run->te_context_offsets == NULL
            ? 0 : run->te_context_offsets[nodes[left]];
        size_t left_end = run->te_context_offsets == NULL
            ? 0 : run->te_context_offsets[nodes[left] + 1];
        for (size_t a_index = left_begin; a_index < left_end; a_index++) {
            const TvTEContext *a_assignment = &run->te_contexts[a_index];
            if (!assignment_strong(run, a_assignment)) {
                continue;
            }
            const TvCopyContext *a =
                &run->contexts[a_assignment->context_index];
            for (int right = 0; right < node_count; right++) {
                if (parent_root_const(parents, right) != right_root) {
                    continue;
                }
                size_t right_begin = run->te_context_offsets == NULL
                    ? 0 : run->te_context_offsets[nodes[right]];
                size_t right_end = run->te_context_offsets == NULL
                    ? 0 : run->te_context_offsets[nodes[right] + 1];
                for (size_t b_index = right_begin; b_index < right_end;
                     b_index++) {
                    const TvTEContext *b_assignment =
                        &run->te_contexts[b_index];
                    if (!assignment_strong(run, b_assignment)) {
                        continue;
                    }
                    const TvCopyContext *b =
                        &run->contexts[b_assignment->context_index];
                    if (a_assignment->context_index
                        == b_assignment->context_index
                        || a->homology_group_id != b->homology_group_id) {
                        return false;
                    }
                }
            }
        }
    }

    for (size_t genome = 0; genome < run->n_genomes; genome++) {
        int copies = 0;
        for (int index = 0; index < node_count; index++) {
            int root = parent_root_const(parents, index);
            if ((root == left_root || root == right_root)
                && run->nodes[nodes[index]].genome == (int)genome
                && !node_has_strong_context(run, nodes[index])) {
                copies++;
            }
        }
        if (copies > run->genomes[genome].max_locus_copies) {
            return false;
        }
    }
    return true;
}

static double partition_objective(const TvLocalEdge *edges, int edge_count,
                                  const int *parents)
{
    double objective = 0.0;

    for (int index = 0; index < edge_count; index++) {
        if (parent_root_const(parents, edges[index].u)
            == parent_root_const(parents, edges[index].v)) {
            objective += edges[index].weight;
        }
    }
    return objective;
}

static void exact_search(TvExactSearch *search, int depth,
                         const int *parents)
{
    if (depth == search->edge_count) {
        double objective = partition_objective(
            search->edges, search->edge_count, parents);
        search->states_explored++;
        if (!search->have_best
            || objective > search->best_objective + 1e-12) {
            memcpy(search->best_parents, parents,
                   (size_t)search->node_count * sizeof(*parents));
            search->best_objective = objective;
            search->have_best = true;
        }
        return;
    }

    const TvLocalEdge *edge = &search->edges[depth];
    int left_root = parent_root_const(parents, edge->u);
    int right_root = parent_root_const(parents, edge->v);

    if (left_root == right_root) {
        exact_search(search, depth + 1, parents);
        return;
    }
    if (local_merge_allowed(search->run, search->nodes,
                            search->node_count, parents,
                            left_root, right_root)) {
        int *next = search->state_rows
            + (size_t)(depth + 1) * (size_t)search->node_count;
        memcpy(next, parents,
               (size_t)search->node_count * sizeof(*parents));
        if (left_root < right_root) {
            next[right_root] = left_root;
        } else {
            next[left_root] = right_root;
        }
        exact_search(search, depth + 1, next);
    }
    exact_search(search, depth + 1, parents);
}

static void solve_exact_component(const TvRun *run, const int *nodes,
                                  int node_count, const TvLocalEdge *edges,
                                  int edge_count, int *parents,
                                  TvSolverComponent *result)
{
    int *state_rows = checked_alloc((size_t)(edge_count + 1)
                                    * (size_t)node_count,
                                    sizeof(*state_rows));
    int *best = checked_alloc((size_t)node_count, sizeof(*best));
    TvExactSearch search = {
        .run = run,
        .nodes = nodes,
        .node_count = node_count,
        .edges = edges,
        .edge_count = edge_count,
        .state_rows = state_rows,
        .best_parents = best,
        .best_objective = -DBL_MAX
    };

    for (int index = 0; index < node_count; index++) {
        state_rows[index] = index;
    }
    exact_search(&search, 0, state_rows);
    memcpy(parents, best, (size_t)node_count * sizeof(*parents));
    result->method = edge_count == 0
        ? TV_SOLVER_TRIVIAL : TV_SOLVER_EXACT_ENUMERATION;
    result->status = TV_SOLVER_OPTIMAL;
    result->objective = search.have_best ? search.best_objective : 0.0;
    result->upper_bound = result->objective;
    result->relative_gap = 0.0;
    result->states_explored = search.states_explored;
    free(best);
    free(state_rows);
}

static void solve_greedy_component(const TvRun *run, const int *nodes,
                                   int node_count, const TvLocalEdge *edges,
                                   int edge_count, int *parents,
                                   TvSolverComponent *result)
{
    double upper_bound = 0.0;

    for (int index = 0; index < node_count; index++) {
        parents[index] = index;
    }
    for (int index = 0; index < edge_count; index++) {
        int left_root = parent_root(parents, edges[index].u);
        int right_root = parent_root(parents, edges[index].v);
        upper_bound += fmax(0.0, edges[index].weight);
        if (left_root != right_root && edges[index].weight > 0.0
            && local_merge_allowed(run, nodes, node_count, parents,
                                   left_root, right_root)) {
            parent_union(parents, left_root, right_root);
        }
    }
    result->method = TV_SOLVER_DETERMINISTIC_GREEDY;
    result->status = TV_SOLVER_HEURISTIC;
    result->objective = partition_objective(edges, edge_count, parents);
    result->upper_bound = upper_bound;
    result->relative_gap = upper_bound <= result->objective + 1e-12
        ? 0.0
        : (upper_bound - result->objective)
          / fmax(fabs(result->objective), 1e-12);
    result->states_explored = (uint64_t)edge_count;
}

static bool edge_active(const TvRun *run, TvEdge *edge)
{
    if (edge->score < run->cfg.min_edge_score) {
        (void)snprintf(edge->selection_reason,
                       sizeof(edge->selection_reason),
                       "BELOW_EDGE_THRESHOLD");
        return false;
    }
    if (!edge->family_compatible) {
        (void)snprintf(edge->selection_reason,
                       sizeof(edge->selection_reason),
                       "DIRECT_FAMILY_CONFLICT");
        return false;
    }
    if (edge->context_relation == TV_CONTEXT_RELATION_CONFLICT) {
        (void)snprintf(edge->selection_reason,
                       sizeof(edge->selection_reason),
                       "BASE_SYNTENY_CONFLICT");
        return false;
    }
    if (edge->context_relation == TV_CONTEXT_RELATION_AMBIGUOUS) {
        (void)snprintf(edge->selection_reason,
                       sizeof(edge->selection_reason),
                       "SYNTENY_CONTEXT_AMBIGUOUS");
        return false;
    }
    if (edge->membership_score < run->cfg.min_membership_score
        || edge->membership_logit <= 0.0) {
        (void)snprintf(edge->selection_reason,
                       sizeof(edge->selection_reason),
                       "BELOW_MEMBERSHIP_THRESHOLD");
        return false;
    }
    if (!edge->matching_selected) {
        (void)snprintf(edge->selection_reason,
                       sizeof(edge->selection_reason),
                       "BLOCK_MATCHING_CONFLICT");
        return false;
    }
    (void)snprintf(edge->selection_reason,
                   sizeof(edge->selection_reason), "SOLVER_PENDING");
    return true;
}

static int compare_component_seed(const void *left, const void *right)
{
    const TvComponentSeed *a = left;
    const TvComponentSeed *b = right;

    if (a->id != b->id) {
        return a->id < b->id ? -1 : 1;
    }
    return a->root == b->root ? 0 : (a->root < b->root ? -1 : 1);
}

static void append_solver(TvRun *run, const TvSolverComponent *component)
{
    if (run->n_solver_components == run->cap_solver_components) {
        run->cap_solver_components = run->cap_solver_components == 0
            ? 32 : run->cap_solver_components * 2;
        run->solver_components = tv_grow(
            run->solver_components, run->cap_solver_components,
            sizeof(*run->solver_components));
    }
    run->solver_components[run->n_solver_components++] = *component;
}

static size_t group_component_nodes(const TvRun *run, int *parents,
                                    const char *kind,
                                    TvComponentSeed **seeds_out,
                                    size_t **offsets_out, int **members_out)
{
    TvComponentSeed *seeds = checked_alloc(run->n_nodes, sizeof(*seeds));
    int *root_to_slot = checked_alloc(run->n_nodes, sizeof(*root_to_slot));
    size_t *counts = checked_alloc(run->n_nodes, sizeof(*counts));
    size_t *offsets = checked_alloc(run->n_nodes + 1, sizeof(*offsets));
    size_t *cursor = checked_alloc(run->n_nodes, sizeof(*cursor));
    int *members = checked_alloc(run->n_nodes, sizeof(*members));
    size_t seed_count = 0;

    for (size_t index = 0; index < run->n_nodes; index++) {
        root_to_slot[index] = -1;
    }
    for (size_t node = 0; node < run->n_nodes; node++) {
        int root = parent_root(parents, (int)node);
        int slot = root_to_slot[root];
        if (slot < 0) {
            slot = (int)seed_count;
            root_to_slot[root] = slot;
            seeds[seed_count] = (TvComponentSeed){
                .root = root,
                .slot = seed_count
            };
            seed_count++;
        }
        counts[slot]++;
    }
    for (size_t slot = 0; slot < seed_count; slot++) {
        offsets[slot + 1] = offsets[slot] + counts[slot];
        cursor[slot] = offsets[slot];
    }
    for (size_t node = 0; node < run->n_nodes; node++) {
        int root = parent_root(parents, (int)node);
        int slot = root_to_slot[root];
        members[cursor[slot]++] = (int)node;
    }
    for (size_t slot = 0; slot < seed_count; slot++) {
        size_t begin = offsets[slot];
        size_t count = offsets[slot + 1] - begin;
        sort_nodes(run, &members[begin], count);
        seeds[slot].id = node_set_hash(run, kind, &members[begin], count);
    }
    qsort(seeds, seed_count, sizeof(*seeds), compare_component_seed);
    free(cursor);
    free(counts);
    free(root_to_slot);
    *seeds_out = seeds;
    *offsets_out = offsets;
    *members_out = members;
    return seed_count;
}

static int assign_locus_components(TvRun *run, int *parents)
{
    TvComponentSeed *seeds = NULL;
    size_t *offsets = NULL;
    int *members = NULL;
    size_t count = group_component_nodes(
        run, parents, "LOCUS", &seeds, &offsets, &members);
    int *root_to_locus = checked_alloc(run->n_nodes,
                                        sizeof(*root_to_locus));

    for (size_t index = 0; index < run->n_nodes; index++) {
        root_to_locus[index] = -1;
    }
    for (size_t index = 0; index < count; index++) {
        root_to_locus[seeds[index].root] = (int)index;
    }
    run->components = tv_grow(NULL, run->n_nodes, sizeof(*run->components));
    for (size_t node = 0; node < run->n_nodes; node++) {
        int root = parent_root(parents, (int)node);
        run->components[node] = root_to_locus[root];
    }
    run->n_loci = (int)count;
    free(root_to_locus);
    free(members);
    free(offsets);
    free(seeds);
    return 0;
}

static int solve_components(TvRun *run)
{
    int *connectivity = checked_alloc(run->n_nodes, sizeof(*connectivity));
    int *final_parents = checked_alloc(run->n_nodes, sizeof(*final_parents));
    bool *active = checked_alloc(run->n_edges, sizeof(*active));
    TvComponentSeed *seeds = NULL;
    size_t *node_offsets = NULL;
    int *component_nodes = NULL;

    for (size_t node = 0; node < run->n_nodes; node++) {
        connectivity[node] = (int)node;
        final_parents[node] = (int)node;
    }
    for (size_t index = 0; index < run->n_edges; index++) {
        active[index] = edge_active(run, &run->edges[index]);
        if (active[index]) {
            parent_union(connectivity, run->edges[index].a,
                         run->edges[index].b);
        }
    }
    size_t seed_count = group_component_nodes(
        run, connectivity, "SOLVER_COMPONENT", &seeds, &node_offsets,
        &component_nodes);
    int *root_to_seed = checked_alloc(run->n_nodes, sizeof(*root_to_seed));
    size_t *edge_offsets = checked_alloc(seed_count + 1,
                                          sizeof(*edge_offsets));
    size_t *edge_cursor = checked_alloc(seed_count, sizeof(*edge_cursor));
    int *component_edges = checked_alloc(run->n_edges,
                                          sizeof(*component_edges));
    int *global_to_local = checked_alloc(run->n_nodes,
                                          sizeof(*global_to_local));

    for (size_t index = 0; index < run->n_nodes; index++) {
        root_to_seed[index] = -1;
        global_to_local[index] = -1;
    }
    for (size_t index = 0; index < seed_count; index++) {
        root_to_seed[seeds[index].root] = (int)index;
    }
    for (size_t index = 0; index < run->n_edges; index++) {
        if (active[index]) {
            int root = parent_root(connectivity, run->edges[index].a);
            int seed = root_to_seed[root];
            edge_offsets[(size_t)seed + 1]++;
        }
    }
    for (size_t index = 0; index < seed_count; index++) {
        edge_offsets[index + 1] += edge_offsets[index];
        edge_cursor[index] = edge_offsets[index];
    }
    for (size_t index = 0; index < run->n_edges; index++) {
        if (active[index]) {
            int root = parent_root(connectivity, run->edges[index].a);
            int seed = root_to_seed[root];
            component_edges[edge_cursor[seed]++] = (int)index;
        }
    }

    for (size_t seed_index = 0; seed_index < seed_count; seed_index++) {
        size_t slot = seeds[seed_index].slot;
        size_t node_begin = node_offsets[slot];
        int node_count = (int)(node_offsets[slot + 1] - node_begin);
        int *nodes = &component_nodes[node_begin];
        size_t edge_begin = edge_offsets[seed_index];
        int edge_count = (int)(edge_offsets[seed_index + 1] - edge_begin);

        for (int index = 0; index < node_count; index++) {
            global_to_local[nodes[index]] = index;
        }
        TvLocalEdge *edges = checked_alloc((size_t)edge_count,
                                           sizeof(*edges));
        for (int index = 0; index < edge_count; index++) {
            int global_edge = component_edges[edge_begin + (size_t)index];
            edges[index] = (TvLocalEdge){
                .edge_index = global_edge,
                .u = global_to_local[run->edges[global_edge].a],
                .v = global_to_local[run->edges[global_edge].b],
                .weight = run->edges[global_edge].membership_logit,
                .edge_id = run->edges[global_edge].edge_id
            };
        }
        qsort(edges, (size_t)edge_count, sizeof(*edges),
              compare_weighted_edge);
        int *local_parents = checked_alloc((size_t)node_count,
                                           sizeof(*local_parents));
        TvSolverComponent result = {
            .solver_id = seeds[seed_index].id,
            .node_count = node_count,
            .edge_count = edge_count
        };

        if (edge_count <= run->cfg.exact_max_edges) {
            solve_exact_component(run, nodes, node_count, edges, edge_count,
                                  local_parents, &result);
        } else {
            solve_greedy_component(run, nodes, node_count, edges, edge_count,
                                   local_parents, &result);
        }
        int solver_index = (int)run->n_solver_components;
        append_solver(run, &result);
        int *local_representative = checked_alloc(
            (size_t)node_count, sizeof(*local_representative));
        for (int index = 0; index < node_count; index++) {
            local_representative[index] = -1;
        }
        for (int index = 0; index < node_count; index++) {
            int root = parent_root_const(local_parents, index);
            if (local_representative[root] < 0) {
                local_representative[root] = nodes[index];
            } else {
                parent_union(final_parents, local_representative[root],
                             nodes[index]);
            }
        }
        for (int index = 0; index < edge_count; index++) {
            TvEdge *edge = &run->edges[edges[index].edge_index];
            edge->solver_component = solver_index;
            edge->selected = parent_root_const(local_parents, edges[index].u)
                == parent_root_const(local_parents, edges[index].v);
            (void)snprintf(
                edge->selection_reason, sizeof(edge->selection_reason), "%s",
                edge->selected
                    ? (result.status == TV_SOLVER_OPTIMAL
                       ? "SELECTED_EXACT" : "SELECTED_HEURISTIC")
                    : "GLOBAL_CONSTRAINT_SEPARATED");
        }
        for (int index = 0; index < node_count; index++) {
            global_to_local[nodes[index]] = -1;
        }
        free(local_representative);
        free(local_parents);
        free(edges);
    }
    int status = assign_locus_components(run, final_parents);
    free(global_to_local);
    free(component_edges);
    free(edge_cursor);
    free(edge_offsets);
    free(root_to_seed);
    free(component_nodes);
    free(node_offsets);
    free(seeds);
    free(active);
    free(final_parents);
    free(connectivity);
    return status;
}

static bool text_known(const char *text)
{
    return text != NULL && text[0] != '\0' && strcmp(text, ".") != 0;
}

static void relation_normalize(TvRelation *relation)
{
    double total = 0.0;
    double entropy = 0.0;

    for (int index = 0; index < TV_RELATION_COUNT; index++) {
        relation->scores[index] = fmax(0.0, relation->scores[index]);
        total += relation->scores[index];
    }
    if (total == 0.0) {
        relation->scores[TV_RELATION_UNKNOWN] = 1.0;
        total = 1.0;
    }
    relation->predicted = TV_RELATION_UNKNOWN;
    for (int index = 0; index < TV_RELATION_COUNT; index++) {
        relation->scores[index] /= total;
        if (relation->scores[index] > 0.0) {
            entropy -= relation->scores[index]
                * log(relation->scores[index]);
        }
        if (relation->scores[index]
            > relation->scores[relation->predicted]) {
            relation->predicted = (TvRelationClass)index;
        }
    }
    relation->entropy = entropy / log((double)TV_RELATION_COUNT);
}

static void classify_relation(const TvRun *run, TvRelation *relation)
{
    const TvTE *a = &run->nodes[relation->a];
    const TvTE *b = &run->nodes[relation->b];
    const TvEdge *edge = relation->edge_index >= 0
        ? &run->edges[relation->edge_index] : NULL;
    bool same_locus = relation->locus >= 0;
    int context_a = strong_context_for_group(run, relation->a, 0);
    int context_b = strong_context_for_group(run, relation->b, 0);

    for (int index = 0; index < TV_RELATION_COUNT; index++) {
        relation->scores[index] = 0.01;
    }
    relation->scores[TV_RELATION_UNKNOWN] = 0.20;

    if (same_locus && a->genome != b->genome) {
        relation->scores[TV_RELATION_ORTHOLOG] = edge == NULL
            ? 0.75 : 0.70 + 0.25 * edge->membership_score;
        relation->out_of_domain = edge != NULL
            && edge->inference_out_of_domain;
    } else if (same_locus && a->genome == b->genome
               && context_a >= 0 && context_b >= 0
               && context_a != context_b
               && run->contexts[context_a].homology_group_id
                  == run->contexts[context_b].homology_group_id) {
        const TvCopyContext *left = &run->contexts[context_a];
        const TvCopyContext *right = &run->contexts[context_b];
        bool allelic = text_known(left->haplotype_id)
            && text_known(right->haplotype_id)
            && strcmp(left->haplotype_id, right->haplotype_id) != 0
            && (!text_known(left->subgenome_id)
                || !text_known(right->subgenome_id)
                || strcmp(left->subgenome_id, right->subgenome_id) == 0);
        bool wgd = text_known(left->wgd_node)
            && text_known(right->wgd_node)
            && strcmp(left->wgd_node, right->wgd_node) == 0;

        if (allelic) {
            relation->scores[TV_RELATION_ALLELIC] = 0.80;
        } else if (wgd) {
            relation->scores[TV_RELATION_WGD_HOMEOLOG] = 0.80;
        } else {
            relation->scores[TV_RELATION_WGD_HOMEOLOG] = 0.40;
            relation->scores[TV_RELATION_UNKNOWN] = 0.50;
            relation->out_of_domain = true;
        }
    } else if (!same_locus && a->genome == b->genome
               && strcmp(a->contig, b->contig) == 0) {
        int64_t distance = a->end <= b->start
            ? b->start - a->end : (b->end <= a->start
                ? a->start - b->end : 0);
        if (distance <= run->cfg.tandem_distance) {
            relation->scores[TV_RELATION_TANDEM_PARALOG] = 0.70;
            relation->scores[TV_RELATION_UNKNOWN] = 0.22;
        }
    } else if (edge != NULL && !edge->selected) {
        if (edge->context_relation == TV_CONTEXT_RELATION_CONFLICT
            && edge->family_compatible && edge->membership_score >= 0.90) {
            relation->scores[TV_RELATION_TRANSPOSED_PARALOG] = 0.55;
            relation->scores[TV_RELATION_UNKNOWN] = 0.35;
            relation->out_of_domain = true;
        }
    }
    relation_normalize(relation);
}

static int edge_index_for_pair(const TvRun *run, int left, int right)
{
    int a = left < right ? left : right;
    int b = left < right ? right : left;
    size_t low = 0;
    size_t high = run->n_edges;

    while (low < high) {
        size_t middle = low + (high - low) / 2;
        const TvEdge *edge = &run->edges[middle];
        if (edge->a < a || (edge->a == a && edge->b < b)) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    return low < run->n_edges && run->edges[low].a == a
        && run->edges[low].b == b ? (int)low : -1;
}

static void append_relation(TvRun *run, int left, int right, int edge_index,
                            int locus)
{
    if (run->n_relations == run->cap_relations) {
        run->cap_relations = run->cap_relations == 0
            ? 64 : run->cap_relations * 2;
        run->relations = tv_grow(run->relations, run->cap_relations,
                                 sizeof(*run->relations));
    }
    TvRelation relation;
    memset(&relation, 0, sizeof(relation));
    relation.relation_id = pair_hash(run, "RELATION", left, right);
    relation.a = left;
    relation.b = right;
    relation.edge_index = edge_index;
    relation.locus = locus;
    classify_relation(run, &relation);
    run->relations[run->n_relations++] = relation;
}

static int compare_relation_id(const void *left, const void *right)
{
    const TvRelation *a = left;
    const TvRelation *b = right;
    return a->relation_id == b->relation_id
        ? 0 : (a->relation_id < b->relation_id ? -1 : 1);
}

static void build_relations(TvRun *run)
{
    for (size_t index = 0; index < run->n_edges; index++) {
        TvEdge *edge = &run->edges[index];
        int locus = edge->selected ? run->components[edge->a] : -1;
        append_relation(run, edge->a, edge->b, (int)index, locus);
    }
    size_t *locus_offsets = checked_alloc((size_t)run->n_loci + 1,
                                           sizeof(*locus_offsets));
    size_t *locus_cursor = checked_alloc((size_t)run->n_loci,
                                          sizeof(*locus_cursor));
    int *locus_members = checked_alloc(run->n_nodes,
                                        sizeof(*locus_members));
    for (size_t node = 0; node < run->n_nodes; node++) {
        locus_offsets[(size_t)run->components[node] + 1]++;
    }
    for (int locus = 0; locus < run->n_loci; locus++) {
        locus_offsets[(size_t)locus + 1] += locus_offsets[locus];
        locus_cursor[locus] = locus_offsets[locus];
    }
    for (size_t node = 0; node < run->n_nodes; node++) {
        int locus = run->components[node];
        locus_members[locus_cursor[locus]++] = (int)node;
    }
    for (int locus = 0; locus < run->n_loci; locus++) {
        size_t begin = locus_offsets[locus];
        size_t end = locus_offsets[(size_t)locus + 1];
        for (size_t left = begin; left < end; left++) {
            for (size_t right = left + 1; right < end; right++) {
                if (edge_index_for_pair(run, locus_members[left],
                                        locus_members[right]) < 0) {
                    append_relation(run, locus_members[left],
                                    locus_members[right], -1, locus);
                }
            }
        }
    }
    for (size_t genome = 0; genome < run->n_genomes; genome++) {
        size_t begin = run->genomes[genome].node_offset;
        size_t end = begin + run->genomes[genome].n_tes;
        for (size_t left = begin; left < end; left++) {
            const TvTE *a = &run->nodes[left];
            for (size_t right = left + 1; right < end; right++) {
                const TvTE *b = &run->nodes[right];
                if (strcmp(a->contig, b->contig) != 0) {
                    break;
                }
                int64_t distance = a->end <= b->start
                    ? b->start - a->end : 0;
                if (distance > run->cfg.tandem_distance) {
                    break;
                }
                bool family_compatible = !family_known_text(a->family)
                    || !family_known_text(b->family)
                    || strcasecmp(a->family, b->family) == 0;
                if (run->components[left] != run->components[right]
                    && family_compatible
                    && edge_index_for_pair(run, (int)left, (int)right) < 0) {
                    append_relation(run, (int)left, (int)right, -1, -1);
                }
            }
        }
    }
    free(locus_members);
    free(locus_cursor);
    free(locus_offsets);
    if (run->n_relations > 1) {
        qsort(run->relations, run->n_relations,
              sizeof(*run->relations), compare_relation_id);
    }
}

int tv_infer_loci(TvRun *run)
{
    if (run == NULL || run->components != NULL
        || run->n_solver_components != 0 || run->n_relations != 0) {
        tv_print_error("v0.5 inference requires one fresh graph");
        return -1;
    }
    finalize_edge_scores(run);
    match_edges(run);
    if (solve_components(run) != 0) {
        return -1;
    }
    build_relations(run);
    return 0;
}
