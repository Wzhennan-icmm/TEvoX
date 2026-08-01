#include "tevox.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

static void *grow(void *pointer, size_t count, size_t width)
{
    void *result;

    if (count == 0) {
        count = 1;
    }
    if (width != 0 && count > SIZE_MAX / width) {
        tv_print_error("allocation overflow");
        exit(EXIT_FAILURE);
    }
    result = realloc(pointer, count * width);
    if (result == NULL) {
        tv_print_error("out of memory");
        exit(EXIT_FAILURE);
    }
    return result;
}

static char *dupstr(const char *text)
{
    char *copy = strdup(text == NULL ? "" : text);

    if (copy == NULL) {
        tv_print_error("out of memory");
        exit(EXIT_FAILURE);
    }
    return copy;
}

static double clamp(double value, double low, double high)
{
    if (value < low) {
        return low;
    }
    if (value > high) {
        return high;
    }
    return value;
}

static int64_t minimum_i64(int64_t left, int64_t right)
{
    return left < right ? left : right;
}

static int64_t maximum_i64(int64_t left, int64_t right)
{
    return left > right ? left : right;
}

static int64_t overlap(int64_t a_start, int64_t a_end,
                       int64_t b_start, int64_t b_end)
{
    int64_t start = maximum_i64(a_start, b_start);
    int64_t end = minimum_i64(a_end, b_end);

    return end > start ? end - start : 0;
}

static uint64_t hash_bytes(uint64_t hash, const void *data, size_t length)
{
    const unsigned char *bytes = data;

    for (size_t index = 0; index < length; index++) {
        hash ^= (uint64_t)bytes[index];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

static uint64_t hash_text(uint64_t hash, const char *text)
{
    return hash_bytes(hash, text, strlen(text));
}

static uint64_t projection_id(const TvRun *run, int source_te,
                              int target_genome, uint64_t group_id)
{
    const TvTE *source = &run->nodes[source_te];
    uint64_t hash = UINT64_C(1469598103934665603);

    hash = hash_bytes(hash, &group_id, sizeof(group_id));
    hash = hash_text(hash, run->genomes[source->genome].id);
    hash = hash_bytes(hash, "\0", 1);
    hash = hash_text(hash, source->id);
    hash = hash_bytes(hash, "\0", 1);
    return hash_text(hash, run->genomes[target_genome].id);
}

static uint64_t decision_id(const TvRun *run, int source_te, int target_genome)
{
    const TvTE *source = &run->nodes[source_te];
    uint64_t hash = UINT64_C(1469598103934665603);

    hash = hash_text(hash, run->genomes[source->genome].id);
    hash = hash_bytes(hash, "\0", 1);
    hash = hash_text(hash, source->id);
    hash = hash_bytes(hash, "\0", 1);
    return hash_text(hash, run->genomes[target_genome].id);
}

static uint64_t candidate_id(const TvRun *run, uint64_t evidence_id,
                             int target_te)
{
    const TvTE *target = &run->nodes[target_te];
    uint64_t hash = UINT64_C(1469598103934665603);

    hash = hash_bytes(hash, &evidence_id, sizeof(evidence_id));
    hash = hash_text(hash, run->genomes[target->genome].id);
    hash = hash_bytes(hash, "\0", 1);
    return hash_text(hash, target->id);
}

static TvContig *contig_at(TvGenome *genome, int index)
{
    if (index < 0 || index >= (int)genome->n_contigs) {
        return NULL;
    }
    return &genome->contigs[index];
}

static double n_fraction(TvRun *run, int genome_index, int contig_index,
                         int64_t start, int64_t end)
{
    TvGenome *genome = &run->genomes[genome_index];
    TvContig *contig = contig_at(genome, contig_index);
    int64_t ambiguous = 0;
    size_t low = 0;
    size_t high;

    if (contig == NULL || contig->line_bases == 0 || start < 0
        || end > contig->length || end <= start) {
        return 1.0;
    }
    run->performance.gap_queries++;
    high = contig->n_gap_runs;
    while (low < high) {
        size_t middle = low + (high - low) / 2;

        if (contig->gap_runs[middle].end <= start) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    for (size_t index = low; index < contig->n_gap_runs; index++) {
        TvGapRun run_interval = contig->gap_runs[index];

        if (run_interval.start >= end) {
            break;
        }
        ambiguous += overlap(run_interval.start, run_interval.end, start, end);
    }
    return (double)ambiguous / (double)(end - start);
}

static bool consumes_query(char code)
{
    return strchr("M=XI", code) != NULL;
}

static bool consumes_target(char code)
{
    return strchr("M=XD", code) != NULL;
}

static bool consumes_both(char code)
{
    return strchr("M=X", code) != NULL;
}

static int64_t aligned_bases(const TvPaf *paf, int64_t start, int64_t end)
{
    int64_t query = paf->qstart;
    int64_t aligned = 0;

    for (size_t index = 0; index < paf->n_ops; index++) {
        TvCigarOp operation = paf->ops[index];

        if (consumes_both(operation.code)) {
            aligned += overlap(query, query + operation.length, start, end);
        }
        if (consumes_query(operation.code)) {
            query += operation.length;
        }
    }
    return aligned;
}

static int64_t inserted_bases(const TvPaf *paf, int64_t start, int64_t end)
{
    int64_t query = paf->qstart;
    int64_t inserted = 0;

    for (size_t index = 0; index < paf->n_ops; index++) {
        TvCigarOp operation = paf->ops[index];

        if (operation.code == 'I') {
            inserted += overlap(query, query + operation.length, start, end);
        }
        if (consumes_query(operation.code)) {
            query += operation.length;
        }
    }
    return inserted;
}

static double local_identity(const TvPaf *paf, int64_t start, int64_t end)
{
    int64_t query = paf->qstart;
    int64_t matches = 0;
    int64_t differences = 0;

    if (paf->identity_method == TV_IDENTITY_MISSING
        || paf->n_identity_ops == 0) {
        return NAN;
    }
    for (size_t index = 0; index < paf->n_identity_ops; index++) {
        TvCigarOp operation = paf->identity_ops[index];

        if (operation.code == '=' || operation.code == 'X') {
            int64_t contribution = overlap(query, query + operation.length,
                                           start, end);

            if (operation.code == '=') {
                matches += contribution;
            } else {
                differences += contribution;
            }
        }
        if (consumes_query(operation.code)) {
            query += operation.length;
        }
    }
    if (matches + differences == 0) {
        return NAN;
    }
    return (double)matches / (double)(matches + differences);
}

static int map_boundary(const TvPaf *paf, int64_t query_position,
                        int64_t *target_position)
{
    int64_t query = paf->qstart;
    int64_t target = paf->strand == '+' ? paf->tstart : paf->tend;

    if (query_position < paf->qstart || query_position > paf->qend) {
        return 0;
    }
    for (size_t index = 0; index < paf->n_ops; index++) {
        TvCigarOp operation = paf->ops[index];
        int64_t next_query = query
            + (consumes_query(operation.code) ? operation.length : 0);

        if (consumes_query(operation.code)
            && query_position >= query && query_position <= next_query) {
            if (consumes_both(operation.code)) {
                int64_t distance = query_position - query;

                *target_position = paf->strand == '+'
                    ? target + distance : target - distance;
            } else {
                *target_position = target;
            }
            return 1;
        }
        if (consumes_query(operation.code)) {
            query = next_query;
        }
        if (consumes_target(operation.code)) {
            target += paf->strand == '+' ? operation.length : -operation.length;
        }
    }
    if (query_position == paf->qend) {
        *target_position = target;
        return 1;
    }
    return 0;
}

static bool family_known(const char *family)
{
    return family != NULL && *family != '\0' && strcmp(family, ".") != 0
        && strcasecmp(family, "unknown") != 0;
}

static TvFamilyRelation family_relation(const TvTE *source, const TvTE *target)
{
    if (!family_known(source->family) || !family_known(target->family)) {
        return TV_FAMILY_UNKNOWN;
    }
    return strcasecmp(source->family, target->family) == 0
        ? TV_FAMILY_MATCH : TV_FAMILY_CONFLICT;
}

static int compare_paf(const void *left, const void *right)
{
    const TvPaf *a = left;
    const TvPaf *b = right;
    int comparison;

    if (a->query_genome != b->query_genome) {
        return a->query_genome < b->query_genome ? -1 : 1;
    }
    if (a->target_genome != b->target_genome) {
        return a->target_genome < b->target_genome ? -1 : 1;
    }
    comparison = strcmp(a->qname, b->qname);
    if (comparison != 0) {
        return comparison;
    }
    if (a->qstart != b->qstart) {
        return a->qstart < b->qstart ? -1 : 1;
    }
    if (a->qend != b->qend) {
        return a->qend < b->qend ? -1 : 1;
    }
    comparison = strcmp(a->tname, b->tname);
    if (comparison != 0) {
        return comparison;
    }
    if (a->tstart != b->tstart) {
        return a->tstart < b->tstart ? -1 : 1;
    }
    if (a->tend != b->tend) {
        return a->tend < b->tend ? -1 : 1;
    }
    if (a->strand != b->strand) {
        return a->strand < b->strand ? -1 : 1;
    }
    if (a->evidence_group_id != b->evidence_group_id) {
        return a->evidence_group_id < b->evidence_group_id ? -1 : 1;
    }
    if (a->origin != b->origin) {
        return a->origin < b->origin ? -1 : 1;
    }
    return 0;
}

static void push_projection(TvRun *run, TvProjection *projection)
{
    if (run->n_projections == run->cap_projections) {
        run->cap_projections = run->cap_projections == 0
            ? 256 : run->cap_projections * 2;
        run->projections = grow(run->projections, run->cap_projections,
                                sizeof(*run->projections));
    }
    run->projections[run->n_projections++] = *projection;
}

static void push_candidate(TvRun *run, TvCandidate *candidate)
{
    if (run->n_candidates == run->cap_candidates) {
        run->cap_candidates = run->cap_candidates == 0
            ? 512 : run->cap_candidates * 2;
        run->candidates = grow(run->candidates, run->cap_candidates,
                               sizeof(*run->candidates));
    }
    run->candidates[run->n_candidates++] = *candidate;
}

static bool candidate_precedes(const TvCandidate *left,
                               const TvCandidate *right)
{
    if (left->eligible != right->eligible) {
        return left->eligible;
    }
    if (left->score != right->score) {
        return left->score > right->score;
    }
    if (left->reciprocal_overlap != right->reciprocal_overlap) {
        return left->reciprocal_overlap > right->reciprocal_overlap;
    }
    if (left->breakpoint_distance != right->breakpoint_distance) {
        return left->breakpoint_distance < right->breakpoint_distance;
    }
    return left->candidate_id < right->candidate_id;
}

static int compare_candidate_order(const void *left, const void *right)
{
    const TvCandidate *a = left;
    const TvCandidate *b = right;

    if (candidate_precedes(a, b)) {
        return -1;
    }
    if (candidate_precedes(b, a)) {
        return 1;
    }
    return 0;
}

static void retain_candidate(TvCandidate **retained, size_t *count,
                             size_t *capacity, int limit,
                             const TvCandidate *candidate)
{
    size_t position;

    if (limit > 0 && *count == (size_t)limit
        && !candidate_precedes(candidate, &(*retained)[*count - 1])) {
        return;
    }
    if (limit == 0 || *count < (size_t)limit) {
        if (*count == *capacity) {
            size_t new_capacity = *capacity == 0 ? 16 : *capacity * 2;

            if (limit > 0 && new_capacity > (size_t)limit) {
                new_capacity = (size_t)limit;
            }
            *retained = grow(*retained, new_capacity, sizeof(**retained));
            *capacity = new_capacity;
        }
        position = (*count)++;
    } else {
        position = *count - 1;
    }
    if (limit == 0) {
        (*retained)[position] = *candidate;
        return;
    }
    while (position > 0
           && candidate_precedes(candidate, &(*retained)[position - 1])) {
        (*retained)[position] = (*retained)[position - 1];
        position--;
    }
    (*retained)[position] = *candidate;
}

static void push_decision(TvRun *run, TvDecision *decision)
{
    if (run->n_decisions == run->cap_decisions) {
        run->cap_decisions = run->cap_decisions == 0
            ? 256 : run->cap_decisions * 2;
        run->decisions = grow(run->decisions, run->cap_decisions,
                              sizeof(*run->decisions));
    }
    run->decisions[run->n_decisions++] = *decision;
}

static int breakpoint_distance(const TvProjection *projection, const TvTE *target)
{
    int64_t start_distance = projection->start >= target->start
        ? projection->start - target->start : target->start - projection->start;
    int64_t end_distance = projection->end >= target->end
        ? projection->end - target->end : target->end - projection->end;

    if (start_distance >= INT32_MAX || end_distance >= INT32_MAX
        || start_distance > INT32_MAX - end_distance) {
        return INT32_MAX;
    }
    return (int)(start_distance + end_distance);
}

static double reciprocal_overlap(const TvProjection *projection,
                                 const TvTE *target)
{
    int64_t projected_length = projection->end - projection->start;
    int64_t target_length = target->end - target->start;
    int64_t shared;

    if (projected_length <= 0 || target_length <= 0) {
        return 0.0;
    }
    shared = overlap(projection->start, projection->end,
                     target->start, target->end);
    return fmin((double)shared / (double)projected_length,
                (double)shared / (double)target_length);
}

static double candidate_score(const TvRun *run, const TvTE *source,
                              const TvTE *target,
                              const TvProjection *projection,
                              double reciprocal, double boundary,
                              TvFamilyRelation relation)
{
    double sum = 0.0;
    double observed_weight = 0.0;
    double family_value = relation == TV_FAMILY_MATCH
        ? 1.0 : relation == TV_FAMILY_CONFLICT ? 0.0 : 0.5;
    double completeness;

    (void)run;
    (void)source;
    (void)target;
    if (projection->left_flank_observed && projection->right_flank_observed) {
        sum += 0.25 * fmin(projection->left_flank, projection->right_flank);
        observed_weight += 0.25;
    }
    if (!isnan(projection->identity)) {
        sum += 0.20 * projection->identity;
        observed_weight += 0.20;
    }
    sum += 0.20 * reciprocal;
    observed_weight += 0.20;
    sum += 0.15 * boundary;
    observed_weight += 0.15;
    if (projection->mapq_observed) {
        sum += 0.10 * clamp((double)projection->mapq / 60.0, 0.0, 1.0);
        observed_weight += 0.10;
    }
    sum += 0.10 * family_value;
    observed_weight += 0.10;
    if (observed_weight == 0.0) {
        return 0.0;
    }
    completeness = observed_weight;
    return 100.0 * (sum / observed_weight) * (0.75 + 0.25 * completeness);
}

static void set_claimability(TvProjection *projection)
{
    projection->claim_type = TV_CLAIM_NONE;
    projection->claimable = false;

    if (projection->technical_state != TV_TECH_CALLABLE) {
        if (projection->technical_state == TV_TECH_GAP) {
            (void)snprintf(projection->claimability_reason,
                           sizeof(projection->claimability_reason),
                           "TECHNICAL_GAP");
        } else if (projection->technical_state == TV_TECH_AMBIGUOUS) {
            (void)snprintf(projection->claimability_reason,
                           sizeof(projection->claimability_reason),
                           "TECHNICAL_AMBIGUITY");
        } else {
            (void)snprintf(projection->claimability_reason,
                           sizeof(projection->claimability_reason),
                           "TECHNICAL_UNCALLABLE");
        }
        return;
    }
    if (projection->mapping_confidence != TV_MAPPING_MAPQ_PASS
        && projection->mapping_confidence != TV_MAPPING_UNIQUE_ALIGNMENT
        && projection->mapping_confidence != TV_MAPPING_UNIQUE_COPY_CONTEXT) {
        (void)snprintf(projection->claimability_reason,
                       sizeof(projection->claimability_reason),
                       "MAPPING_CONFIDENCE_NOT_ESTABLISHED");
        return;
    }
    if (projection->biological_state == TV_BIO_PRESENT
        && projection->annotation_state == TV_ANN_MATCHED) {
        projection->claim_type = TV_CLAIM_PRESENCE;
        projection->claimable = true;
        (void)snprintf(projection->claimability_reason,
                       sizeof(projection->claimability_reason),
                       "PASS_ANNOTATED_PRESENCE");
        return;
    }
    if (projection->biological_state == TV_BIO_PRESENT
        && projection->annotation_state == TV_ANN_FAMILY_CONFLICT) {
        projection->claim_type = TV_CLAIM_ANNOTATION_DISCORDANCE;
        projection->claimable = true;
        (void)snprintf(projection->claimability_reason,
                       sizeof(projection->claimability_reason),
                       "PASS_ANNOTATION_DISCORDANCE");
        return;
    }
    if (isnan(projection->identity)) {
        (void)snprintf(projection->claimability_reason,
                       sizeof(projection->claimability_reason),
                       "LOCAL_IDENTITY_MISSING");
        return;
    }
    if (projection->identity < 0.50) {
        (void)snprintf(projection->claimability_reason,
                       sizeof(projection->claimability_reason),
                       "LOCAL_IDENTITY_LT_0.50");
        return;
    }
    if (projection->biological_state == TV_BIO_PRESENT
        && projection->annotation_state == TV_ANN_MISSING
        && projection->te_aligned_fraction >= 0.50) {
        projection->claim_type = TV_CLAIM_PRESENCE;
        projection->claimable = true;
        (void)snprintf(projection->claimability_reason,
                       sizeof(projection->claimability_reason),
                       "PASS_UNANNOTATED_SEQUENCE_PRESENCE");
    } else if (projection->biological_state == TV_BIO_EMPTY
               && projection->insertion_fraction >= 0.70) {
        projection->claim_type = TV_CLAIM_EMPTY_SITE;
        projection->claimable = true;
        (void)snprintf(projection->claimability_reason,
                       sizeof(projection->claimability_reason),
                       "PASS_PAIRED_FLANK_EMPTY_SITE");
    } else if (projection->biological_state == TV_BIO_STRUCTURAL_ALTERNATIVE) {
        projection->claim_type = TV_CLAIM_STRUCTURAL_ALTERNATIVE;
        projection->claimable = true;
        (void)snprintf(projection->claimability_reason,
                       sizeof(projection->claimability_reason),
                       "PASS_STRUCTURAL_ALTERNATIVE");
    } else {
        (void)snprintf(projection->claimability_reason,
                       sizeof(projection->claimability_reason),
                       "BIOLOGICAL_EVIDENCE_INSUFFICIENT");
    }
}

static TvState legacy_state(const TvProjection *projection)
{
    if (projection->technical_state == TV_TECH_GAP) {
        return TV_ASSEMBLY_GAP;
    }
    if (projection->technical_state == TV_TECH_AMBIGUOUS) {
        return TV_PROJECTION_AMBIGUOUS;
    }
    if (projection->technical_state == TV_TECH_UNCALLABLE) {
        return TV_UNCALLABLE;
    }
    if (projection->biological_state == TV_BIO_PRESENT
        && projection->annotation_state == TV_ANN_MATCHED) {
        return TV_PRESENT_ANNOTATED;
    }
    if (projection->biological_state == TV_BIO_PRESENT
        && projection->annotation_state == TV_ANN_FAMILY_CONFLICT) {
        return TV_FAMILY_OR_BOUNDARY_DISCORDANCE;
    }
    if (projection->biological_state == TV_BIO_PRESENT
        && projection->annotation_state == TV_ANN_MISSING
        && projection->claimable) {
        return TV_PRESENT_UNANNOTATED;
    }
    if (projection->biological_state == TV_BIO_EMPTY && projection->claimable) {
        return TV_EMPTY_SITE_CONFIRMED;
    }
    if (projection->biological_state == TV_BIO_STRUCTURAL_ALTERNATIVE
        && projection->claimable) {
        return TV_STRUCTURAL_ALTERNATIVE;
    }
    return TV_UNCALLABLE;
}

static void calculate_quality(TvProjection *projection, double uniqueness)
{
    double sum = 0.0;
    double observed_weight = 0.0;

    if (projection->left_flank_observed && projection->right_flank_observed) {
        sum += 0.35 * fmin(projection->left_flank, projection->right_flank);
        observed_weight += 0.35;
    }
    if (!isnan(projection->identity)) {
        sum += 0.25 * projection->identity;
        observed_weight += 0.25;
    }
    if (projection->mapq_observed) {
        sum += 0.15 * clamp((double)projection->mapq / 60.0, 0.0, 1.0);
        observed_weight += 0.15;
    }
    if (!isnan(projection->n_fraction)) {
        sum += 0.15 * (1.0 - projection->n_fraction);
        observed_weight += 0.15;
    }
    if (!isnan(uniqueness)) {
        sum += 0.10 * uniqueness;
        observed_weight += 0.10;
    }
    projection->evidence_completeness = observed_weight;
    if (observed_weight == 0.0) {
        projection->quality = 0.0;
    } else {
        projection->quality = 100.0 * (sum / observed_weight)
            * (0.50 + 0.50 * observed_weight);
    }
}

static TvProjection evaluate(TvRun *run, int source_index, int target_genome,
                             int paf_index, int decision_index)
{
    TvTE *source = &run->nodes[source_index];
    TvPaf *paf = &run->pafs[paf_index];
    TvProjection projection = {
        .source_te = source_index,
        .target_genome = target_genome,
        .paf_index = paf_index,
        .decision_index = decision_index,
        .state = TV_UNCALLABLE,
        .technical_state = TV_TECH_UNCALLABLE,
        .biological_state = TV_BIO_UNKNOWN,
        .annotation_state = TV_ANN_UNKNOWN,
        .claim_type = TV_CLAIM_NONE,
        .mapping_confidence = paf->provider == TV_ALIGNMENT_PAF
            && paf->mapq_observed && paf->mapq >= run->cfg.min_mapq
            ? TV_MAPPING_MAPQ_PASS : TV_MAPPING_NOT_ESTABLISHED,
        .identity = NAN,
        .left_flank = NAN,
        .right_flank = NAN,
        .n_fraction = NAN,
        .mapq = paf->mapq,
        .mapq_observed = paf->mapq_observed,
        .contig = dupstr(paf->tname),
        .start = -1,
        .end = -1,
        .target_te = -1,
        .target_contig_index = paf->tcontig_index,
        .selected_candidate_index = -1,
        .candidate_start = run->n_candidates,
        .evidence_id = projection_id(run, source_index, target_genome,
                                     paf->evidence_group_id)
    };
    TvContig *source_contig = contig_at(&run->genomes[source->genome],
                                        source->contig_index);
    int flank_length = run->cfg.flank;
    int64_t left_start;
    int64_t right_end;
    int64_t mapped_start;
    int64_t mapped_end;
    TvContig *target_contig;
    int64_t n_start;
    int64_t n_end;
    int64_t candidate_low;
    int64_t candidate_high;
    TvCandidate *retained = NULL;
    size_t retained_count = 0;
    size_t retained_capacity = 0;
    int best_target = -1;
    uint64_t best_candidate_id = 0;
    TvCandidate best_candidate;
    bool have_best_candidate = false;
    double best_score = -1.0;
    double second_score = -1.0;
    double uniqueness = NAN;

    (void)snprintf(projection.decision_code, sizeof(projection.decision_code),
                   "UNINITIALIZED");
    (void)snprintf(projection.claimability_reason,
                   sizeof(projection.claimability_reason), "NOT_EVALUATED");
    if (source_contig == NULL) {
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code),
                       "SOURCE_CONTIG_NOT_INDEXED");
        calculate_quality(&projection, NAN);
        set_claimability(&projection);
        projection.state = legacy_state(&projection);
        return projection;
    }
    left_start = maximum_i64(0, source->start - flank_length);
    right_end = minimum_i64(source_contig->length, source->end + flank_length);
    projection.left_flank_observed = source->start - left_start == flank_length;
    projection.right_flank_observed = right_end - source->end == flank_length;
    if (projection.left_flank_observed) {
        projection.left_flank = clamp(
            (double)aligned_bases(paf, left_start, source->start)
            / (double)flank_length, 0.0, 1.0);
    }
    if (projection.right_flank_observed) {
        projection.right_flank = clamp(
            (double)aligned_bases(paf, source->end, right_end)
            / (double)flank_length, 0.0, 1.0);
    }
    projection.identity = local_identity(paf, left_start, right_end);
    projection.te_aligned_fraction = clamp(
        (double)aligned_bases(paf, source->start, source->end)
        / (double)(source->end - source->start), 0.0, 1.0);
    projection.insertion_fraction = clamp(
        (double)inserted_bases(paf, source->start, source->end)
        / (double)(source->end - source->start), 0.0, 1.0);
    if (!map_boundary(paf, source->start, &mapped_start)
        || !map_boundary(paf, source->end, &mapped_end)) {
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code),
                       "TE_BOUNDARY_UNMAPPED");
        calculate_quality(&projection, NAN);
        set_claimability(&projection);
        projection.state = legacy_state(&projection);
        return projection;
    }
    projection.start = minimum_i64(mapped_start, mapped_end);
    projection.end = maximum_i64(mapped_start, mapped_end);
    target_contig = contig_at(&run->genomes[target_genome],
                              projection.target_contig_index);
    if (target_contig == NULL) {
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code),
                       "TARGET_CONTIG_NOT_INDEXED");
        calculate_quality(&projection, NAN);
        set_claimability(&projection);
        projection.state = legacy_state(&projection);
        return projection;
    }
    n_start = maximum_i64(0, projection.start
                          - (projection.start == projection.end ? 10 : 0));
    n_end = minimum_i64(target_contig->length, projection.end
                        + (projection.start == projection.end ? 10 : 0));
    if (n_end <= n_start) {
        n_end = minimum_i64(target_contig->length, n_start + 1);
    }
    projection.n_fraction = n_fraction(run, target_genome,
                                       projection.target_contig_index,
                                       n_start, n_end);

    candidate_low = projection.start >= run->cfg.candidate_window
        ? projection.start - run->cfg.candidate_window : 0;
    candidate_high = projection.end
        <= target_contig->length - run->cfg.candidate_window
        ? projection.end + run->cfg.candidate_window : target_contig->length;
    TvIndexRange te_range = tv_te_index_range(
        run, target_genome, projection.target_contig_index,
        candidate_low, candidate_high);

    for (size_t interval_index = te_range.begin;
         interval_index < te_range.end; interval_index++) {
        int target_node = run->te_interval_index[interval_index].value;
        TvTE *target = &run->nodes[target_node];
        double reciprocal;
        int distance;
        double boundary;
        TvFamilyRelation relation;
        double score;
        TvCandidate candidate;

        run->performance.te_records_examined++;
        if (target->contig_index != projection.target_contig_index
            || target->end < candidate_low || target->start > candidate_high) {
            continue;
        }
        projection.nearby_candidate_count++;
        reciprocal = reciprocal_overlap(&projection, target);
        distance = breakpoint_distance(&projection, target);
        boundary = 1.0 - clamp(
            (double)distance
            / (double)(2 * run->cfg.candidate_window + 1), 0.0, 1.0);
        relation = family_relation(source, target);
        score = candidate_score(run, source, target, &projection,
                                reciprocal, boundary, relation);
        memset(&candidate, 0, sizeof(candidate));
        candidate.projection_index = (int)run->n_projections;
        candidate.target_te = target_node;
        candidate.score = score;
        candidate.reciprocal_overlap = reciprocal;
        candidate.boundary_score = boundary;
        candidate.breakpoint_distance = distance;
        candidate.family_relation = relation;
        candidate.context_relation = tv_candidate_context_relation(
            run, source_index, target_node,
            &candidate.shared_homology_group_id);
        candidate.context_compatible =
            candidate.context_relation != TV_CONTEXT_RELATION_CONFLICT
            && candidate.context_relation != TV_CONTEXT_RELATION_AMBIGUOUS;
        candidate.eligible = reciprocal >= run->cfg.min_reciprocal_overlap;
        candidate.candidate_id = candidate_id(run, projection.evidence_id,
                                              target_node);
        (void)snprintf(
            candidate.decision_code, sizeof(candidate.decision_code), "%s",
            !candidate.eligible ? "RECIPROCAL_OVERLAP_LT_0.50"
            : candidate.context_relation == TV_CONTEXT_RELATION_CONFLICT
                ? "BASE_SYNTENY_CONFLICT"
            : candidate.context_relation == TV_CONTEXT_RELATION_AMBIGUOUS
                ? "SYNTENY_CONTEXT_AMBIGUOUS"
            : "PASS_RECIPROCAL_OVERLAP");
        retain_candidate(&retained, &retained_count, &retained_capacity,
                         run->cfg.max_candidates, &candidate);
        if (!candidate.eligible) {
            continue;
        }
        projection.eligible_candidate_count++;
        if (!have_best_candidate
            || candidate_precedes(&candidate, &best_candidate)) {
            if (have_best_candidate && best_candidate.score > second_score) {
                second_score = best_candidate.score;
            }
            best_candidate = candidate;
            have_best_candidate = true;
            best_score = candidate.score;
            best_target = candidate.target_te;
            best_candidate_id = candidate.candidate_id;
        } else if (score > second_score) {
            second_score = score;
        }
    }
    if (retained_count > 1) {
        qsort(retained, retained_count, sizeof(*retained),
              compare_candidate_order);
    }
    for (size_t index = 0; index < retained_count; index++) {
        retained[index].rank = (int)index + 1;
        if (retained[index].candidate_id == best_candidate_id
            && retained[index].target_te == best_target) {
            retained[index].selected = true;
            projection.selected_candidate_index = (int)run->n_candidates;
        }
        push_candidate(run, &retained[index]);
    }
    free(retained);
    if (best_target >= 0) {
        projection.target_te = best_target;
        if (projection.eligible_candidate_count == 1) {
            uniqueness = 1.0;
        } else if (projection.eligible_candidate_count > 1) {
            uniqueness = clamp((best_score - second_score) / 20.0, 0.0, 1.0);
        }
    }
    projection.candidate_count = run->n_candidates - projection.candidate_start;
    calculate_quality(&projection, uniqueness);

    if (projection.target_te >= 0) {
        TvFamilyRelation relation = family_relation(
            source, &run->nodes[projection.target_te]);

        projection.biological_state = TV_BIO_PRESENT;
        projection.annotation_state = relation == TV_FAMILY_CONFLICT
            ? TV_ANN_FAMILY_CONFLICT : TV_ANN_MATCHED;
    } else if (projection.insertion_fraction >= 0.70) {
        projection.biological_state = TV_BIO_EMPTY;
        projection.annotation_state = TV_ANN_NOT_APPLICABLE;
    } else if (projection.te_aligned_fraction >= 0.50
               && !isnan(projection.identity)) {
        projection.biological_state = TV_BIO_PRESENT;
        projection.annotation_state = TV_ANN_MISSING;
    } else if (projection.end - projection.start
               < (source->end - source->start) / 2) {
        projection.biological_state = TV_BIO_STRUCTURAL_ALTERNATIVE;
        projection.annotation_state = TV_ANN_NOT_APPLICABLE;
    } else {
        projection.biological_state = TV_BIO_UNKNOWN;
        projection.annotation_state = TV_ANN_MISSING;
    }

    if (!projection.left_flank_observed || !projection.right_flank_observed) {
        projection.technical_state = TV_TECH_UNCALLABLE;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code),
                       "SOURCE_CONTIG_EDGE");
    } else if (paf->provider == TV_ALIGNMENT_PAF
               && !projection.mapq_observed) {
        projection.technical_state = TV_TECH_UNCALLABLE;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code), "MAPQ_MISSING");
    } else if (paf->provider == TV_ALIGNMENT_PAF
               && projection.mapq < run->cfg.min_mapq) {
        projection.technical_state = TV_TECH_UNCALLABLE;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code), "MAPQ_BELOW_THRESHOLD");
    } else if (paf->provider == TV_ALIGNMENT_MUMMER_DELTA
               && (isnan(paf->aggregate_identity)
                   || paf->aggregate_identity < run->cfg.min_delta_identity)) {
        projection.technical_state = TV_TECH_UNCALLABLE;
        projection.mapping_confidence = TV_MAPPING_FAILED;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code),
                       "DELTA_IDENTITY_BELOW_THRESHOLD");
    } else if (projection.left_flank < run->cfg.min_flank_fraction
               || projection.right_flank < run->cfg.min_flank_fraction) {
        projection.technical_state = TV_TECH_UNCALLABLE;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code),
                       "FLANK_SUPPORT_BELOW_THRESHOLD");
    } else if (projection.n_fraction > run->cfg.max_n_fraction) {
        projection.technical_state = TV_TECH_GAP;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code), "TARGET_SEQUENCE_GAP");
    } else if (projection.selected_candidate_index >= 0
               && run->candidates[projection.selected_candidate_index].context_relation
                  == TV_CONTEXT_RELATION_CONFLICT) {
        projection.technical_state = TV_TECH_AMBIGUOUS;
        projection.mapping_confidence = TV_MAPPING_AMBIGUOUS;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code),
                       "BASE_SYNTENY_CONFLICT");
    } else if (projection.selected_candidate_index >= 0
               && run->candidates[projection.selected_candidate_index].context_relation
                  == TV_CONTEXT_RELATION_AMBIGUOUS) {
        projection.technical_state = TV_TECH_AMBIGUOUS;
        projection.mapping_confidence = TV_MAPPING_AMBIGUOUS;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code),
                       "SYNTENY_CONTEXT_AMBIGUOUS");
    } else if (best_target >= 0
               && second_score >= best_score - run->cfg.near_best_delta) {
        projection.technical_state = TV_TECH_AMBIGUOUS;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code),
                       "MULTIPLE_NEAR_BEST_CANDIDATES");
    } else {
        projection.technical_state = TV_TECH_CALLABLE;
        if (projection.biological_state == TV_BIO_PRESENT
            && projection.annotation_state == TV_ANN_MATCHED) {
            (void)snprintf(projection.decision_code,
                           sizeof(projection.decision_code),
                           "ANNOTATED_LOCUS_MATCH");
        } else if (projection.annotation_state == TV_ANN_FAMILY_CONFLICT) {
            (void)snprintf(projection.decision_code,
                           sizeof(projection.decision_code),
                           "ANNOTATED_FAMILY_CONFLICT");
        } else if (projection.biological_state == TV_BIO_EMPTY) {
            (void)snprintf(projection.decision_code,
                           sizeof(projection.decision_code),
                           "TE_SPANS_QUERY_INSERTION");
        } else if (projection.biological_state == TV_BIO_PRESENT) {
            (void)snprintf(projection.decision_code,
                           sizeof(projection.decision_code),
                           "ALIGNED_TE_SEQUENCE_WITHOUT_ANNOTATION");
        } else if (projection.biological_state
                   == TV_BIO_STRUCTURAL_ALTERNATIVE) {
            (void)snprintf(projection.decision_code,
                           sizeof(projection.decision_code),
                           "NON_INSERTION_LENGTH_CHANGE");
        } else {
            (void)snprintf(projection.decision_code,
                           sizeof(projection.decision_code),
                           "LOCAL_IDENTITY_REQUIRED_FOR_SEQUENCE_CLAIM");
        }
    }
    set_claimability(&projection);
    projection.state = legacy_state(&projection);
    return projection;
}

static bool projection_precedes(const TvRun *run, int left_index,
                                int right_index)
{
    const TvProjection *left = &run->projections[left_index];
    const TvProjection *right = &run->projections[right_index];
    const TvPaf *left_paf = &run->pafs[left->paf_index];
    const TvPaf *right_paf = &run->pafs[right->paf_index];

    if (left->quality != right->quality) {
        return left->quality > right->quality;
    }
    if (left_paf->evidence_group_id != right_paf->evidence_group_id) {
        return left_paf->evidence_group_id < right_paf->evidence_group_id;
    }
    if (strcmp(left->contig, right->contig) != 0) {
        return strcmp(left->contig, right->contig) < 0;
    }
    if (left->start != right->start) {
        return left->start < right->start;
    }
    return left->end < right->end;
}

static bool incompatible_projection(const TvProjection *left,
                                    const TvProjection *right,
                                    int candidate_window)
{
    if (left->start < 0 || right->start < 0) {
        return true;
    }
    if (strcmp(left->contig, right->contig) != 0) {
        return true;
    }
    if (llabs(left->start - right->start) > candidate_window
        || llabs(left->end - right->end) > candidate_window) {
        return true;
    }
    if (left->target_te != right->target_te) {
        return true;
    }
    if (left->biological_state != right->biological_state
        || left->annotation_state != right->annotation_state) {
        return true;
    }
    return false;
}

static void finalize_decision(TvRun *run, TvDecision *decision)
{
    size_t start = decision->projection_start;
    size_t end = start + decision->projection_count;
    int winner = -1;

    if (decision->projection_count == 0) {
        decision->winner_projection = -1;
        decision->technical_state = TV_TECH_UNCALLABLE;
        decision->biological_state = TV_BIO_UNKNOWN;
        decision->annotation_state = TV_ANN_UNKNOWN;
        decision->claim_type = TV_CLAIM_NONE;
        decision->state = TV_UNCALLABLE;
        (void)snprintf(decision->decision_code,
                       sizeof(decision->decision_code),
                       "NO_SPANNING_ALIGNMENT");
        (void)snprintf(decision->claimability_reason,
                       sizeof(decision->claimability_reason),
                       "NO_SPANNING_ALIGNMENT");
        return;
    }
    for (size_t index = start; index < end; index++) {
        if (winner < 0 || projection_precedes(run, (int)index, winner)) {
            winner = (int)index;
        }
    }
    decision->winner_projection = winner;
    run->projections[winner].primary = true;
    for (size_t index = start; index < end; index++) {
        TvProjection *projection = &run->projections[index];
        int rank = 1;

        for (size_t other = start; other < end; other++) {
            if (other != index
                && projection_precedes(run, (int)other, (int)index)) {
                rank++;
            }
        }
        projection->rank = rank;
        if (projection->quality
            >= run->projections[winner].quality - run->cfg.near_best_delta) {
            projection->near_best = true;
            decision->near_best_count++;
        }
    }
    for (size_t left = start; left < end; left++) {
        if (!run->projections[left].near_best) {
            continue;
        }
        for (size_t right = left + 1; right < end; right++) {
            if (run->projections[right].near_best
                && incompatible_projection(&run->projections[left],
                                           &run->projections[right],
                                           run->cfg.candidate_window)) {
                decision->ambiguous = true;
            }
        }
    }
    if (run->projections[winner].technical_state == TV_TECH_AMBIGUOUS) {
        decision->ambiguous = true;
    }
    if (decision->ambiguous) {
        decision->technical_state = TV_TECH_AMBIGUOUS;
        decision->biological_state = TV_BIO_UNKNOWN;
        decision->annotation_state = TV_ANN_UNKNOWN;
        decision->claim_type = TV_CLAIM_NONE;
        decision->claimable = false;
        decision->state = TV_PROJECTION_AMBIGUOUS;
        decision->quality = run->projections[winner].quality * 0.75;
        decision->evidence_completeness =
            run->projections[winner].evidence_completeness;
        (void)snprintf(decision->decision_code,
                       sizeof(decision->decision_code),
                       "MULTIPLE_NEAR_BEST_PROJECTIONS");
        (void)snprintf(decision->claimability_reason,
                       sizeof(decision->claimability_reason),
                       "TECHNICAL_AMBIGUITY");
        for (size_t index = start; index < end; index++) {
            if (run->projections[index].near_best) {
                run->projections[index].decision_ambiguous = true;
                run->projections[index].mapping_confidence =
                    TV_MAPPING_AMBIGUOUS;
            }
        }
    } else {
        TvProjection *projection = &run->projections[winner];
        TvPaf *alignment = &run->pafs[projection->paf_index];

        if (alignment->provider == TV_ALIGNMENT_MUMMER_DELTA
            && projection->technical_state == TV_TECH_CALLABLE
            && projection->mapping_confidence
               == TV_MAPPING_NOT_ESTABLISHED) {
            TvContextRelation relation = TV_CONTEXT_RELATION_UNKNOWN;

            if (projection->selected_candidate_index >= 0) {
                relation = run->candidates[
                    projection->selected_candidate_index].context_relation;
            }
            projection->mapping_confidence =
                relation == TV_CONTEXT_RELATION_SUPPORTED
                ? TV_MAPPING_UNIQUE_COPY_CONTEXT
                : TV_MAPPING_UNIQUE_ALIGNMENT;
            set_claimability(projection);
            projection->state = legacy_state(projection);
        }

        decision->technical_state = projection->technical_state;
        decision->biological_state = projection->biological_state;
        decision->annotation_state = projection->annotation_state;
        decision->claim_type = projection->claim_type;
        decision->claimable = projection->claimable;
        decision->state = projection->state;
        decision->quality = projection->quality;
        decision->evidence_completeness = projection->evidence_completeness;
        (void)snprintf(decision->decision_code,
                       sizeof(decision->decision_code), "%s",
                       projection->decision_code);
        (void)snprintf(decision->claimability_reason,
                       sizeof(decision->claimability_reason), "%s",
                       projection->claimability_reason);
    }
}

static void build_projections(TvRun *run)
{
    for (size_t source_index = 0; source_index < run->n_nodes; source_index++) {
        TvTE *source = &run->nodes[source_index];

        for (size_t target = 0; target < run->n_genomes; target++) {
            TvDecision decision;
            int current_decision = (int)run->n_decisions;

            if ((int)target == source->genome) {
                continue;
            }
            memset(&decision, 0, sizeof(decision));
            decision.source_te = (int)source_index;
            decision.target_genome = (int)target;
            decision.projection_start = run->n_projections;
            decision.winner_projection = -1;
            decision.state = TV_UNCALLABLE;
            decision.technical_state = TV_TECH_UNCALLABLE;
            decision.biological_state = TV_BIO_UNKNOWN;
            decision.annotation_state = TV_ANN_UNKNOWN;
            decision.claim_type = TV_CLAIM_NONE;
            decision.decision_id = decision_id(run, (int)source_index,
                                               (int)target);
            TvIndexRange paf_range = tv_paf_index_range(
                run, source->genome, (int)target, source->contig_index,
                source->start, source->end);

            for (size_t interval_index = paf_range.begin;
                 interval_index < paf_range.end; interval_index++) {
                int paf_index = run->paf_interval_index[interval_index].value;
                TvPaf *paf = &run->pafs[paf_index];
                TvProjection projection;

                run->performance.paf_records_examined++;
                if (paf->query_genome != source->genome
                    || paf->target_genome != (int)target
                    || paf->qcontig_index != source->contig_index
                    || paf->qend <= source->start
                    || paf->qstart >= source->end) {
                    continue;
                }
                projection = evaluate(run, (int)source_index, (int)target,
                                      paf_index, current_decision);
                push_projection(run, &projection);
            }
            decision.projection_count = run->n_projections
                - decision.projection_start;
            finalize_decision(run, &decision);
            push_decision(run, &decision);
        }
    }
}

static TvCandidate *selected_candidate(TvRun *run, int projection_index)
{
    TvProjection *projection = &run->projections[projection_index];

    if (projection->selected_candidate_index >= 0
        && projection->selected_candidate_index < (int)run->n_candidates) {
        return &run->candidates[projection->selected_candidate_index];
    }
    return NULL;
}

static void append_edge(TvRun *run, int left, int right, double score,
                        bool compatible, int distance)
{
    int a;
    int b;

    if (left == right) {
        return;
    }
    a = left < right ? left : right;
    b = left < right ? right : left;
    if (run->n_edges == run->cap_edges) {
        run->cap_edges = run->cap_edges == 0 ? 256 : run->cap_edges * 2;
        run->edges = grow(run->edges, run->cap_edges, sizeof(*run->edges));
    }
    run->edges[run->n_edges++] = (TvEdge){
        .a = a,
        .b = b,
        .score = score,
        .family_compatible = compatible,
        .breakpoint_distance = distance
    };
}

static int compare_edge_pair(const void *left, const void *right)
{
    const TvEdge *a = left;
    const TvEdge *b = right;

    if (a->a != b->a) {
        return a->a < b->a ? -1 : 1;
    }
    if (a->b != b->b) {
        return a->b < b->b ? -1 : 1;
    }
    if (a->score != b->score) {
        return a->score > b->score ? -1 : 1;
    }
    if (a->breakpoint_distance != b->breakpoint_distance) {
        return a->breakpoint_distance < b->breakpoint_distance ? -1 : 1;
    }
    return 0;
}

static void reduce_edges(TvRun *run)
{
    size_t write = 0;

    if (run->n_edges > 1) {
        qsort(run->edges, run->n_edges, sizeof(*run->edges), compare_edge_pair);
    }
    for (size_t index = 0; index < run->n_edges; index++) {
        TvEdge *edge = &run->edges[index];

        if (write > 0 && run->edges[write - 1].a == edge->a
            && run->edges[write - 1].b == edge->b) {
            run->edges[write - 1].family_compatible =
                run->edges[write - 1].family_compatible
                && edge->family_compatible;
            continue;
        }
        if (write != index) {
            run->edges[write] = *edge;
        }
        write++;
    }
    run->n_edges = write;
}

static int find_edge(const TvRun *run, int a, int b)
{
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
    if (low < run->n_edges && run->edges[low].a == a
        && run->edges[low].b == b) {
        return (int)low;
    }
    return -1;
}

static void push_edge_support(TvRun *run, const TvEdgeSupport *support)
{
    if (run->n_edge_support == run->cap_edge_support) {
        run->cap_edge_support = run->cap_edge_support == 0
            ? 256 : run->cap_edge_support * 2;
        run->edge_support = grow(run->edge_support, run->cap_edge_support,
                                 sizeof(*run->edge_support));
    }
    run->edge_support[run->n_edge_support++] = *support;
}

static int compare_edge_support(const void *left, const void *right)
{
    const TvEdgeSupport *a = left;
    const TvEdgeSupport *b = right;

    if (a->a != b->a) {
        return a->a < b->a ? -1 : 1;
    }
    if (a->b != b->b) {
        return a->b < b->b ? -1 : 1;
    }
    if (a->evidence_group_id != b->evidence_group_id) {
        return a->evidence_group_id < b->evidence_group_id ? -1 : 1;
    }
    if (a->evidence_id != b->evidence_id) {
        return a->evidence_id < b->evidence_id ? -1 : 1;
    }
    if (a->direction_ab != b->direction_ab) {
        return a->direction_ab ? -1 : 1;
    }
    if (a->native != b->native) {
        return a->native ? -1 : 1;
    }
    return 0;
}

static void aggregate_edge_support(TvRun *run)
{
    size_t cursor = 0;

    for (size_t index = 0; index < run->n_projections; index++) {
        TvProjection *projection = &run->projections[index];
        TvDecision *decision = &run->decisions[projection->decision_index];
        TvPaf *alignment = &run->pafs[projection->paf_index];
        int a;
        int b;

        if (decision->ambiguous || !projection->near_best
            || projection->technical_state != TV_TECH_CALLABLE
            || projection->target_te < 0) {
            continue;
        }
        a = projection->source_te < projection->target_te
            ? projection->source_te : projection->target_te;
        b = projection->source_te < projection->target_te
            ? projection->target_te : projection->source_te;
        if (find_edge(run, a, b) < 0) {
            continue;
        }
        TvEdgeSupport support = {
            .a = a,
            .b = b,
            .evidence_group_id = alignment->evidence_group_id,
            .evidence_id = projection->evidence_id,
            .direction_ab = projection->source_te == a,
            .native = alignment->origin == TV_EVIDENCE_NATIVE
        };

        push_edge_support(run, &support);
    }
    if (run->n_edge_support > 1) {
        qsort(run->edge_support, run->n_edge_support,
              sizeof(*run->edge_support), compare_edge_support);
    }
    run->performance.edge_support_records = run->n_edge_support;
    for (size_t edge_index = 0; edge_index < run->n_edges; edge_index++) {
        TvEdge *edge = &run->edges[edge_index];
        uint64_t previous_group = 0;
        bool have_previous_group = false;

        while (cursor < run->n_edge_support
               && (run->edge_support[cursor].a < edge->a
                   || (run->edge_support[cursor].a == edge->a
                       && run->edge_support[cursor].b < edge->b))) {
            cursor++;
        }
        edge->support_start = cursor;
        while (cursor < run->n_edge_support
               && run->edge_support[cursor].a == edge->a
               && run->edge_support[cursor].b == edge->b) {
            TvEdgeSupport *current = &run->edge_support[cursor];

            if (!have_previous_group
                || current->evidence_group_id != previous_group) {
                edge->support_count++;
                previous_group = current->evidence_group_id;
                have_previous_group = true;
            }
            cursor++;
        }
        edge->support_record_count = cursor - edge->support_start;
        uint64_t first_ab = 0;
        uint64_t first_ba = 0;
        bool have_ab = false;
        bool have_ba = false;
        bool alternate_ab = false;
        bool alternate_ba = false;

        for (size_t support_index = edge->support_start;
             support_index < edge->support_start + edge->support_record_count;
             support_index++) {
            const TvEdgeSupport *support = &run->edge_support[support_index];

            if (!support->native) {
                continue;
            }
            if (support->direction_ab) {
                if (!have_ab) {
                    first_ab = support->evidence_group_id;
                    have_ab = true;
                } else if (support->evidence_group_id != first_ab) {
                    alternate_ab = true;
                }
            } else if (!have_ba) {
                first_ba = support->evidence_group_id;
                have_ba = true;
            } else if (support->evidence_group_id != first_ba) {
                alternate_ba = true;
            }
        }
        if (have_ab && have_ba
            && (first_ab != first_ba || alternate_ab || alternate_ba)) {
            edge->independent_reciprocal = true;
            edge->score = clamp(edge->score + 5.0, 0.0, 100.0);
        }
    }
}

static void build_edges(TvRun *run)
{
    for (size_t index = 0; index < run->n_decisions; index++) {
        TvDecision *decision = &run->decisions[index];
        TvProjection *projection;
        TvCandidate *candidate;

        if (decision->ambiguous || decision->winner_projection < 0
            || decision->technical_state != TV_TECH_CALLABLE) {
            continue;
        }
        projection = &run->projections[decision->winner_projection];
        if (projection->target_te < 0
            || projection->biological_state != TV_BIO_PRESENT) {
            continue;
        }
        candidate = selected_candidate(run, decision->winner_projection);
        if (candidate == NULL || !candidate->eligible) {
            continue;
        }
        append_edge(run, projection->source_te, projection->target_te,
                    candidate->score,
                    candidate->family_relation != TV_FAMILY_CONFLICT,
                    candidate->breakpoint_distance);
    }
    reduce_edges(run);
    aggregate_edge_support(run);
}

static int root(int *parents, int node)
{
    if (parents[node] != node) {
        parents[node] = root(parents, parents[node]);
    }
    return parents[node];
}

static bool node_has_copy_context(const TvRun *run, int node)
{
    if (run->te_context_offsets == NULL || node < 0
        || (size_t)node >= run->n_nodes) {
        return false;
    }
    for (size_t index = run->te_context_offsets[node];
         index < run->te_context_offsets[node + 1]; index++) {
        if (run->te_contexts[index].assignment == TV_CONTEXT_BRACKETED
            && run->te_contexts[index].context_index >= 0
            && (size_t)run->te_contexts[index].context_index < run->n_contexts
            && strcmp(run->contexts[run->te_contexts[index].context_index].status,
                      "PASS") == 0) {
            return true;
        }
    }
    return false;
}

static bool quota_allows(const TvRun *run, int *parents, int left_root,
                         int right_root)
{
    for (size_t genome = 0; genome < run->n_genomes; genome++) {
        int count = 0;

        for (size_t node = 0; node < run->n_nodes; node++) {
            int component = root(parents, (int)node);

            if ((component == left_root || component == right_root)
                && run->nodes[node].genome == (int)genome
                && !node_has_copy_context(run, (int)node)) {
                count++;
            }
        }
        if (count > run->genomes[genome].max_locus_copies) {
            return false;
        }
    }
    return true;
}

static bool component_families_compatible(const TvRun *run, int *parents,
                                          int left_root, int right_root)
{
    for (size_t left = 0; left < run->n_nodes; left++) {
        if (root(parents, (int)left) != left_root
            || !family_known(run->nodes[left].family)) {
            continue;
        }
        for (size_t right = 0; right < run->n_nodes; right++) {
            if (root(parents, (int)right) == right_root
                && family_known(run->nodes[right].family)
                && strcasecmp(run->nodes[left].family,
                              run->nodes[right].family) != 0) {
                return false;
            }
        }
    }
    return true;
}

typedef struct {
    int index;
    double score;
    int a;
    int b;
} TvEdgeOrder;

static int compare_edge_order(const void *left, const void *right)
{
    const TvEdgeOrder *a = left;
    const TvEdgeOrder *b = right;

    if (a->score != b->score) {
        return a->score > b->score ? -1 : 1;
    }
    if (a->a != b->a) {
        return a->a < b->a ? -1 : 1;
    }
    if (a->b != b->b) {
        return a->b < b->b ? -1 : 1;
    }
    return a->index == b->index ? 0 : (a->index < b->index ? -1 : 1);
}

static void select_graph(TvRun *run)
{
    int *parents = grow(NULL, run->n_nodes, sizeof(*parents));
    int *ranks = calloc(run->n_nodes, sizeof(*ranks));
    TvEdgeOrder *order = grow(NULL, run->n_edges, sizeof(*order));
    int *component_map;

    if (ranks == NULL) {
        tv_print_error("out of memory");
        exit(EXIT_FAILURE);
    }
    for (size_t index = 0; index < run->n_nodes; index++) {
        parents[index] = (int)index;
    }
    for (size_t index = 0; index < run->n_edges; index++) {
        order[index] = (TvEdgeOrder){
            .index = (int)index,
            .score = run->edges[index].score,
            .a = run->edges[index].a,
            .b = run->edges[index].b
        };
        (void)snprintf(run->edges[index].selection_reason,
                       sizeof(run->edges[index].selection_reason),
                       "NOT_EVALUATED");
    }
    if (run->n_edges > 1) {
        qsort(order, run->n_edges, sizeof(*order), compare_edge_order);
    }
    for (size_t index = 0; index < run->n_edges; index++) {
        TvEdge *edge = &run->edges[order[index].index];
        int left_root;
        int right_root;

        if (edge->score < run->cfg.min_edge_score) {
            (void)snprintf(edge->selection_reason,
                           sizeof(edge->selection_reason),
                           "BELOW_EDGE_THRESHOLD");
            continue;
        }
        if (!edge->family_compatible) {
            (void)snprintf(edge->selection_reason,
                           sizeof(edge->selection_reason),
                           "DIRECT_FAMILY_CONFLICT");
            continue;
        }
        left_root = root(parents, edge->a);
        right_root = root(parents, edge->b);
        if (left_root == right_root) {
            edge->selected = true;
            (void)snprintf(edge->selection_reason,
                           sizeof(edge->selection_reason), "CYCLE_SUPPORT");
            continue;
        }
        if (!component_families_compatible(run, parents, left_root, right_root)) {
            (void)snprintf(edge->selection_reason,
                           sizeof(edge->selection_reason),
                           "COMPONENT_FAMILY_CONFLICT");
            continue;
        }
        const char *context_reason = tv_context_merge_reason(
            run, parents, left_root, right_root);

        if (context_reason != NULL) {
            (void)snprintf(edge->selection_reason,
                           sizeof(edge->selection_reason), "%s",
                           context_reason);
            continue;
        }
        if (!quota_allows(run, parents, left_root, right_root)) {
            (void)snprintf(edge->selection_reason,
                           sizeof(edge->selection_reason), "COPY_QUOTA");
            continue;
        }
        if (ranks[left_root] < ranks[right_root]) {
            parents[left_root] = right_root;
        } else {
            parents[right_root] = left_root;
            if (ranks[left_root] == ranks[right_root]) {
                ranks[left_root]++;
            }
        }
        edge->selected = true;
        (void)snprintf(edge->selection_reason,
                       sizeof(edge->selection_reason), "SELECTED");
    }
    component_map = grow(NULL, run->n_nodes, sizeof(*component_map));
    run->components = grow(NULL, run->n_nodes, sizeof(*run->components));
    for (size_t index = 0; index < run->n_nodes; index++) {
        component_map[index] = -1;
    }
    run->n_loci = 0;
    for (size_t index = 0; index < run->n_nodes; index++) {
        int component_root = root(parents, (int)index);

        if (component_map[component_root] < 0) {
            component_map[component_root] = run->n_loci++;
        }
        run->components[index] = component_map[component_root];
    }
    free(component_map);
    free(order);
    free(ranks);
    free(parents);
}

int tv_analyze(TvRun *run)
{
    size_t node_count = 0;
    size_t node = 0;

    if (run->n_genomes < 2
        || (run->n_pafs == 0 && run->n_synteny_blocks == 0)
        || run->nodes != NULL) {
        tv_print_error(
            "analysis requires genomes plus alignment or synteny evidence and runs once");
        return -1;
    }
    if (run->n_genomes > (size_t)INT_MAX
        || run->n_pafs > (size_t)INT_MAX
        || run->n_synteny_blocks > (size_t)INT_MAX) {
        tv_print_error("input contains too many indexed records");
        return -1;
    }
    if (run->n_pafs > 1) {
        qsort(run->pafs, run->n_pafs, sizeof(*run->pafs), compare_paf);
    }
    for (size_t genome = 0; genome < run->n_genomes; genome++) {
        if (node_count > SIZE_MAX - run->genomes[genome].n_tes) {
            return -1;
        }
        node_count += run->genomes[genome].n_tes;
    }
    if (node_count == 0) {
        tv_print_error("no TE annotations were loaded");
        return -1;
    }
    if (node_count > (size_t)INT_MAX) {
        tv_print_error("too many TE nodes for graph indexing");
        return -1;
    }
    run->n_nodes = node_count;
    run->nodes = grow(NULL, run->n_nodes, sizeof(*run->nodes));
    for (size_t genome = 0; genome < run->n_genomes; genome++) {
        run->genomes[genome].node_offset = node;
        for (size_t te = 0; te < run->genomes[genome].n_tes; te++) {
            run->nodes[node++] = run->genomes[genome].tes[te];
        }
    }
    if (tv_build_synteny_contexts(run) != 0
        || tv_build_interval_indexes(run) != 0) {
        return -1;
    }
    build_projections(run);
    build_edges(run);
    select_graph(run);
    if (run->cfg.verbose) {
        fprintf(stderr,
                "tevox: %zu TEs, %zu evidence observations, %zu candidates, "
                "%zu decisions, %zu edges, %zu copy contexts, %d loci\n",
                run->n_nodes, run->n_projections, run->n_candidates,
                run->n_decisions, run->n_edges, run->n_contexts, run->n_loci);
    }
    return 0;
}
