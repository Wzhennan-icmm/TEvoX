#include "tevox.h"

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

static TvContig *find_contig(TvGenome *genome, const char *name)
{
    for (size_t index = 0; index < genome->n_contigs; index++) {
        if (strcmp(genome->contigs[index].name, name) == 0) {
            return &genome->contigs[index];
        }
    }
    return NULL;
}

static double n_fraction(TvGenome *genome, const char *name,
                         int64_t start, int64_t end)
{
    TvContig *contig = find_contig(genome, name);
    FILE *stream;
    int64_t position;
    int64_t bases = 0;
    int64_t ambiguous = 0;

    if (contig == NULL || contig->line_bases == 0 || start < 0
        || end > contig->length || end <= start) {
        return 1.0;
    }
    stream = fopen(genome->fasta_path, "rb");
    if (stream == NULL) {
        return 1.0;
    }
    position = start;
    while (position < end) {
        int64_t row = position / contig->line_bases;
        int64_t column = position % contig->line_bases;
        int64_t take = minimum_i64(contig->line_bases - column, end - position);
        int64_t offset = contig->seq_offset + row * contig->line_bytes + column;

        if (fseeko(stream, (off_t)offset, SEEK_SET) != 0) {
            fclose(stream);
            return 1.0;
        }
        for (int64_t index = 0; index < take; index++) {
            int character = fgetc(stream);

            if (character == EOF) {
                fclose(stream);
                return 1.0;
            }
            bases++;
            if (character == 'N' || character == 'n' || character == '-'
                || character == '.') {
                ambiguous++;
            }
        }
        position += take;
    }
    fclose(stream);
    return bases == 0 ? 1.0 : (double)ambiguous / (double)bases;
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

static int node_for(const TvRun *run, int genome, size_t local_index)
{
    size_t offset = 0;

    for (int index = 0; index < genome; index++) {
        offset += run->genomes[index].n_tes;
    }
    return (int)(offset + local_index);
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
    int64_t distance = llabs(projection->start - target->start)
        + llabs(projection->end - target->end);

    return distance > INT32_MAX ? INT32_MAX : (int)distance;
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
        .evidence_id = projection_id(run, source_index, target_genome,
                                     paf->evidence_group_id)
    };
    TvContig *source_contig = find_contig(&run->genomes[source->genome],
                                          source->contig);
    int flank_length = run->cfg.flank;
    int64_t left_start;
    int64_t right_end;
    int64_t mapped_start;
    int64_t mapped_end;
    TvContig *target_contig;
    int64_t n_start;
    int64_t n_end;
    int best_candidate = -1;
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
    target_contig = find_contig(&run->genomes[target_genome], projection.contig);
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
    projection.n_fraction = n_fraction(&run->genomes[target_genome],
                                       projection.contig, n_start, n_end);

    for (size_t local = 0; local < run->genomes[target_genome].n_tes; local++) {
        TvTE *target = &run->genomes[target_genome].tes[local];
        int target_node;
        double reciprocal;
        int distance;
        double boundary;
        TvFamilyRelation relation;
        double score;
        TvCandidate candidate;

        if (strcmp(target->contig, projection.contig) != 0
            || target->end < projection.start - run->cfg.candidate_window
            || target->start > projection.end + run->cfg.candidate_window) {
            continue;
        }
        projection.nearby_candidate_count++;
        target_node = node_for(run, target_genome, local);
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
        candidate.eligible = reciprocal >= run->cfg.min_reciprocal_overlap;
        candidate.candidate_id = candidate_id(run, projection.evidence_id,
                                              target_node);
        (void)snprintf(candidate.decision_code,
                       sizeof(candidate.decision_code), "%s",
                       candidate.eligible ? "PASS_RECIPROCAL_OVERLAP"
                       : "RECIPROCAL_OVERLAP_LT_0.50");
        push_candidate(run, &candidate);
        if (!candidate.eligible) {
            continue;
        }
        projection.eligible_candidate_count++;
        if (score > best_score
            || (score == best_score
                && (best_candidate < 0
                    || target_node
                       < run->candidates[(size_t)best_candidate].target_te))) {
            second_score = best_score;
            best_score = score;
            best_candidate = (int)(run->n_candidates - 1);
        } else if (score > second_score) {
            second_score = score;
        }
    }
    if (best_candidate >= 0) {
        TvCandidate *chosen = &run->candidates[(size_t)best_candidate];

        chosen->selected = true;
        projection.target_te = chosen->target_te;
        if (projection.eligible_candidate_count == 1) {
            uniqueness = 1.0;
        } else if (projection.eligible_candidate_count > 1) {
            uniqueness = clamp((best_score - second_score) / 20.0, 0.0, 1.0);
        }
    }
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
    } else if (!projection.mapq_observed) {
        projection.technical_state = TV_TECH_UNCALLABLE;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code), "MAPQ_MISSING");
    } else if (projection.mapq < run->cfg.min_mapq) {
        projection.technical_state = TV_TECH_UNCALLABLE;
        (void)snprintf(projection.decision_code,
                       sizeof(projection.decision_code), "MAPQ_BELOW_THRESHOLD");
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
    } else if (best_candidate >= 0
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
            }
        }
    } else {
        TvProjection *projection = &run->projections[winner];

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
            for (size_t paf_index = 0; paf_index < run->n_pafs; paf_index++) {
                TvPaf *paf = &run->pafs[paf_index];
                TvProjection projection;

                if (paf->query_genome != source->genome
                    || paf->target_genome != (int)target
                    || strcmp(paf->qname, source->contig) != 0
                    || paf->qend <= source->start
                    || paf->qstart >= source->end) {
                    continue;
                }
                projection = evaluate(run, (int)source_index, (int)target,
                                      (int)paf_index, current_decision);
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
    for (size_t index = 0; index < run->n_candidates; index++) {
        if (run->candidates[index].projection_index == projection_index
            && run->candidates[index].selected) {
            return &run->candidates[index];
        }
    }
    return NULL;
}

static void add_edge(TvRun *run, int left, int right, double score,
                     bool compatible, int distance)
{
    int a;
    int b;

    if (left == right) {
        return;
    }
    a = left < right ? left : right;
    b = left < right ? right : left;
    for (size_t index = 0; index < run->n_edges; index++) {
        TvEdge *edge = &run->edges[index];

        if (edge->a == a && edge->b == b) {
            edge->family_compatible = edge->family_compatible && compatible;
            if (score > edge->score) {
                edge->score = score;
                edge->breakpoint_distance = distance;
            }
            return;
        }
    }
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

static bool projection_supports_pair(const TvProjection *projection,
                                     int source, int target)
{
    return projection->near_best
        && projection->source_te == source
        && projection->target_te == target
        && projection->technical_state == TV_TECH_CALLABLE;
}

static void finalize_reciprocal(TvRun *run, TvEdge *edge)
{
    edge->support_count = 0;
    for (size_t index = 0; index < run->n_projections; index++) {
        TvProjection *projection = &run->projections[index];
        TvDecision *decision = &run->decisions[projection->decision_index];
        uint64_t group = run->pafs[projection->paf_index].evidence_group_id;
        bool supports = !decision->ambiguous
            && (projection_supports_pair(projection, edge->a, edge->b)
                || projection_supports_pair(projection, edge->b, edge->a));
        bool seen = false;

        if (!supports) {
            continue;
        }
        for (size_t previous = 0; previous < index; previous++) {
            TvProjection *prior = &run->projections[previous];
            TvDecision *prior_decision =
                &run->decisions[prior->decision_index];
            uint64_t previous_group =
                run->pafs[prior->paf_index].evidence_group_id;
            bool previous_supports = !prior_decision->ambiguous
                && (projection_supports_pair(prior, edge->a, edge->b)
                    || projection_supports_pair(prior, edge->b, edge->a));

            if (previous_supports && previous_group == group) {
                seen = true;
                break;
            }
        }
        if (!seen) {
            edge->support_count++;
        }
    }
    for (size_t left = 0; left < run->n_projections; left++) {
        TvProjection *projection_ab = &run->projections[left];
        TvDecision *decision_ab =
            &run->decisions[projection_ab->decision_index];
        TvPaf *paf_ab = &run->pafs[projection_ab->paf_index];

        if (decision_ab->ambiguous
            || paf_ab->origin != TV_EVIDENCE_NATIVE
            || !projection_supports_pair(projection_ab, edge->a, edge->b)) {
            continue;
        }
        for (size_t right = 0; right < run->n_projections; right++) {
            TvProjection *projection_ba = &run->projections[right];
            TvDecision *decision_ba =
                &run->decisions[projection_ba->decision_index];
            TvPaf *paf_ba = &run->pafs[projection_ba->paf_index];

            if (!decision_ba->ambiguous
                && paf_ba->origin == TV_EVIDENCE_NATIVE
                && projection_supports_pair(projection_ba, edge->b, edge->a)
                && paf_ab->evidence_group_id != paf_ba->evidence_group_id) {
                edge->independent_reciprocal = true;
                edge->score = clamp(edge->score + 5.0, 0.0, 100.0);
                return;
            }
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
        add_edge(run, projection->source_te, projection->target_te,
                 candidate->score,
                 candidate->family_relation != TV_FAMILY_CONFLICT,
                 candidate->breakpoint_distance);
    }
    for (size_t index = 0; index < run->n_edges; index++) {
        finalize_reciprocal(run, &run->edges[index]);
    }
}

static int root(int *parents, int node)
{
    if (parents[node] != node) {
        parents[node] = root(parents, parents[node]);
    }
    return parents[node];
}

static bool quota_allows(const TvRun *run, int *parents, int left_root,
                         int right_root)
{
    for (size_t genome = 0; genome < run->n_genomes; genome++) {
        int count = 0;

        for (size_t node = 0; node < run->n_nodes; node++) {
            int component = root(parents, (int)node);

            if ((component == left_root || component == right_root)
                && run->nodes[node].genome == (int)genome) {
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

static bool edge_precedes(const TvRun *run, int left_index, int right_index)
{
    const TvEdge *left = &run->edges[left_index];
    const TvEdge *right = &run->edges[right_index];

    if (left->score != right->score) {
        return left->score > right->score;
    }
    if (left->a != right->a) {
        return left->a < right->a;
    }
    return left->b < right->b;
}

static void select_graph(TvRun *run)
{
    int *parents = grow(NULL, run->n_nodes, sizeof(*parents));
    int *ranks = calloc(run->n_nodes, sizeof(*ranks));
    int *order = grow(NULL, run->n_edges, sizeof(*order));
    int *component_map;

    if (ranks == NULL) {
        tv_print_error("out of memory");
        exit(EXIT_FAILURE);
    }
    for (size_t index = 0; index < run->n_nodes; index++) {
        parents[index] = (int)index;
    }
    for (size_t index = 0; index < run->n_edges; index++) {
        order[index] = (int)index;
        (void)snprintf(run->edges[index].selection_reason,
                       sizeof(run->edges[index].selection_reason),
                       "NOT_EVALUATED");
    }
    for (size_t index = 1; index < run->n_edges; index++) {
        int value = order[index];
        size_t position = index;

        while (position > 0
               && edge_precedes(run, value, order[position - 1])) {
            order[position] = order[position - 1];
            position--;
        }
        order[position] = value;
    }
    for (size_t index = 0; index < run->n_edges; index++) {
        TvEdge *edge = &run->edges[order[index]];
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

    if (run->n_genomes < 2 || run->n_pafs == 0 || run->nodes != NULL) {
        tv_print_error("analysis requires genomes/alignments and runs once");
        return -1;
    }
    qsort(run->pafs, run->n_pafs, sizeof(*run->pafs), compare_paf);
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
    run->n_nodes = node_count;
    run->nodes = grow(NULL, run->n_nodes, sizeof(*run->nodes));
    for (size_t genome = 0; genome < run->n_genomes; genome++) {
        for (size_t te = 0; te < run->genomes[genome].n_tes; te++) {
            run->nodes[node++] = run->genomes[genome].tes[te];
        }
    }
    build_projections(run);
    build_edges(run);
    select_graph(run);
    if (run->cfg.verbose) {
        fprintf(stderr,
                "tevox: %zu TEs, %zu evidence observations, %zu candidates, "
                "%zu decisions, %zu edges, %d loci\n",
                run->n_nodes, run->n_projections, run->n_candidates,
                run->n_decisions, run->n_edges, run->n_loci);
    }
    return 0;
}
