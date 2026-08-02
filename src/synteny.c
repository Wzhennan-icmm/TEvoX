#include "tevox.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>

typedef struct {
    char *source_id;
    char *collinearity_path;
    char *genes_path;
    char *wgd_node;
} TvSyntenySource;

typedef struct {
    int gene_a;
    int gene_b;
    int provider_rank;
    double reported_evalue;
} TvPendingAnchor;

typedef struct {
    bool active;
    int alignment_id;
    int source_record;
    int source_line;
    int declared_anchor_count;
    double reported_score;
    double reported_evalue;
    char *header_contig_a;
    char *header_contig_b;
    char orientation;
    TvPendingAnchor *anchors;
    size_t n_anchors;
    size_t cap_anchors;
} TvPendingBlock;

typedef struct {
    int block_index;
    bool side_b;
    int genome;
    int contig_index;
    int64_t start;
    int64_t end;
} TvBlockSide;

static int compare_block_side(const void *left, const void *right)
{
    const TvBlockSide *a = left;
    const TvBlockSide *b = right;

    if (a->genome != b->genome) {
        return a->genome < b->genome ? -1 : 1;
    }
    if (a->contig_index != b->contig_index) {
        return a->contig_index < b->contig_index ? -1 : 1;
    }
    if (a->start != b->start) {
        return a->start < b->start ? -1 : 1;
    }
    if (a->end != b->end) {
        return a->end < b->end ? -1 : 1;
    }
    if (a->block_index != b->block_index) {
        return a->block_index < b->block_index ? -1 : 1;
    }
    return a->side_b == b->side_b ? 0 : (a->side_b ? 1 : -1);
}

typedef struct {
    int old_index;
    uint64_t id;
} TvIdOrder;

static int64_t min_i64(int64_t left, int64_t right)
{
    return left < right ? left : right;
}

static int64_t max_i64(int64_t left, int64_t right)
{
    return left > right ? left : right;
}

static int64_t interval_overlap(int64_t a_start, int64_t a_end,
                                int64_t b_start, int64_t b_end)
{
    int64_t start = max_i64(a_start, b_start);
    int64_t end = min_i64(a_end, b_end);

    return end > start ? end - start : 0;
}

static bool known_text(const char *text)
{
    return text != NULL && *text != '\0' && strcmp(text, TEVOX_UNKNOWN) != 0;
}

static int parse_finite_double(const char *text, double *value)
{
    char *end = NULL;

    errno = 0;
    *value = strtod(text, &end);
    return errno == 0 && end != text && *tv_strip(end) == '\0'
        && isfinite(*value);
}

static uint64_t hash_separator(uint64_t hash)
{
    return tv_hash_bytes(hash, "\0", 1);
}

static uint64_t hash_i64_text(uint64_t hash, int64_t value)
{
    char buffer[64];

    (void)snprintf(buffer, sizeof(buffer), "%lld", (long long)value);
    return tv_hash_text(hash, buffer);
}

static uint64_t hash_u64_text(uint64_t hash, uint64_t value)
{
    char buffer[32];

    (void)snprintf(buffer, sizeof(buffer), "%016llx",
                   (unsigned long long)value);
    return tv_hash_text(hash, buffer);
}

static int gene_index_by_id(const TvRun *run, const char *id)
{
    size_t low = 0;
    size_t high = run->n_genes;

    while (low < high) {
        size_t middle = low + (high - low) / 2;
        int comparison = strcmp(run->genes[middle].id, id);

        if (comparison < 0) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    if (low >= run->n_genes || strcmp(run->genes[low].id, id) != 0) {
        return -1;
    }
    if ((low > 0 && strcmp(run->genes[low - 1].id, id) == 0)
        || (low + 1 < run->n_genes
            && strcmp(run->genes[low + 1].id, id) == 0)) {
        return -2;
    }
    return (int)low;
}

static int compare_gene(const void *left, const void *right)
{
    const TvGene *a = left;
    const TvGene *b = right;
    int comparison;

    if (a->genome != b->genome) {
        return a->genome < b->genome ? -1 : 1;
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

static bool same_gene_record(const TvGene *left, const TvGene *right)
{
    return left->genome == right->genome
        && strcmp(left->id, right->id) == 0
        && strcmp(left->contig, right->contig) == 0
        && left->start == right->start && left->end == right->end
        && left->strand == right->strand
        && strcmp(left->subgenome_id, right->subgenome_id) == 0
        && strcmp(left->haplotype_id, right->haplotype_id) == 0;
}

static int push_gene(TvRun *run, TvGene *gene, const char *path,
                     int line_number)
{
    (void)path;
    (void)line_number;
    if (run->n_genes == run->cap_genes) {
        run->cap_genes = run->cap_genes == 0 ? 128 : run->cap_genes * 2;
        run->genes = tv_grow(run->genes, run->cap_genes,
                             sizeof(*run->genes));
    }
    run->genes[run->n_genes++] = *gene;
    return 0;
}

static int compare_gene_global_id(const void *left, const void *right)
{
    const TvGene *a = left;
    const TvGene *b = right;
    int comparison = strcmp(a->id, b->id);

    if (comparison != 0) {
        return comparison;
    }
    if (a->genome != b->genome) {
        return a->genome < b->genome ? -1 : 1;
    }
    return compare_gene(left, right);
}

static void free_gene_fields(TvGene *gene)
{
    free(gene->id);
    free(gene->contig);
    free(gene->subgenome_id);
    free(gene->haplotype_id);
    memset(gene, 0, sizeof(*gene));
}

static int load_gene_table(TvRun *run, const char *path)
{
    FILE *stream = fopen(path, "r");
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length;
    int line_number = 0;
    bool saw_header = false;
    int records = 0;

    if (stream == NULL) {
        tv_print_error("cannot open synteny gene table '%s'", path);
        return -1;
    }
    while ((length = getline(&line, &capacity, stream)) >= 0) {
        char *row;
        char *fields[9] = {0};
        int count;
        int genome;
        int64_t start;
        int64_t end;
        TvContig *sequence;
        TvGene gene;

        (void)length;
        line_number++;
        row = tv_strip(line);
        if (*row == '\0' || *row == '#') {
            continue;
        }
        count = tv_split_tabs(row, fields, 9);
        if (!saw_header) {
            static const char *expected[] = {
                "genome_id", "gene_id", "contig", "start", "end", "strand",
                "subgenome_id", "haplotype_id"
            };

            if (count != 8) {
                tv_print_error("%s:%d expected eight-column gene header", path,
                               line_number);
                goto fail;
            }
            for (int index = 0; index < 8; index++) {
                if (strcmp(fields[index], expected[index]) != 0) {
                    tv_print_error("%s:%d invalid gene-table header", path,
                                   line_number);
                    goto fail;
                }
            }
            saw_header = true;
            continue;
        }
        if (count != 8 || *fields[0] == '\0' || *fields[1] == '\0'
            || *fields[2] == '\0' || *fields[5] == '\0'
            || *fields[6] == '\0' || *fields[7] == '\0') {
            tv_print_error("%s:%d expected eight nonempty gene fields", path,
                           line_number);
            goto fail;
        }
        genome = tv_genome_id(run, fields[0]);
        if (genome < 0) {
            tv_print_error("%s:%d unknown genome '%s'", path, line_number,
                           fields[0]);
            goto fail;
        }
        if (!tv_parse_i64(fields[3], &start)
            || !tv_parse_i64(fields[4], &end)) {
            tv_print_error("%s:%d invalid gene coordinates", path, line_number);
            goto fail;
        }
        sequence = tv_find_contig(&run->genomes[genome], fields[2]);
        if (sequence == NULL || start < 0 || end <= start
            || end > sequence->length) {
            tv_print_error("%s:%d gene interval outside FASTA", path,
                           line_number);
            goto fail;
        }
        if (!(fields[5][0] == '+' || fields[5][0] == '-'
              || fields[5][0] == '.') || fields[5][1] != '\0') {
            tv_print_error("%s:%d invalid gene strand", path, line_number);
            goto fail;
        }
        memset(&gene, 0, sizeof(gene));
        gene.genome = genome;
        gene.id = tv_dupstr(fields[1]);
        gene.contig = tv_dupstr(fields[2]);
        gene.contig_index = sequence->index;
        gene.start = start;
        gene.end = end;
        gene.strand = fields[5][0];
        gene.subgenome_id = tv_dupstr(fields[6]);
        gene.haplotype_id = tv_dupstr(fields[7]);
        if (push_gene(run, &gene, path, line_number) != 0) {
            free(gene.id);
            free(gene.contig);
            free(gene.subgenome_id);
            free(gene.haplotype_id);
            goto fail;
        }
        records++;
    }
    free(line);
    fclose(stream);
    if (!saw_header || records == 0) {
        tv_print_error("synteny gene table '%s' contains no records", path);
        return -1;
    }
    return 0;

fail:
    free(line);
    fclose(stream);
    return -1;
}

static int compare_source(const void *left, const void *right)
{
    const TvSyntenySource *a = left;
    const TvSyntenySource *b = right;
    int comparison = strcmp(a->source_id, b->source_id);

    if (comparison != 0) {
        return comparison;
    }
    comparison = strcmp(a->collinearity_path, b->collinearity_path);
    if (comparison != 0) {
        return comparison;
    }
    return strcmp(a->genes_path, b->genes_path);
}

static void free_sources(TvSyntenySource *sources, size_t count)
{
    for (size_t index = 0; index < count; index++) {
        free(sources[index].source_id);
        free(sources[index].collinearity_path);
        free(sources[index].genes_path);
        free(sources[index].wgd_node);
    }
    free(sources);
}

static int load_source_specs(const char *path, TvSyntenySource **output,
                             size_t *output_count)
{
    FILE *stream = fopen(path, "r");
    char *base;
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length;
    int line_number = 0;
    bool saw_header = false;
    TvSyntenySource *sources = NULL;
    size_t count = 0;
    size_t source_capacity = 0;

    if (stream == NULL) {
        tv_print_error("cannot open synteny sources table '%s'", path);
        return -1;
    }
    base = tv_directory(path);
    while ((length = getline(&line, &capacity, stream)) >= 0) {
        char *row;
        char *fields[6] = {0};
        int field_count;

        (void)length;
        line_number++;
        row = tv_strip(line);
        if (*row == '\0' || *row == '#') {
            continue;
        }
        field_count = tv_split_tabs(row, fields, 6);
        if (!saw_header) {
            static const char *expected[] = {
                "source_id", "format", "collinearity", "genes", "wgd_node"
            };

            if (field_count != 5) {
                tv_print_error("%s:%d expected five-column synteny header",
                               path, line_number);
                goto fail;
            }
            for (int index = 0; index < 5; index++) {
                if (strcmp(fields[index], expected[index]) != 0) {
                    tv_print_error("%s:%d invalid synteny-sources header", path,
                                   line_number);
                    goto fail;
                }
            }
            saw_header = true;
            continue;
        }
        if (field_count != 5 || *fields[0] == '\0' || *fields[1] == '\0'
            || *fields[2] == '\0' || *fields[3] == '\0'
            || *fields[4] == '\0') {
            tv_print_error("%s:%d expected five nonempty synteny fields", path,
                           line_number);
            goto fail;
        }
        if (strcasecmp(fields[1], "mcscanx") != 0) {
            tv_print_error("%s:%d unsupported synteny format '%s'", path,
                           line_number, fields[1]);
            goto fail;
        }
        if (count == source_capacity) {
            source_capacity = source_capacity == 0 ? 8 : source_capacity * 2;
            sources = tv_grow(sources, source_capacity, sizeof(*sources));
        }
        memset(&sources[count], 0, sizeof(sources[count]));
        sources[count].source_id = tv_dupstr(fields[0]);
        sources[count].collinearity_path = tv_resolve(base, fields[2]);
        sources[count].genes_path = tv_resolve(base, fields[3]);
        sources[count].wgd_node = tv_dupstr(fields[4]);
        count++;
    }
    free(base);
    free(line);
    fclose(stream);
    if (!saw_header || count == 0) {
        tv_print_error("synteny sources table '%s' contains no records", path);
        free_sources(sources, count);
        return -1;
    }
    qsort(sources, count, sizeof(*sources), compare_source);
    for (size_t index = 1; index < count; index++) {
        if (strcmp(sources[index - 1].source_id, sources[index].source_id) == 0) {
            tv_print_error("synteny source IDs must be unique: '%s'",
                           sources[index].source_id);
            free_sources(sources, count);
            return -1;
        }
    }
    *output = sources;
    *output_count = count;
    return 0;

fail:
    free(base);
    free(line);
    fclose(stream);
    free_sources(sources, count);
    return -1;
}

static void pending_block_reset(TvPendingBlock *block)
{
    free(block->header_contig_a);
    free(block->header_contig_b);
    free(block->anchors);
    memset(block, 0, sizeof(*block));
}

static int parse_alignment_header(const char *row, TvPendingBlock *block,
                                  const char *path, int line_number)
{
    int record;
    int anchors;
    double score;
    double evalue;
    char pair[1024];
    char orientation[32];
    char extra;
    char *separator;
    int parsed;

    parsed = sscanf(row,
                    "## Alignment %d: score=%lf e_value=%lf N=%d %1023s %31s %c",
                    &record, &score, &evalue, &anchors, pair, orientation,
                    &extra);
    if (parsed != 6 || record < 0 || anchors < 2 || !isfinite(score)
        || !isfinite(evalue) || score < 0.0 || evalue < 0.0) {
        tv_print_error("%s:%d malformed MCScanX alignment header", path,
                       line_number);
        return -1;
    }
    if (strcmp(orientation, "plus") != 0
        && strcmp(orientation, "minus") != 0) {
        tv_print_error("%s:%d invalid MCScanX orientation", path, line_number);
        return -1;
    }
    separator = strchr(pair, '&');
    if (separator == NULL || separator == pair || separator[1] == '\0'
        || strchr(separator + 1, '&') != NULL) {
        tv_print_error("%s:%d invalid MCScanX contig pair", path, line_number);
        return -1;
    }
    *separator = '\0';
    memset(block, 0, sizeof(*block));
    block->active = true;
    block->alignment_id = record;
    block->source_record = record;
    block->source_line = line_number;
    block->declared_anchor_count = anchors;
    block->reported_score = score;
    block->reported_evalue = evalue;
    block->header_contig_a = tv_dupstr(pair);
    block->header_contig_b = tv_dupstr(separator + 1);
    block->orientation = strcmp(orientation, "plus") == 0 ? '+' : '-';
    return 0;
}

static int compare_pending_anchor(const void *left, const void *right,
                                  void *run_pointer)
{
    const TvPendingAnchor *a = left;
    const TvPendingAnchor *b = right;
    const TvRun *run = run_pointer;
    const TvGene *gene_a = &run->genes[a->gene_a];
    const TvGene *gene_b = &run->genes[b->gene_a];

    if (gene_a->start != gene_b->start) {
        return gene_a->start < gene_b->start ? -1 : 1;
    }
    if (gene_a->end != gene_b->end) {
        return gene_a->end < gene_b->end ? -1 : 1;
    }
    return strcmp(gene_a->id, gene_b->id);
}

static void sort_pending_anchors(TvPendingAnchor *anchors, size_t count,
                                 const TvRun *run)
{
    /* The test sets are small, and insertion sort avoids nonportable qsort_r. */
    for (size_t index = 1; index < count; index++) {
        TvPendingAnchor value = anchors[index];
        size_t position = index;

        while (position > 0
               && compare_pending_anchor(&value, &anchors[position - 1],
                                         (void *)run) < 0) {
            anchors[position] = anchors[position - 1];
            position--;
        }
        anchors[position] = value;
    }
}

static void normalize_pending_sides(TvPendingBlock *pending,
                                    const TvRun *run)
{
    const TvGene *left = &run->genes[pending->anchors[0].gene_a];
    const TvGene *right = &run->genes[pending->anchors[0].gene_b];
    int comparison = strcmp(run->genomes[left->genome].id,
                            run->genomes[right->genome].id);

    if (comparison == 0) {
        comparison = strcmp(left->contig, right->contig);
    }
    if (comparison == 0 && left->start != right->start) {
        comparison = left->start < right->start ? -1 : 1;
    }
    if (comparison <= 0) {
        return;
    }
    for (size_t index = 0; index < pending->n_anchors; index++) {
        int temporary = pending->anchors[index].gene_a;

        pending->anchors[index].gene_a = pending->anchors[index].gene_b;
        pending->anchors[index].gene_b = temporary;
    }
    char *temporary_contig = pending->header_contig_a;

    pending->header_contig_a = pending->header_contig_b;
    pending->header_contig_b = temporary_contig;
}

static uint64_t normalized_block_hash(const TvRun *run,
                                      const TvSyntenySource *source,
                                      const TvPendingBlock *pending,
                                      const char *kind)
{
    uint64_t hash = UINT64_C(1469598103934665603);

    hash = tv_hash_text(hash, kind);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, source->source_id);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, source->wgd_node);
    hash = hash_separator(hash);
    hash = tv_hash_bytes(hash, &pending->orientation,
                         sizeof(pending->orientation));
    for (size_t index = 0; index < pending->n_anchors; index++) {
        const TvGene *a = &run->genes[pending->anchors[index].gene_a];
        const TvGene *b = &run->genes[pending->anchors[index].gene_b];

        hash = hash_separator(hash);
        hash = tv_hash_text(hash, run->genomes[a->genome].id);
        hash = hash_separator(hash);
        hash = tv_hash_text(hash, a->id);
        hash = hash_separator(hash);
        hash = tv_hash_text(hash, run->genomes[b->genome].id);
        hash = hash_separator(hash);
        hash = tv_hash_text(hash, b->id);
    }
    return hash;
}

static void push_synteny_block(TvRun *run, const TvSyntenyBlock *block)
{
    if (run->n_synteny_blocks == run->cap_synteny_blocks) {
        run->cap_synteny_blocks = run->cap_synteny_blocks == 0
            ? 32 : run->cap_synteny_blocks * 2;
        run->synteny_blocks = tv_grow(
            run->synteny_blocks, run->cap_synteny_blocks,
            sizeof(*run->synteny_blocks));
    }
    run->synteny_blocks[run->n_synteny_blocks++] = *block;
}

static void push_synteny_anchor(TvRun *run, const TvSyntenyAnchor *anchor)
{
    if (run->n_synteny_anchors == run->cap_synteny_anchors) {
        run->cap_synteny_anchors = run->cap_synteny_anchors == 0
            ? 128 : run->cap_synteny_anchors * 2;
        run->synteny_anchors = tv_grow(
            run->synteny_anchors, run->cap_synteny_anchors,
            sizeof(*run->synteny_anchors));
    }
    run->synteny_anchors[run->n_synteny_anchors++] = *anchor;
}

static int finalize_pending_block(TvRun *run,
                                  const TvSyntenySource *source,
                                  TvPendingBlock *pending)
{
    TvGene *first_a;
    TvGene *first_b;
    TvSyntenyBlock block;
    size_t block_index;
    int64_t previous_a = -1;
    int64_t previous_b = -1;

    if (!pending->active) {
        return 0;
    }
    if (pending->n_anchors != (size_t)pending->declared_anchor_count) {
        tv_print_error("%s:%d MCScanX N=%d but parsed %zu anchors",
                       source->collinearity_path, pending->source_line,
                       pending->declared_anchor_count, pending->n_anchors);
        return -1;
    }
    for (int expected = 0; expected < pending->declared_anchor_count;
         expected++) {
        bool found = false;

        for (size_t index = 0; index < pending->n_anchors; index++) {
            if (pending->anchors[index].provider_rank == expected) {
                found = true;
                break;
            }
        }
        if (!found) {
            tv_print_error("%s:%d MCScanX anchor ranks must cover 0..N-1",
                           source->collinearity_path, pending->source_line);
            return -1;
        }
    }
    normalize_pending_sides(pending, run);
    sort_pending_anchors(pending->anchors, pending->n_anchors, run);
    first_a = &run->genes[pending->anchors[0].gene_a];
    first_b = &run->genes[pending->anchors[0].gene_b];
    memset(&block, 0, sizeof(block));
    block.source_id = tv_dupstr(source->source_id);
    block.source_path = tv_dupstr(source->collinearity_path);
    block.source_record = pending->source_record;
    block.source_line = pending->source_line;
    block.genome_a = first_a->genome;
    block.genome_b = first_b->genome;
    block.contig_a = tv_dupstr(first_a->contig);
    block.contig_b = tv_dupstr(first_b->contig);
    block.contig_a_index = first_a->contig_index;
    block.contig_b_index = first_b->contig_index;
    block.start_a = first_a->start;
    block.end_a = first_a->end;
    block.start_b = first_b->start;
    block.end_b = first_b->end;
    block.orientation = pending->orientation;
    block.anchor_count = pending->declared_anchor_count;
    block.anchor_start = run->n_synteny_anchors;
    block.reported_score = pending->reported_score;
    block.reported_evalue = pending->reported_evalue;
    block.wgd_node = tv_dupstr(source->wgd_node);
    block.context_a = -1;
    block.context_b = -1;
    (void)snprintf(block.status, sizeof(block.status), "PASS");
    (void)snprintf(block.reason, sizeof(block.reason), "VALIDATED");

    for (size_t index = 0; index < pending->n_anchors; index++) {
        TvGene *a = &run->genes[pending->anchors[index].gene_a];
        TvGene *b = &run->genes[pending->anchors[index].gene_b];

        if (a->genome != block.genome_a || b->genome != block.genome_b
            || a->contig_index != block.contig_a_index
            || b->contig_index != block.contig_b_index) {
            tv_print_error("%s:%d MCScanX block mixes genome or contig sides",
                           source->collinearity_path, pending->source_line);
            goto fail;
        }
        if (strcmp(a->contig, pending->header_contig_a) != 0
            || strcmp(b->contig, pending->header_contig_b) != 0) {
            tv_print_error("%s:%d MCScanX header/anchor contigs disagree",
                           source->collinearity_path, pending->source_line);
            goto fail;
        }
        if (index > 0) {
            if (a->start <= previous_a
                || (pending->orientation == '+' && b->start <= previous_b)
                || (pending->orientation == '-' && b->start >= previous_b)) {
                tv_print_error("%s:%d anchor order contradicts orientation",
                               source->collinearity_path,
                               pending->source_line);
                goto fail;
            }
        }
        previous_a = a->start;
        previous_b = b->start;
        block.start_a = min_i64(block.start_a, a->start);
        block.end_a = max_i64(block.end_a, a->end);
        block.start_b = min_i64(block.start_b, b->start);
        block.end_b = max_i64(block.end_b, b->end);
    }
    if (block.genome_a == block.genome_b
        && block.contig_a_index == block.contig_b_index
        && block.start_a == block.start_b && block.end_a == block.end_b) {
        tv_print_error("%s:%d self-synteny block maps an interval to itself",
                       source->collinearity_path, pending->source_line);
        goto fail;
    }
    block.block_id = normalized_block_hash(run, source, pending,
                                            "MCSCANX_BLOCK");
    block.evidence_group_id = normalized_block_hash(
        run, source, pending, "MCSCANX_EVIDENCE_GROUP");
    block_index = run->n_synteny_blocks;
    push_synteny_block(run, &block);
    for (size_t index = 0; index < pending->n_anchors; index++) {
        TvSyntenyAnchor anchor;
        uint64_t hash = UINT64_C(1469598103934665603);
        const TvGene *a = &run->genes[pending->anchors[index].gene_a];
        const TvGene *b = &run->genes[pending->anchors[index].gene_b];

        memset(&anchor, 0, sizeof(anchor));
        hash = tv_hash_text(hash, "MCSCANX_ANCHOR");
        hash = hash_separator(hash);
        hash = hash_u64_text(hash, block.block_id);
        hash = hash_separator(hash);
        hash = tv_hash_text(hash, a->id);
        hash = hash_separator(hash);
        hash = tv_hash_text(hash, b->id);
        anchor.anchor_id = hash;
        anchor.block_index = (int)block_index;
        anchor.rank = (int)index + 1;
        anchor.provider_rank = pending->anchors[index].provider_rank;
        anchor.gene_a = pending->anchors[index].gene_a;
        anchor.gene_b = pending->anchors[index].gene_b;
        anchor.reported_evalue = pending->anchors[index].reported_evalue;
        push_synteny_anchor(run, &anchor);
    }
    return 0;

fail:
    free(block.source_id);
    free(block.source_path);
    free(block.contig_a);
    free(block.contig_b);
    free(block.wgd_node);
    return -1;
}

static int push_pending_anchor(TvPendingBlock *block, int gene_a, int gene_b,
                               int provider_rank, double evalue, const char *path,
                               int line_number)
{
    for (size_t index = 0; index < block->n_anchors; index++) {
        if (block->anchors[index].provider_rank == provider_rank
            || block->anchors[index].gene_a == gene_a
            || block->anchors[index].gene_b == gene_b) {
            tv_print_error("%s:%d duplicate MCScanX anchor rank or gene", path,
                           line_number);
            return -1;
        }
    }
    if (block->n_anchors == block->cap_anchors) {
        block->cap_anchors = block->cap_anchors == 0
            ? 16 : block->cap_anchors * 2;
        block->anchors = tv_grow(block->anchors, block->cap_anchors,
                                 sizeof(*block->anchors));
    }
    block->anchors[block->n_anchors++] = (TvPendingAnchor){
        .gene_a = gene_a,
        .gene_b = gene_b,
        .provider_rank = provider_rank,
        .reported_evalue = evalue
    };
    return 0;
}

static int parse_anchor_row(TvRun *run, TvPendingBlock *block,
                            const char *row, const char *path, int line_number)
{
    char gene_a[1024];
    char gene_b[1024];
    char evalue_text[128];
    char extra;
    double evalue;
    int provider_block;
    int provider_rank;
    int consumed = 0;
    int parsed;
    int index_a;
    int index_b;

    if (!block->active) {
        tv_print_error("%s:%d anchor outside an MCScanX block", path,
                       line_number);
        return -1;
    }
    parsed = sscanf(row, " %d- %d: %n", &provider_block, &provider_rank,
                    &consumed);
    if (parsed != 2 || consumed <= 0 || provider_block != block->alignment_id
        || provider_rank < 0) {
        tv_print_error("%s:%d invalid MCScanX anchor prefix", path,
                       line_number);
        return -1;
    }
    parsed = sscanf(row + consumed, " %1023s %1023s %127s %c", gene_a,
                    gene_b, evalue_text, &extra);
    if (parsed != 3 || !parse_finite_double(evalue_text, &evalue)
        || evalue < 0.0) {
        tv_print_error("%s:%d malformed MCScanX anchor", path, line_number);
        return -1;
    }
    index_a = gene_index_by_id(run, gene_a);
    index_b = gene_index_by_id(run, gene_b);
    if (index_a == -2 || index_b == -2) {
        tv_print_error("%s:%d ambiguous unqualified MCScanX gene ID", path,
                       line_number);
        return -1;
    }
    if (index_a < 0 || index_b < 0) {
        tv_print_error("%s:%d unknown MCScanX anchor gene", path, line_number);
        return -1;
    }
    if (index_a == index_b) {
        tv_print_error("%s:%d MCScanX anchor maps a gene to itself", path,
                       line_number);
        return -1;
    }
    return push_pending_anchor(block, index_a, index_b, provider_rank,
                               evalue, path, line_number);
}

static int load_mcscanx(TvRun *run, const TvSyntenySource *source)
{
    FILE *stream = fopen(source->collinearity_path, "r");
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length;
    int line_number = 0;
    int blocks = 0;
    int *alignment_ids = NULL;
    size_t n_alignment_ids = 0;
    size_t cap_alignment_ids = 0;
    TvPendingBlock pending;

    if (stream == NULL) {
        tv_print_error("cannot open MCScanX collinearity '%s'",
                       source->collinearity_path);
        return -1;
    }
    memset(&pending, 0, sizeof(pending));
    while ((length = getline(&line, &capacity, stream)) >= 0) {
        char *row;

        (void)length;
        line_number++;
        row = tv_strip(line);
        if (*row == '\0') {
            continue;
        }
        if (strncmp(row, "## Alignment ", 13) == 0) {
            if (pending.active) {
                if (finalize_pending_block(run, source, &pending) != 0) {
                    goto fail;
                }
                blocks++;
                pending_block_reset(&pending);
            }
            if (parse_alignment_header(row, &pending,
                                       source->collinearity_path,
                                       line_number) != 0) {
                goto fail;
            }
            for (size_t index = 0; index < n_alignment_ids; index++) {
                if (alignment_ids[index] == pending.alignment_id) {
                    tv_print_error("%s:%d duplicate MCScanX Alignment ID %d",
                                   source->collinearity_path, line_number,
                                   pending.alignment_id);
                    goto fail;
                }
            }
            if (n_alignment_ids == cap_alignment_ids) {
                cap_alignment_ids = cap_alignment_ids == 0
                    ? 16 : cap_alignment_ids * 2;
                alignment_ids = tv_grow(alignment_ids, cap_alignment_ids,
                                        sizeof(*alignment_ids));
            }
            alignment_ids[n_alignment_ids++] = pending.alignment_id;
        } else if (*row == '#') {
            continue;
        } else if (parse_anchor_row(run, &pending, row,
                                    source->collinearity_path,
                                    line_number) != 0) {
            goto fail;
        }
    }
    if (pending.active) {
        if (finalize_pending_block(run, source, &pending) != 0) {
            goto fail;
        }
        blocks++;
    }
    pending_block_reset(&pending);
    free(alignment_ids);
    free(line);
    fclose(stream);
    if (blocks == 0) {
        tv_print_error("MCScanX file '%s' contains no alignment blocks",
                       source->collinearity_path);
        return -1;
    }
    return 0;

fail:
    pending_block_reset(&pending);
    free(alignment_ids);
    free(line);
    fclose(stream);
    return -1;
}

static int compare_id_order(const void *left, const void *right)
{
    const TvIdOrder *a = left;
    const TvIdOrder *b = right;

    if (a->id != b->id) {
        return a->id < b->id ? -1 : 1;
    }
    return a->old_index == b->old_index
        ? 0 : (a->old_index < b->old_index ? -1 : 1);
}

static int canonicalize_blocks(TvRun *run)
{
    TvIdOrder *order;
    TvSyntenyBlock *blocks;
    TvSyntenyAnchor *anchors;
    size_t anchor_count = 0;

    if (run->n_synteny_blocks == 0) {
        return 0;
    }
    order = tv_grow(NULL, run->n_synteny_blocks, sizeof(*order));
    for (size_t index = 0; index < run->n_synteny_blocks; index++) {
        order[index] = (TvIdOrder){(int)index,
                                   run->synteny_blocks[index].block_id};
    }
    qsort(order, run->n_synteny_blocks, sizeof(*order), compare_id_order);
    for (size_t index = 1; index < run->n_synteny_blocks; index++) {
        if (order[index - 1].id == order[index].id) {
            tv_print_error("duplicate normalized MCScanX block");
            free(order);
            return -1;
        }
    }
    blocks = tv_grow(NULL, run->n_synteny_blocks, sizeof(*blocks));
    anchors = tv_grow(NULL, run->n_synteny_anchors, sizeof(*anchors));
    for (size_t index = 0; index < run->n_synteny_blocks; index++) {
        int old_index = order[index].old_index;

        blocks[index] = run->synteny_blocks[old_index];
        blocks[index].anchor_start = anchor_count;
        for (size_t anchor_index = 0;
             anchor_index < run->n_synteny_anchors; anchor_index++) {
            if (run->synteny_anchors[anchor_index].block_index != old_index) {
                continue;
            }
            anchors[anchor_count] = run->synteny_anchors[anchor_index];
            anchors[anchor_count].block_index = (int)index;
            anchor_count++;
        }
    }
    free(order);
    free(run->synteny_blocks);
    free(run->synteny_anchors);
    run->synteny_blocks = blocks;
    run->synteny_anchors = anchors;
    run->cap_synteny_blocks = run->n_synteny_blocks;
    run->cap_synteny_anchors = run->n_synteny_anchors;
    return anchor_count == run->n_synteny_anchors ? 0 : -1;
}

int tv_load_synteny_sources(TvRun *run, const char *path)
{
    TvSyntenySource *sources = NULL;
    size_t source_count = 0;

    if (run == NULL || path == NULL || *path == '\0' || run->n_genomes < 1
        || run->n_synteny_blocks != 0 || run->n_contexts != 0) {
        tv_print_error("synteny sources require loaded genomes and one load");
        return -1;
    }
    if (tv_register_input(run, "synteny_manifest", path) != 0) {
        return -1;
    }
    if (load_source_specs(path, &sources, &source_count) != 0) {
        return -1;
    }
    for (size_t index = 0; index < source_count; index++) {
        if (tv_register_input(run, "synteny_gene_table",
                              sources[index].genes_path) != 0
            || tv_register_input(run, "synteny_collinearity",
                                 sources[index].collinearity_path) != 0
            || load_gene_table(run, sources[index].genes_path) != 0) {
            free_sources(sources, source_count);
            return -1;
        }
    }
    qsort(run->genes, run->n_genes, sizeof(*run->genes), compare_gene);
    size_t write = 0;
    for (size_t index = 0; index < run->n_genes; index++) {
        if (write > 0
            && run->genes[write - 1].genome == run->genes[index].genome
            && strcmp(run->genes[write - 1].id, run->genes[index].id) == 0) {
            if (!same_gene_record(&run->genes[write - 1],
                                  &run->genes[index])) {
                tv_print_error("conflicting synteny gene key '%s:%s'",
                               run->genomes[run->genes[index].genome].id,
                               run->genes[index].id);
                free_sources(sources, source_count);
                return -1;
            }
            free_gene_fields(&run->genes[index]);
            continue;
        }
        if (write != index) {
            run->genes[write] = run->genes[index];
        }
        write++;
    }
    run->n_genes = write;
    qsort(run->genes, run->n_genes, sizeof(*run->genes),
          compare_gene_global_id);
    for (size_t index = 1; index < run->n_genes; index++) {
        if (strcmp(run->genes[index - 1].id, run->genes[index].id) == 0) {
            tv_print_error("MCScanX gene IDs must be globally unique: '%s'",
                           run->genes[index].id);
            free_sources(sources, source_count);
            return -1;
        }
    }
    for (size_t index = 0; index < source_count; index++) {
        if (load_mcscanx(run, &sources[index]) != 0) {
            free_sources(sources, source_count);
            return -1;
        }
        if (run->n_synteny_source_paths == run->cap_synteny_source_paths) {
            run->cap_synteny_source_paths = run->cap_synteny_source_paths == 0
                ? 8 : run->cap_synteny_source_paths * 2;
            run->synteny_source_paths = tv_grow(
                run->synteny_source_paths, run->cap_synteny_source_paths,
                sizeof(*run->synteny_source_paths));
            run->synteny_gene_paths = tv_grow(
                run->synteny_gene_paths, run->cap_synteny_source_paths,
                sizeof(*run->synteny_gene_paths));
        }
        run->synteny_source_paths[run->n_synteny_source_paths] =
            tv_dupstr(sources[index].collinearity_path);
        run->synteny_gene_paths[run->n_synteny_source_paths] =
            tv_dupstr(sources[index].genes_path);
        run->n_synteny_source_paths++;
    }
    free_sources(sources, source_count);
    return canonicalize_blocks(run);
}

static int parent_root(int *parents, int value)
{
    if (parents[value] != value) {
        parents[value] = parent_root(parents, parents[value]);
    }
    return parents[value];
}

static void parent_union(int *parents, int *ranks, int left, int right)
{
    int a = parent_root(parents, left);
    int b = parent_root(parents, right);

    if (a == b) {
        return;
    }
    if (ranks[a] < ranks[b]) {
        parents[a] = b;
    } else {
        parents[b] = a;
        if (ranks[a] == ranks[b]) {
            ranks[a]++;
        }
    }
}

static int side_gene(const TvSyntenyAnchor *anchor, bool side_b)
{
    return side_b ? anchor->gene_b : anchor->gene_a;
}

static size_t collect_side_genes(const TvRun *run, int block_index,
                                 bool side_b, int *genes)
{
    size_t count = 0;

    const TvSyntenyBlock *block = &run->synteny_blocks[block_index];
    size_t end = block->anchor_start + (size_t)block->anchor_count;

    for (size_t index = block->anchor_start; index < end; index++) {
        const TvSyntenyAnchor *anchor = &run->synteny_anchors[index];

        genes[count++] = side_gene(anchor, side_b);
    }
    return count;
}

static bool side_metadata(const TvRun *run, const TvBlockSide *side,
                          const char **subgenome, const char **haplotype)
{
    const TvSyntenyBlock *block = &run->synteny_blocks[side->block_index];
    bool subgenome_conflict = false;
    bool haplotype_conflict = false;

    *subgenome = NULL;
    *haplotype = NULL;
    size_t end = block->anchor_start + (size_t)block->anchor_count;

    for (size_t index = block->anchor_start; index < end; index++) {
        const TvSyntenyAnchor *anchor = &run->synteny_anchors[index];
        const TvGene *gene;
        gene = &run->genes[side_gene(anchor, side->side_b)];
        if (known_text(gene->subgenome_id)) {
            if (*subgenome != NULL
                && strcmp(*subgenome, gene->subgenome_id) != 0) {
                subgenome_conflict = true;
            } else if (!subgenome_conflict) {
                *subgenome = gene->subgenome_id;
            }
        }
        if (known_text(gene->haplotype_id)) {
            if (*haplotype != NULL
                && strcmp(*haplotype, gene->haplotype_id) != 0) {
                haplotype_conflict = true;
            } else if (!haplotype_conflict) {
                *haplotype = gene->haplotype_id;
            }
        }
    }
    if (subgenome_conflict) {
        *subgenome = NULL;
    }
    if (haplotype_conflict) {
        *haplotype = NULL;
    }
    return block->anchor_count >= 2;
}

static bool side_mergeable(const TvRun *run, const TvBlockSide *left,
                           const TvBlockSide *right)
{
    int64_t left_length;
    int64_t right_length;
    int64_t shared;
    int *left_genes;
    int *right_genes;
    size_t left_count;
    size_t right_count;
    size_t intersection = 0;
    size_t union_count;
    bool result;
    const TvSyntenyBlock *left_block;
    const TvSyntenyBlock *right_block;
    const char *left_subgenome;
    const char *right_subgenome;
    const char *left_haplotype;
    const char *right_haplotype;

    if (left->genome != right->genome
        || left->contig_index != right->contig_index
        || left->block_index == right->block_index) {
        return false;
    }
    left_block = &run->synteny_blocks[left->block_index];
    right_block = &run->synteny_blocks[right->block_index];
    if ((known_text(left_block->wgd_node)
         && known_text(right_block->wgd_node)
         && strcmp(left_block->wgd_node, right_block->wgd_node) != 0)
        || !side_metadata(run, left, &left_subgenome, &left_haplotype)
        || !side_metadata(run, right, &right_subgenome, &right_haplotype)
        || (left_subgenome != NULL && right_subgenome != NULL
            && strcmp(left_subgenome, right_subgenome) != 0)
        || (left_haplotype != NULL && right_haplotype != NULL
            && strcmp(left_haplotype, right_haplotype) != 0)) {
        return false;
    }
    left_length = left->end - left->start;
    right_length = right->end - right->start;
    shared = interval_overlap(left->start, left->end,
                              right->start, right->end);
    if (left_length <= 0 || right_length <= 0
        || (double)shared / (double)left_length < 0.80
        || (double)shared / (double)right_length < 0.80) {
        return false;
    }
    left_genes = tv_grow(NULL, (size_t)run->synteny_blocks[left->block_index]
                         .anchor_count, sizeof(*left_genes));
    right_genes = tv_grow(NULL, (size_t)run->synteny_blocks[right->block_index]
                          .anchor_count, sizeof(*right_genes));
    left_count = collect_side_genes(run, left->block_index, left->side_b,
                                    left_genes);
    right_count = collect_side_genes(run, right->block_index, right->side_b,
                                     right_genes);
    for (size_t a = 0; a < left_count; a++) {
        for (size_t b = 0; b < right_count; b++) {
            if (left_genes[a] == right_genes[b]) {
                intersection++;
                break;
            }
        }
    }
    union_count = left_count + right_count - intersection;
    result = union_count > 0
        && (double)intersection / (double)union_count >= 0.50;
    free(left_genes);
    free(right_genes);
    return result;
}

static void merge_metadata_value(char **current, bool *conflict,
                                 const char *value)
{
    if (!known_text(value)) {
        return;
    }
    if (*current == NULL) {
        *current = tv_dupstr(value);
    } else if (strcmp(*current, value) != 0) {
        *conflict = true;
    }
}

static uint64_t context_hash(const TvRun *run, int genome, int contig,
                             int64_t start, int64_t end,
                             const int *genes, size_t gene_count,
                             const char *subgenome, const char *haplotype,
                             const char *wgd_node)
{
    uint64_t hash = UINT64_C(1469598103934665603);

    hash = tv_hash_text(hash, "COPY_CONTEXT");
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, run->genomes[genome].id);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, run->genomes[genome].contigs[contig].name);
    hash = hash_separator(hash);
    hash = hash_i64_text(hash, start);
    hash = hash_separator(hash);
    hash = hash_i64_text(hash, end);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, subgenome);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, haplotype);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, wgd_node);
    for (size_t index = 0; index < gene_count; index++) {
        hash = hash_separator(hash);
        hash = tv_hash_text(hash, run->genes[genes[index]].id);
    }
    return hash;
}

static int compare_int_by_gene_id(const void *left, const void *right,
                                  void *run_pointer)
{
    int a = *(const int *)left;
    int b = *(const int *)right;
    const TvRun *run = run_pointer;

    return strcmp(run->genes[a].id, run->genes[b].id);
}

static void sort_gene_indices(int *values, size_t count, const TvRun *run)
{
    for (size_t index = 1; index < count; index++) {
        int value = values[index];
        size_t position = index;

        while (position > 0
               && compare_int_by_gene_id(&value, &values[position - 1],
                                         (void *)run) < 0) {
            values[position] = values[position - 1];
            position--;
        }
        values[position] = value;
    }
}

static bool append_unique_int(int **values, size_t *count, size_t *capacity,
                              int value)
{
    for (size_t index = 0; index < *count; index++) {
        if ((*values)[index] == value) {
            return false;
        }
    }
    if (*count == *capacity) {
        *capacity = *capacity == 0 ? 16 : *capacity * 2;
        *values = tv_grow(*values, *capacity, sizeof(**values));
    }
    (*values)[(*count)++] = value;
    return true;
}

static int build_contexts_from_sides(TvRun *run, TvBlockSide *sides,
                                     size_t side_count, int *side_parents)
{
    int *roots = tv_grow(NULL, side_count, sizeof(*roots));
    size_t root_count = 0;

    for (size_t side = 0; side < side_count; side++) {
        int root = parent_root(side_parents, (int)side);
        bool seen = false;

        for (size_t index = 0; index < root_count; index++) {
            if (roots[index] == root) {
                seen = true;
                break;
            }
        }
        if (!seen) {
            roots[root_count++] = root;
        }
    }
    for (size_t root_index = 0; root_index < root_count; root_index++) {
        int root = roots[root_index];
        TvCopyContext context;
        int *genes = NULL;
        size_t gene_count = 0;
        size_t gene_capacity = 0;
        char *subgenome = NULL;
        char *haplotype = NULL;
        char *wgd_node = NULL;
        bool subgenome_conflict = false;
        bool haplotype_conflict = false;
        bool wgd_node_conflict = false;
        bool initialized = false;

        memset(&context, 0, sizeof(context));
        for (size_t side = 0; side < side_count; side++) {
            TvBlockSide *item;
            TvSyntenyBlock *block;
            int *side_genes;
            size_t count;

            if (parent_root(side_parents, (int)side) != root) {
                continue;
            }
            item = &sides[side];
            block = &run->synteny_blocks[item->block_index];
            if (!initialized) {
                context.genome = item->genome;
                context.contig_index = item->contig_index;
                context.contig = tv_dupstr(
                    run->genomes[item->genome].contigs[item->contig_index].name);
                context.start = item->start;
                context.end = item->end;
                initialized = true;
            } else {
                context.start = min_i64(context.start, item->start);
                context.end = max_i64(context.end, item->end);
            }
            merge_metadata_value(&wgd_node, &wgd_node_conflict,
                                 block->wgd_node);
            side_genes = tv_grow(NULL, (size_t)block->anchor_count,
                                 sizeof(*side_genes));
            count = collect_side_genes(run, item->block_index, item->side_b,
                                       side_genes);
            for (size_t index = 0; index < count; index++) {
                TvGene *gene = &run->genes[side_genes[index]];

                (void)append_unique_int(&genes, &gene_count, &gene_capacity,
                                        side_genes[index]);
                merge_metadata_value(&subgenome, &subgenome_conflict,
                                     gene->subgenome_id);
                merge_metadata_value(&haplotype, &haplotype_conflict,
                                     gene->haplotype_id);
            }
            free(side_genes);
        }
        sort_gene_indices(genes, gene_count, run);
        context.subgenome_id = subgenome_conflict
            ? tv_dupstr(TEVOX_UNKNOWN)
            : (subgenome == NULL ? tv_dupstr(TEVOX_UNKNOWN) : subgenome);
        context.haplotype_id = haplotype_conflict
            ? tv_dupstr(TEVOX_UNKNOWN)
            : (haplotype == NULL ? tv_dupstr(TEVOX_UNKNOWN) : haplotype);
        context.wgd_node = wgd_node_conflict
            ? tv_dupstr(TEVOX_UNKNOWN)
            : (wgd_node == NULL ? tv_dupstr(TEVOX_UNKNOWN) : wgd_node);
        if (subgenome_conflict) {
            free(subgenome);
        }
        if (haplotype_conflict) {
            free(haplotype);
        }
        if (wgd_node_conflict) {
            free(wgd_node);
        }
        context.context_id = context_hash(
            run, context.genome, context.contig_index, context.start,
            context.end, genes, gene_count, context.subgenome_id,
            context.haplotype_id, context.wgd_node);
        context.anchor_count = (int)gene_count;
        (void)snprintf(context.status, sizeof(context.status), "%s",
                       subgenome_conflict || haplotype_conflict
                           || wgd_node_conflict
                       ? "METADATA_AMBIGUOUS" : "PASS");
        context.syntenic_copy_id = tv_dupstr(TEVOX_UNKNOWN);
        if (run->n_contexts == run->cap_contexts) {
            run->cap_contexts = run->cap_contexts == 0
                ? 32 : run->cap_contexts * 2;
            run->contexts = tv_grow(run->contexts, run->cap_contexts,
                                    sizeof(*run->contexts));
        }
        int context_index = (int)run->n_contexts;
        run->contexts[run->n_contexts++] = context;
        for (size_t side = 0; side < side_count; side++) {
            TvBlockSide *item;

            if (parent_root(side_parents, (int)side) != root) {
                continue;
            }
            item = &sides[side];
            if (item->side_b) {
                run->synteny_blocks[item->block_index].context_b = context_index;
            } else {
                run->synteny_blocks[item->block_index].context_a = context_index;
            }
        }
        free(genes);
    }
    free(roots);
    return 0;
}

static uint64_t homology_group_hash(const TvRun *run, int *parents, int root)
{
    uint64_t *ids = tv_grow(NULL, run->n_contexts, sizeof(*ids));
    size_t count = 0;
    uint64_t hash = UINT64_C(1469598103934665603);

    for (size_t index = 0; index < run->n_contexts; index++) {
        if (parent_root(parents, (int)index) == root) {
            ids[count++] = run->contexts[index].context_id;
        }
    }
    for (size_t index = 1; index < count; index++) {
        uint64_t value = ids[index];
        size_t position = index;

        while (position > 0 && value < ids[position - 1]) {
            ids[position] = ids[position - 1];
            position--;
        }
        ids[position] = value;
    }
    hash = tv_hash_text(hash, "HOMOLOGY_GROUP");
    for (size_t index = 0; index < count; index++) {
        hash = hash_separator(hash);
        hash = hash_u64_text(hash, ids[index]);
    }
    free(ids);
    return hash;
}

static int compare_context_index(const TvRun *run, int left, int right)
{
    const TvCopyContext *a = &run->contexts[left];
    const TvCopyContext *b = &run->contexts[right];
    int comparison;

    if (a->homology_group_id != b->homology_group_id) {
        return a->homology_group_id < b->homology_group_id ? -1 : 1;
    }
    if (a->genome != b->genome) {
        return a->genome < b->genome ? -1 : 1;
    }
    comparison = strcmp(a->contig, b->contig);
    if (comparison != 0) {
        return comparison;
    }
    if (a->start != b->start) {
        return a->start < b->start ? -1 : 1;
    }
    if (a->end != b->end) {
        return a->end < b->end ? -1 : 1;
    }
    return a->context_id == b->context_id
        ? 0 : (a->context_id < b->context_id ? -1 : 1);
}

static void assign_homology_groups_and_copies(TvRun *run)
{
    int *parents = tv_grow(NULL, run->n_contexts, sizeof(*parents));
    int *ranks = calloc(run->n_contexts, sizeof(*ranks));
    int *order = tv_grow(NULL, run->n_contexts, sizeof(*order));

    if (ranks == NULL) {
        tv_print_error("out of memory");
        exit(EXIT_FAILURE);
    }
    for (size_t index = 0; index < run->n_contexts; index++) {
        parents[index] = (int)index;
        order[index] = (int)index;
    }
    for (size_t index = 0; index < run->n_synteny_blocks; index++) {
        TvSyntenyBlock *block = &run->synteny_blocks[index];

        parent_union(parents, ranks, block->context_a, block->context_b);
    }
    for (size_t index = 0; index < run->n_contexts; index++) {
        int root = parent_root(parents, (int)index);

        run->contexts[index].homology_group_id =
            homology_group_hash(run, parents, root);
    }
    for (size_t index = 0; index < run->n_synteny_blocks; index++) {
        run->synteny_blocks[index].homology_group_id =
            run->contexts[run->synteny_blocks[index].context_a]
            .homology_group_id;
    }
    for (size_t index = 1; index < run->n_contexts; index++) {
        int value = order[index];
        size_t position = index;

        while (position > 0
               && compare_context_index(run, value, order[position - 1]) < 0) {
            order[position] = order[position - 1];
            position--;
        }
        order[position] = value;
    }
    uint64_t prior_group = 0;
    int prior_genome = -1;
    int copy = 0;
    for (size_t index = 0; index < run->n_contexts; index++) {
        TvCopyContext *context = &run->contexts[order[index]];
        char copy_id[32];

        if (index == 0 || context->homology_group_id != prior_group
            || context->genome != prior_genome) {
            copy = 1;
        } else {
            copy++;
        }
        prior_group = context->homology_group_id;
        prior_genome = context->genome;
        (void)snprintf(copy_id, sizeof(copy_id), "copy%03d", copy);
        free(context->syntenic_copy_id);
        context->syntenic_copy_id = tv_dupstr(copy_id);
    }
    free(order);
    free(ranks);
    free(parents);
}

static int compare_context_id_order(const void *left, const void *right)
{
    return compare_id_order(left, right);
}

static int canonicalize_contexts(TvRun *run)
{
    TvIdOrder *order = tv_grow(NULL, run->n_contexts, sizeof(*order));
    TvCopyContext *contexts = tv_grow(NULL, run->n_contexts,
                                      sizeof(*contexts));
    int *remap = tv_grow(NULL, run->n_contexts, sizeof(*remap));

    for (size_t index = 0; index < run->n_contexts; index++) {
        order[index] = (TvIdOrder){(int)index, run->contexts[index].context_id};
    }
    qsort(order, run->n_contexts, sizeof(*order), compare_context_id_order);
    for (size_t index = 1; index < run->n_contexts; index++) {
        if (order[index - 1].id == order[index].id) {
            tv_print_error("duplicate normalized copy context ID");
            free(remap);
            free(contexts);
            free(order);
            return -1;
        }
    }
    for (size_t index = 0; index < run->n_contexts; index++) {
        contexts[index] = run->contexts[order[index].old_index];
        remap[order[index].old_index] = (int)index;
    }
    for (size_t index = 0; index < run->n_synteny_blocks; index++) {
        run->synteny_blocks[index].context_a =
            remap[run->synteny_blocks[index].context_a];
        run->synteny_blocks[index].context_b =
            remap[run->synteny_blocks[index].context_b];
    }
    free(run->contexts);
    run->contexts = contexts;
    run->cap_contexts = run->n_contexts;
    free(remap);
    free(order);
    return 0;
}

static bool block_uses_context(const TvSyntenyBlock *block, int context,
                               bool *side_b)
{
    if (block->context_a == context) {
        *side_b = false;
        return true;
    }
    if (block->context_b == context) {
        *side_b = true;
        return true;
    }
    return false;
}

static int compare_gene_coordinate(const TvRun *run, int left, int right)
{
    const TvGene *a = &run->genes[left];
    const TvGene *b = &run->genes[right];

    if (a->start != b->start) {
        return a->start < b->start ? -1 : 1;
    }
    if (a->end != b->end) {
        return a->end < b->end ? -1 : 1;
    }
    return strcmp(a->id, b->id);
}

static void sort_genes_by_coordinate(int *genes, size_t count,
                                     const TvRun *run)
{
    for (size_t index = 1; index < count; index++) {
        int value = genes[index];
        size_t position = index;

        while (position > 0
               && compare_gene_coordinate(run, value,
                                          genes[position - 1]) < 0) {
            genes[position] = genes[position - 1];
            position--;
        }
        genes[position] = value;
    }
}

static bool bracket_in_context(const TvRun *run, int te_node, int context,
                               int *left_gene, int *right_gene)
{
    const TvTE *te = &run->nodes[te_node];

    for (size_t block_index = 0; block_index < run->n_synteny_blocks;
         block_index++) {
        const TvSyntenyBlock *block = &run->synteny_blocks[block_index];
        bool side_b;
        int *genes;
        size_t count;

        if (!block_uses_context(block, context, &side_b)) {
            continue;
        }
        genes = tv_grow(NULL, (size_t)block->anchor_count, sizeof(*genes));
        count = collect_side_genes(run, (int)block_index, side_b, genes);
        sort_genes_by_coordinate(genes, count, run);
        for (size_t index = 0; index + 1 < count; index++) {
            const TvGene *left = &run->genes[genes[index]];
            const TvGene *right = &run->genes[genes[index + 1]];

            if (left->end <= te->start && right->start >= te->end) {
                *left_gene = genes[index];
                *right_gene = genes[index + 1];
                free(genes);
                return true;
            }
        }
        free(genes);
    }
    return false;
}

static uint64_t te_context_hash(const TvRun *run, int te_node,
                                int context_index)
{
    const TvTE *te = &run->nodes[te_node];
    uint64_t hash = UINT64_C(1469598103934665603);

    hash = tv_hash_text(hash, "TE_CONTEXT");
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, run->genomes[te->genome].id);
    hash = hash_separator(hash);
    hash = tv_hash_text(hash, te->id);
    hash = hash_separator(hash);
    return hash_u64_text(hash, run->contexts[context_index].context_id);
}

static void push_te_context(TvRun *run, const TvTEContext *assignment)
{
    if (run->n_te_contexts == run->cap_te_contexts) {
        run->cap_te_contexts = run->cap_te_contexts == 0
            ? 64 : run->cap_te_contexts * 2;
        run->te_contexts = tv_grow(run->te_contexts,
                                   run->cap_te_contexts,
                                   sizeof(*run->te_contexts));
    }
    run->te_contexts[run->n_te_contexts++] = *assignment;
}

static void assign_te_contexts(TvRun *run)
{
    run->te_context_offsets = tv_grow(NULL, run->n_nodes + 1,
                                      sizeof(*run->te_context_offsets));
    for (size_t node = 0; node < run->n_nodes; node++) {
        const TvTE *te = &run->nodes[node];
        size_t assignment_start = run->n_te_contexts;
        int strong_count = 0;

        run->te_context_offsets[node] = assignment_start;

        for (size_t context_index = 0; context_index < run->n_contexts;
             context_index++) {
            const TvCopyContext *context = &run->contexts[context_index];
            int64_t shared;
            int left_gene = -1;
            int right_gene = -1;
            TvTEContext assignment;

            if (te->genome != context->genome
                || te->contig_index != context->contig_index) {
                continue;
            }
            shared = interval_overlap(te->start, te->end,
                                      context->start, context->end);
            if (shared <= 0) {
                continue;
            }
            memset(&assignment, 0, sizeof(assignment));
            assignment.te_node = (int)node;
            assignment.context_index = (int)context_index;
            assignment.overlap_fraction =
                (double)shared / (double)(te->end - te->start);
            assignment.left_anchor = -1;
            assignment.right_anchor = -1;
            if (bracket_in_context(run, (int)node, (int)context_index,
                                   &left_gene, &right_gene)) {
                assignment.assignment = TV_CONTEXT_BRACKETED;
                assignment.left_anchor = left_gene;
                assignment.right_anchor = right_gene;
                strong_count++;
            } else {
                assignment.assignment = TV_CONTEXT_BLOCK_INTERIOR;
            }
            assignment.te_context_id = te_context_hash(
                run, (int)node, (int)context_index);
            push_te_context(run, &assignment);
        }
        if (strong_count > 1) {
            for (size_t index = assignment_start;
                 index < run->n_te_contexts; index++) {
                if (run->te_contexts[index].assignment
                    == TV_CONTEXT_BRACKETED) {
                    run->te_contexts[index].assignment = TV_CONTEXT_AMBIGUOUS;
                }
            }
        }
        run->te_context_offsets[node + 1] = run->n_te_contexts;
    }
    if (run->n_te_contexts > 0) {
        run->context_scratch_left = tv_grow(
            NULL, run->n_te_contexts, sizeof(*run->context_scratch_left));
        run->context_scratch_right = tv_grow(
            NULL, run->n_te_contexts, sizeof(*run->context_scratch_right));
        run->context_scratch_capacity = run->n_te_contexts;
    }
}

int tv_build_synteny_contexts(TvRun *run)
{
    size_t side_count;
    TvBlockSide *sides;
    int *parents;
    int *ranks;

    if (run == NULL || run->n_contexts != 0 || run->n_te_contexts != 0) {
        tv_print_error("synteny contexts can be built only once");
        return -1;
    }
    if (run->n_synteny_blocks == 0) {
        return 0;
    }
    if (run->nodes == NULL || run->n_nodes == 0) {
        tv_print_error("TE nodes must exist before building synteny contexts");
        return -1;
    }
    if (run->n_synteny_blocks > SIZE_MAX / 2
        || run->n_synteny_blocks > (size_t)INT_MAX / 2) {
        tv_print_error("too many synteny blocks for context indexing");
        return -1;
    }
    side_count = run->n_synteny_blocks * 2;
    sides = tv_grow(NULL, side_count, sizeof(*sides));
    parents = tv_grow(NULL, side_count, sizeof(*parents));
    ranks = calloc(side_count, sizeof(*ranks));
    if (ranks == NULL) {
        tv_print_error("out of memory");
        exit(EXIT_FAILURE);
    }
    for (size_t block_index = 0; block_index < run->n_synteny_blocks;
         block_index++) {
        const TvSyntenyBlock *block = &run->synteny_blocks[block_index];
        size_t a = block_index * 2;
        size_t b = a + 1;

        sides[a] = (TvBlockSide){
            .block_index = (int)block_index,
            .side_b = false,
            .genome = block->genome_a,
            .contig_index = block->contig_a_index,
            .start = block->start_a,
            .end = block->end_a
        };
        sides[b] = (TvBlockSide){
            .block_index = (int)block_index,
            .side_b = true,
            .genome = block->genome_b,
            .contig_index = block->contig_b_index,
            .start = block->start_b,
            .end = block->end_b
        };
    }
    qsort(sides, side_count, sizeof(*sides), compare_block_side);
    for (size_t side = 0; side < side_count; side++) {
        parents[side] = (int)side;
    }
    for (size_t left = 0; left < side_count; left++) {
        for (size_t right = left + 1; right < side_count; right++) {
            if (sides[right].genome != sides[left].genome
                || sides[right].contig_index != sides[left].contig_index) {
                break;
            }
            if (sides[right].start >= sides[left].end) {
                break;
            }
            if (side_mergeable(run, &sides[left], &sides[right])) {
                parent_union(parents, ranks, (int)left, (int)right);
            }
        }
    }
    if (build_contexts_from_sides(run, sides, side_count, parents) != 0) {
        free(ranks);
        free(parents);
        free(sides);
        return -1;
    }
    assign_homology_groups_and_copies(run);
    if (canonicalize_contexts(run) != 0) {
        free(ranks);
        free(parents);
        free(sides);
        return -1;
    }
    assign_te_contexts(run);
    free(ranks);
    free(parents);
    free(sides);
    return 0;
}

static bool assignment_is_strong(const TvRun *run,
                                 const TvTEContext *assignment)
{
    return assignment->assignment == TV_CONTEXT_BRACKETED
        && assignment->context_index >= 0
        && (size_t)assignment->context_index < run->n_contexts
        && strcmp(run->contexts[assignment->context_index].status, "PASS") == 0;
}

static bool assignment_is_ambiguous(const TvRun *run,
                                    const TvTEContext *assignment)
{
    return assignment->assignment == TV_CONTEXT_AMBIGUOUS
        || (assignment->assignment == TV_CONTEXT_BRACKETED
            && !assignment_is_strong(run, assignment));
}

TvContextRelation tv_candidate_context_relation(
    const TvRun *run, int source_te, int target_te,
    uint64_t *shared_homology_group_id)
{
    bool source_known = false;
    bool target_known = false;
    uint64_t shared_group = 0;
    size_t source_begin;
    size_t source_end;
    size_t target_begin;
    size_t target_end;

    if (shared_homology_group_id != NULL) {
        *shared_homology_group_id = 0;
    }
    if (run == NULL || source_te < 0 || target_te < 0
        || (size_t)source_te >= run->n_nodes
        || (size_t)target_te >= run->n_nodes
        || run->te_context_offsets == NULL) {
        return TV_CONTEXT_RELATION_UNKNOWN;
    }
    source_begin = run->te_context_offsets[source_te];
    source_end = run->te_context_offsets[source_te + 1];
    target_begin = run->te_context_offsets[target_te];
    target_end = run->te_context_offsets[target_te + 1];
    for (size_t index = source_begin; index < source_end; index++) {
        if (assignment_is_ambiguous(run, &run->te_contexts[index])) {
            return TV_CONTEXT_RELATION_AMBIGUOUS;
        }
    }
    for (size_t index = target_begin; index < target_end; index++) {
        if (assignment_is_ambiguous(run, &run->te_contexts[index])) {
            return TV_CONTEXT_RELATION_AMBIGUOUS;
        }
    }
    for (size_t left = source_begin; left < source_end; left++) {
        const TvTEContext *source = &run->te_contexts[left];

        if (!assignment_is_strong(run, source)) {
            continue;
        }
        source_known = true;
        for (size_t right = target_begin; right < target_end; right++) {
            const TvTEContext *target = &run->te_contexts[right];
            uint64_t source_group;
            uint64_t target_group;

            if (!assignment_is_strong(run, target)) {
                continue;
            }
            target_known = true;
            source_group = run->contexts[source->context_index]
                .homology_group_id;
            target_group = run->contexts[target->context_index]
                .homology_group_id;
            if (source_group == target_group) {
                if (shared_group == 0 || source_group < shared_group) {
                    shared_group = source_group;
                }
            }
        }
    }
    if (!source_known || !target_known) {
        return TV_CONTEXT_RELATION_UNKNOWN;
    }
    if (shared_group != 0) {
        if (shared_homology_group_id != NULL) {
            *shared_homology_group_id = shared_group;
        }
        return TV_CONTEXT_RELATION_SUPPORTED;
    }
    return TV_CONTEXT_RELATION_CONFLICT;
}

static bool node_has_strong_context(const TvRun *run, int node,
                                    int *contexts, size_t *count,
                                    size_t capacity)
{
    bool ambiguous = false;

    size_t begin = run->te_context_offsets == NULL
        ? 0 : run->te_context_offsets[node];
    size_t end = run->te_context_offsets == NULL
        ? 0 : run->te_context_offsets[node + 1];

    for (size_t index = begin; index < end; index++) {
        const TvTEContext *assignment = &run->te_contexts[index];
        if (assignment_is_ambiguous(run, assignment)) {
            ambiguous = true;
        } else if (assignment_is_strong(run, assignment)
                   && *count < capacity) {
            contexts[(*count)++] = assignment->context_index;
        }
    }
    return !ambiguous;
}

const char *tv_context_merge_reason(const TvRun *run, int *parents,
                                    int left_root, int right_root)
{
    int *left_contexts;
    int *right_contexts;
    size_t left_count = 0;
    size_t right_count = 0;
    size_t assignment_capacity;
    bool left_usable = true;
    bool right_usable = true;

    if (run == NULL || parents == NULL || run->n_contexts == 0
        || run->context_scratch_capacity == 0) {
        return NULL;
    }
    assignment_capacity = run->context_scratch_capacity;
    left_contexts = run->context_scratch_left;
    right_contexts = run->context_scratch_right;
    for (size_t node = 0; node < run->n_nodes; node++) {
        int component = parent_root(parents, (int)node);

        if (component == left_root) {
            left_usable = node_has_strong_context(
                run, (int)node, left_contexts, &left_count,
                assignment_capacity) && left_usable;
        } else if (component == right_root) {
            right_usable = node_has_strong_context(
                run, (int)node, right_contexts, &right_count,
                assignment_capacity) && right_usable;
        }
    }
    if (!left_usable || !right_usable || left_count == 0 || right_count == 0) {
        return NULL;
    }
    for (size_t left = 0; left < left_count; left++) {
        for (size_t right = 0; right < right_count; right++) {
            if (left_contexts[left] == right_contexts[right]) {
                return "COPY_CONTEXT_CONFLICT";
            }
        }
    }
    for (size_t left = 0; left < left_count; left++) {
        for (size_t right = 0; right < right_count; right++) {
            if (run->contexts[left_contexts[left]].homology_group_id
                != run->contexts[right_contexts[right]].homology_group_id) {
                return "HOMOLOGY_GROUP_CONFLICT";
            }
        }
    }
    return NULL;
}
