#include "tevox.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>

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

static char *strip(char *text)
{
    char *end;

    while (*text != '\0' && isspace((unsigned char)*text)) {
        text++;
    }
    end = text + strlen(text);
    while (end > text && isspace((unsigned char)end[-1])) {
        *--end = '\0';
    }
    return text;
}

static int split_tabs(char *text, char **fields, int capacity)
{
    int count = 0;

    while (count < capacity) {
        char *tab;

        fields[count++] = text;
        tab = strchr(text, '\t');
        if (tab == NULL) {
            break;
        }
        *tab = '\0';
        text = tab + 1;
    }
    if (count > 0) {
        fields[count - 1] = strip(fields[count - 1]);
    }
    return count;
}

static int parse_i64(const char *text, int64_t *value)
{
    char *end = NULL;
    long long parsed;

    errno = 0;
    parsed = strtoll(text, &end, 10);
    if (errno != 0 || end == text || *strip(end) != '\0') {
        return 0;
    }
    *value = (int64_t)parsed;
    return 1;
}

static int parse_int(const char *text, int *value)
{
    int64_t parsed;

    if (!parse_i64(text, &parsed) || parsed < INT_MIN || parsed > INT_MAX) {
        return 0;
    }
    *value = (int)parsed;
    return 1;
}

static char *directory(const char *path)
{
    const char *slash = strrchr(path, '/');
    size_t length;
    char *result;

    if (slash == NULL) {
        return dupstr(".");
    }
    if (slash == path) {
        return dupstr("/");
    }
    length = (size_t)(slash - path);
    result = grow(NULL, length + 1, 1);
    memcpy(result, path, length);
    result[length] = '\0';
    return result;
}

static char *resolve(const char *base, const char *path)
{
    size_t length;
    char *result;

    if (path[0] == '/') {
        return dupstr(path);
    }
    length = strlen(base) + strlen(path) + 2;
    result = grow(NULL, length, 1);
    (void)snprintf(result, length, "%s/%s", base, path);
    return result;
}

static TvContig *find_contig(TvGenome *genome, const char *name)
{
    size_t index;

    for (index = 0; index < genome->n_contigs; index++) {
        if (strcmp(genome->contigs[index].name, name) == 0) {
            return &genome->contigs[index];
        }
    }
    return NULL;
}

static int genome_id(TvRun *run, const char *id)
{
    size_t index;

    for (index = 0; index < run->n_genomes; index++) {
        if (strcmp(run->genomes[index].id, id) == 0) {
            return (int)index;
        }
    }
    return -1;
}

static uint64_t hash_bytes(uint64_t hash, const void *data, size_t length)
{
    const unsigned char *bytes = data;
    size_t index;

    for (index = 0; index < length; index++) {
        hash ^= (uint64_t)bytes[index];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

static uint64_t hash_text(uint64_t hash, const char *text)
{
    return hash_bytes(hash, text, strlen(text));
}

static uint64_t evidence_group_hash(const TvRun *run, int query_genome,
                                    int target_genome, const char *record)
{
    uint64_t hash = UINT64_C(1469598103934665603);

    hash = hash_text(hash, run->genomes[query_genome].id);
    hash = hash_bytes(hash, "\0", 1);
    hash = hash_text(hash, run->genomes[target_genome].id);
    hash = hash_bytes(hash, "\0", 1);
    return hash_text(hash, record);
}

static int fasta_index(TvGenome *genome)
{
    FILE *stream = fopen(genome->fasta_path, "rb");
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length;
    TvContig *current = NULL;
    int previous_bases = 0;
    int line_number = 0;

    if (stream == NULL) {
        tv_print_error("cannot open FASTA '%s': %s", genome->fasta_path,
                       strerror(errno));
        return -1;
    }
    while ((length = getline(&line, &capacity, stream)) >= 0) {
        int bases;

        line_number++;
        if (length > 0 && line[0] == '>') {
            char *name;
            char *space;

            if (current != NULL && previous_bases > current->line_bases) {
                tv_print_error("%s:%d invalid FASTA line width",
                               genome->fasta_path, line_number - 1);
                goto fail;
            }
            name = strip(line + 1);
            space = strpbrk(name, " \t");
            if (space != NULL) {
                *space = '\0';
            }
            if (*name == '\0' || find_contig(genome, name) != NULL) {
                tv_print_error("%s:%d duplicate or empty contig",
                               genome->fasta_path, line_number);
                goto fail;
            }
            if (genome->n_contigs == genome->cap_contigs) {
                genome->cap_contigs = genome->cap_contigs == 0
                    ? 16 : genome->cap_contigs * 2;
                genome->contigs = grow(genome->contigs, genome->cap_contigs,
                                       sizeof(*genome->contigs));
            }
            current = &genome->contigs[genome->n_contigs++];
            memset(current, 0, sizeof(*current));
            current->name = dupstr(name);
            current->seq_offset = (int64_t)ftello(stream);
            previous_bases = 0;
            continue;
        }
        if (current == NULL) {
            if (*strip(line) != '\0') {
                tv_print_error("%s:%d sequence before header",
                               genome->fasta_path, line_number);
                goto fail;
            }
            continue;
        }
        bases = (int)length;
        while (bases > 0 && (line[bases - 1] == '\n' || line[bases - 1] == '\r')) {
            bases--;
        }
        if (bases == 0) {
            continue;
        }
        if (current->line_bases == 0) {
            current->line_bases = bases;
            current->line_bytes = (int)length;
        } else if (previous_bases != current->line_bases) {
            tv_print_error("%s:%d variable-width nonterminal FASTA",
                           genome->fasta_path, line_number - 1);
            goto fail;
        } else if (bases == current->line_bases
                   && (int)length != current->line_bytes) {
            tv_print_error("%s:%d inconsistent newlines", genome->fasta_path,
                           line_number);
            goto fail;
        }
        for (int index = 0; index < bases; index++) {
            if (!isalpha((unsigned char)line[index])
                && strchr("*-.", line[index]) == NULL) {
                tv_print_error("%s:%d invalid FASTA character",
                               genome->fasta_path, line_number);
                goto fail;
            }
        }
        current->length += bases;
        previous_bases = bases;
    }
    free(line);
    fclose(stream);
    if (genome->n_contigs == 0) {
        tv_print_error("FASTA '%s' is empty", genome->fasta_path);
        return -1;
    }
    return 0;

fail:
    free(line);
    fclose(stream);
    return -1;
}

static char *attribute(const char *text, const char *key)
{
    size_t key_length = strlen(key);

    while (*text != '\0') {
        const char *end;
        const char *equals;

        while (*text == ';' || isspace((unsigned char)*text)) {
            text++;
        }
        end = strchr(text, ';');
        if (end == NULL) {
            end = text + strlen(text);
        }
        equals = memchr(text, '=', (size_t)(end - text));
        if (equals != NULL && (size_t)(equals - text) == key_length
            && strncasecmp(text, key, key_length) == 0) {
            size_t value_length = (size_t)(end - equals - 1);
            char *value = grow(NULL, value_length + 1, 1);

            memcpy(value, equals + 1, value_length);
            value[value_length] = '\0';
            return value;
        }
        text = *end == '\0' ? end : end + 1;
    }
    return NULL;
}

static int push_te(TvGenome *genome, TvTE *te, int line_number)
{
    TvContig *sequence = find_contig(genome, te->contig);

    if (sequence == NULL || te->start < 0 || te->end <= te->start
        || te->end > sequence->length) {
        tv_print_error("%s:%d annotation interval outside FASTA",
                       genome->te_path, line_number);
        return -1;
    }
    if (genome->n_tes == genome->cap_tes) {
        genome->cap_tes = genome->cap_tes == 0 ? 64 : genome->cap_tes * 2;
        genome->tes = grow(genome->tes, genome->cap_tes,
                           sizeof(*genome->tes));
    }
    genome->tes[genome->n_tes++] = *te;
    return 0;
}

static int compare_te(const void *left, const void *right)
{
    const TvTE *a = left;
    const TvTE *b = right;
    int comparison = strcmp(a->contig, b->contig);

    if (comparison != 0) {
        return comparison;
    }
    if (a->start != b->start) {
        return a->start < b->start ? -1 : 1;
    }
    if (a->end != b->end) {
        return a->end < b->end ? -1 : 1;
    }
    return strcmp(a->id, b->id);
}

static int annotation(TvGenome *genome, int genome_index)
{
    FILE *stream = fopen(genome->te_path, "r");
    int is_gff = strstr(genome->te_path, ".gff") != NULL
        || strstr(genome->te_path, ".GFF") != NULL;
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length;
    int line_number = 0;

    if (stream == NULL) {
        tv_print_error("cannot open annotation '%s': %s", genome->te_path,
                       strerror(errno));
        return -1;
    }
    while ((length = getline(&line, &capacity, stream)) >= 0) {
        char *row;
        char *fields[16] = {0};
        int count;
        TvTE te = {.genome = genome_index, .strand = '.'};
        int64_t start;
        int64_t end;

        (void)length;
        line_number++;
        row = strip(line);
        if (*row == '\0' || *row == '#') {
            continue;
        }
        count = split_tabs(row, fields, 16);
        if ((is_gff && count < 9) || (!is_gff && count < 3)) {
            tv_print_error("%s:%d malformed annotation", genome->te_path,
                           line_number);
            goto fail;
        }
        te.contig = dupstr(fields[0]);
        if (is_gff) {
            if (!parse_i64(fields[3], &start) || !parse_i64(fields[4], &end)) {
                goto record_bad;
            }
            te.start = start - 1;
            te.end = end;
            te.strand = *fields[6] == '\0' ? '.' : *fields[6];
            te.type = dupstr(fields[2]);
            te.id = attribute(fields[8], "ID");
            if (te.id == NULL) {
                te.id = attribute(fields[8], "Name");
            }
            te.family = attribute(fields[8], "family");
            if (te.family == NULL) {
                te.family = attribute(fields[8], "Family");
            }
            if (te.family == NULL) {
                te.family = attribute(fields[8], "classification");
            }
        } else {
            if (!parse_i64(fields[1], &start) || !parse_i64(fields[2], &end)) {
                goto record_bad;
            }
            te.start = start;
            te.end = end;
            te.id = count > 3 && *fields[3] != '\0' ? dupstr(fields[3]) : NULL;
            te.strand = count > 5 && *fields[5] != '\0' ? *fields[5] : '.';
            te.family = count > 6 && *fields[6] != '\0'
                ? dupstr(fields[6]) : NULL;
            te.type = count > 7 && *fields[7] != '\0'
                ? dupstr(fields[7]) : NULL;
        }
        if (te.id == NULL) {
            char id[192];

            (void)snprintf(id, sizeof(id), "%s:%s:%lld-%lld", genome->id,
                           te.contig, (long long)te.start, (long long)te.end);
            te.id = dupstr(id);
        }
        if (te.type == NULL) {
            te.type = dupstr("transposable_element");
        }
        if (te.family == NULL) {
            te.family = dupstr(TEVOX_UNKNOWN);
        }
        if (push_te(genome, &te, line_number) != 0) {
            goto record_bad;
        }
        continue;

record_bad:
        free(te.id);
        free(te.contig);
        free(te.type);
        free(te.family);
        tv_print_error("%s:%d invalid annotation record", genome->te_path,
                       line_number);
        goto fail;
    }
    free(line);
    fclose(stream);
    if (genome->n_tes == 0) {
        tv_print_error("annotation '%s' contains no records", genome->te_path);
        return -1;
    }
    qsort(genome->tes, genome->n_tes, sizeof(*genome->tes), compare_te);
    for (size_t index = 1; index < genome->n_tes; index++) {
        if (strcmp(genome->tes[index - 1].id, genome->tes[index].id) == 0) {
            tv_print_error("annotation '%s' contains duplicate TE ID '%s'",
                           genome->te_path, genome->tes[index].id);
            return -1;
        }
    }
    return 0;

fail:
    free(line);
    fclose(stream);
    return -1;
}

void tv_run_init(TvRun *run)
{
    memset(run, 0, sizeof(*run));
    run->cfg.flank = 100;
    run->cfg.candidate_window = 100;
    run->cfg.min_mapq = 20;
    run->cfg.min_flank_fraction = 0.60;
    run->cfg.min_edge_score = 45.0;
    run->cfg.max_n_fraction = 0.25;
    run->cfg.min_reciprocal_overlap = 0.50;
    run->cfg.near_best_delta = 5.0;
}

void tv_run_free(TvRun *run)
{
    size_t genome_index;
    size_t index;

    for (genome_index = 0; genome_index < run->n_genomes; genome_index++) {
        TvGenome *genome = &run->genomes[genome_index];

        free(genome->id);
        free(genome->fasta_path);
        free(genome->te_path);
        for (index = 0; index < genome->n_contigs; index++) {
            free(genome->contigs[index].name);
        }
        free(genome->contigs);
        for (index = 0; index < genome->n_tes; index++) {
            free(genome->tes[index].id);
            free(genome->tes[index].contig);
            free(genome->tes[index].type);
            free(genome->tes[index].family);
        }
        free(genome->tes);
    }
    for (index = 0; index < run->n_pafs; index++) {
        free(run->pafs[index].qname);
        free(run->pafs[index].tname);
        free(run->pafs[index].ops);
        free(run->pafs[index].identity_ops);
        free(run->pafs[index].source_path);
    }
    for (index = 0; index < run->n_projections; index++) {
        free(run->projections[index].contig);
    }
    free(run->genomes);
    free(run->nodes);
    free(run->pafs);
    free(run->projections);
    free(run->candidates);
    free(run->decisions);
    free(run->edges);
    free(run->components);
    memset(run, 0, sizeof(*run));
}

int tv_add_genome(TvRun *run, const char *id, const char *fasta,
                  const char *annotation_path, int copies)
{
    int index;
    TvGenome *genome;

    if (id == NULL || *id == '\0' || genome_id(run, id) >= 0 || copies < 1) {
        tv_print_error("genome IDs must be unique and copy quota positive");
        return -1;
    }
    if (run->n_genomes == run->cap_genomes) {
        run->cap_genomes = run->cap_genomes == 0 ? 8 : run->cap_genomes * 2;
        run->genomes = grow(run->genomes, run->cap_genomes,
                            sizeof(*run->genomes));
    }
    index = (int)run->n_genomes;
    genome = &run->genomes[run->n_genomes++];
    memset(genome, 0, sizeof(*genome));
    genome->id = dupstr(id);
    genome->fasta_path = dupstr(fasta);
    genome->te_path = dupstr(annotation_path);
    genome->max_locus_copies = copies;
    if (fasta_index(genome) != 0 || annotation(genome, index) != 0) {
        return -1;
    }
    return index;
}

static int compare_genome(const void *left, const void *right)
{
    const TvGenome *a = left;
    const TvGenome *b = right;

    return strcmp(a->id, b->id);
}

int tv_load_manifest(TvRun *run, const char *path)
{
    FILE *stream = fopen(path, "r");
    char *base;
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length;
    int line_number = 0;

    if (stream == NULL) {
        tv_print_error("cannot open manifest '%s'", path);
        return -1;
    }
    base = directory(path);
    while ((length = getline(&line, &capacity, stream)) >= 0) {
        char *row;
        char *fields[5] = {0};
        int count;
        int copies = 1;
        char *fasta;
        char *annotations;
        int added;

        (void)length;
        line_number++;
        row = strip(line);
        if (*row == '\0' || *row == '#') {
            continue;
        }
        count = split_tabs(row, fields, 5);
        if (strcmp(fields[0], "genome_id") == 0) {
            continue;
        }
        if (count < 3 || count > 4) {
            tv_print_error("%s:%d expected 3 or 4 columns", path, line_number);
            goto fail;
        }
        if (count == 4 && (!parse_int(fields[3], &copies) || copies < 1)) {
            tv_print_error("%s:%d invalid quota", path, line_number);
            goto fail;
        }
        fasta = resolve(base, fields[1]);
        annotations = resolve(base, fields[2]);
        added = tv_add_genome(run, fields[0], fasta, annotations, copies);
        free(fasta);
        free(annotations);
        if (added < 0) {
            goto fail;
        }
    }
    free(base);
    free(line);
    fclose(stream);
    if (run->n_genomes < 2) {
        tv_print_error("manifest needs at least two genomes");
        return -1;
    }
    qsort(run->genomes, run->n_genomes, sizeof(*run->genomes), compare_genome);
    for (size_t genome_index = 0; genome_index < run->n_genomes; genome_index++) {
        for (size_t te_index = 0; te_index < run->genomes[genome_index].n_tes;
             te_index++) {
            run->genomes[genome_index].tes[te_index].genome = (int)genome_index;
        }
    }
    return 0;

fail:
    free(base);
    free(line);
    fclose(stream);
    return -1;
}

static int parse_cigar(const char *text, TvCigarOp **output,
                       size_t *output_count, int64_t query_span,
                       int64_t target_span)
{
    TvCigarOp *operations = NULL;
    size_t count = 0;
    size_t capacity = 0;
    int64_t query_consumed = 0;
    int64_t target_consumed = 0;

    while (*text != '\0') {
        char *end = NULL;
        long long length;
        char code;

        errno = 0;
        length = strtoll(text, &end, 10);
        if (errno != 0 || end == text || length <= 0 || *end == '\0') {
            goto invalid;
        }
        code = *end++;
        if (strchr("M=XID", code) == NULL) {
            goto invalid;
        }
        if (count == capacity) {
            capacity = capacity == 0 ? 16 : capacity * 2;
            operations = grow(operations, capacity, sizeof(*operations));
        }
        operations[count++] = (TvCigarOp){code, (int64_t)length};
        if (strchr("M=XI", code) != NULL) {
            query_consumed += length;
        }
        if (strchr("M=XD", code) != NULL) {
            target_consumed += length;
        }
        text = end;
    }
    if (count == 0 || query_consumed != query_span
        || target_consumed != target_span) {
        goto invalid;
    }
    *output = operations;
    *output_count = count;
    return 0;

invalid:
    free(operations);
    return -1;
}

static int push_identity_op(TvCigarOp **operations, size_t *count,
                            size_t *capacity, char code, int64_t length)
{
    if (length <= 0) {
        return -1;
    }
    if (*count > 0 && (*operations)[*count - 1].code == code) {
        if ((*operations)[*count - 1].length > INT64_MAX - length) {
            return -1;
        }
        (*operations)[*count - 1].length += length;
        return 0;
    }
    if (*count == *capacity) {
        *capacity = *capacity == 0 ? 16 : *capacity * 2;
        *operations = grow(*operations, *capacity, sizeof(**operations));
    }
    (*operations)[(*count)++] = (TvCigarOp){code, length};
    return 0;
}

static int parse_cs(const char *text, TvCigarOp **output,
                    size_t *output_count, int64_t query_span,
                    int64_t target_span)
{
    TvCigarOp *operations = NULL;
    size_t count = 0;
    size_t capacity = 0;
    int64_t query_consumed = 0;
    int64_t target_consumed = 0;

    while (*text != '\0') {
        char operation = *text++;
        int64_t length = 0;
        char code;

        if (operation == ':') {
            char *end = NULL;
            long long parsed;

            errno = 0;
            parsed = strtoll(text, &end, 10);
            if (errno != 0 || end == text || parsed <= 0) {
                goto invalid;
            }
            length = (int64_t)parsed;
            text = end;
            code = '=';
        } else if (operation == '=') {
            const char *start = text;

            while (isalpha((unsigned char)*text)) {
                text++;
            }
            length = (int64_t)(text - start);
            code = '=';
        } else if (operation == '*') {
            if (!isalpha((unsigned char)text[0])
                || !isalpha((unsigned char)text[1])) {
                goto invalid;
            }
            text += 2;
            length = 1;
            code = 'X';
        } else if (operation == '+' || operation == '-') {
            const char *start = text;

            while (isalpha((unsigned char)*text)) {
                text++;
            }
            length = (int64_t)(text - start);
            code = operation == '+' ? 'I' : 'D';
        } else {
            /* Splice (~) and unknown cs operations are outside assembly PAF. */
            goto invalid;
        }
        if (push_identity_op(&operations, &count, &capacity, code, length) != 0) {
            goto invalid;
        }
        if (strchr("=XI", code) != NULL) {
            query_consumed += length;
        }
        if (strchr("=XD", code) != NULL) {
            target_consumed += length;
        }
    }
    if (count == 0 || query_consumed != query_span
        || target_consumed != target_span) {
        goto invalid;
    }
    *output = operations;
    *output_count = count;
    return 0;

invalid:
    free(operations);
    return -1;
}

static bool cigar_has_exact_identity(const TvCigarOp *operations, size_t count)
{
    bool has_exact = false;

    for (size_t index = 0; index < count; index++) {
        if (operations[index].code == 'M') {
            return false;
        }
        if (operations[index].code == '=' || operations[index].code == 'X') {
            has_exact = true;
        }
    }
    return has_exact;
}

static bool cigar_cs_compatible(const TvCigarOp *cigar, size_t cigar_count,
                                const TvCigarOp *identity,
                                size_t identity_count)
{
    size_t cigar_index = 0;
    size_t identity_index = 0;
    int64_t cigar_remaining = cigar_count > 0 ? cigar[0].length : 0;
    int64_t identity_remaining = identity_count > 0 ? identity[0].length : 0;

    while (cigar_index < cigar_count && identity_index < identity_count) {
        char cigar_code = cigar[cigar_index].code;
        char identity_code = identity[identity_index].code;
        bool compatible;
        int64_t consumed;

        if (cigar_code == 'M') {
            compatible = identity_code == '=' || identity_code == 'X';
        } else {
            compatible = cigar_code == identity_code;
        }
        if (!compatible) {
            return false;
        }
        consumed = cigar_remaining < identity_remaining
            ? cigar_remaining : identity_remaining;
        cigar_remaining -= consumed;
        identity_remaining -= consumed;
        if (cigar_remaining == 0) {
            cigar_index++;
            if (cigar_index < cigar_count) {
                cigar_remaining = cigar[cigar_index].length;
            }
        }
        if (identity_remaining == 0) {
            identity_index++;
            if (identity_index < identity_count) {
                identity_remaining = identity[identity_index].length;
            }
        }
    }
    return cigar_index == cigar_count && identity_index == identity_count;
}

static void push_paf(TvRun *run, TvPaf *paf)
{
    if (run->n_pafs == run->cap_pafs) {
        run->cap_pafs = run->cap_pafs == 0 ? 128 : run->cap_pafs * 2;
        run->pafs = grow(run->pafs, run->cap_pafs, sizeof(*run->pafs));
    }
    run->pafs[run->n_pafs++] = *paf;
}

static void reverse_operations(const TvCigarOp *source, size_t count,
                               char strand, TvCigarOp *target)
{
    for (size_t index = 0; index < count; index++) {
        size_t source_index = strand == '-' ? count - index - 1 : index;

        target[index] = source[source_index];
        if (target[index].code == 'I') {
            target[index].code = 'D';
        } else if (target[index].code == 'D') {
            target[index].code = 'I';
        }
    }
}

static void reverse_paf(TvRun *run, const TvPaf *source)
{
    TvPaf paf = {
        .query_genome = source->target_genome,
        .target_genome = source->query_genome,
        .qname = dupstr(source->tname),
        .qlen = source->tlen,
        .qstart = source->tstart,
        .qend = source->tend,
        .strand = source->strand,
        .tname = dupstr(source->qname),
        .tlen = source->qlen,
        .tstart = source->qstart,
        .tend = source->qend,
        .matches = source->matches,
        .block_len = source->block_len,
        .mapq = source->mapq,
        .mapq_observed = source->mapq_observed,
        .n_ops = source->n_ops,
        .n_identity_ops = source->n_identity_ops,
        .identity_method = source->identity_method,
        .origin = TV_EVIDENCE_DERIVED_REVERSE,
        .source_path = dupstr(source->source_path),
        .source_line = source->source_line,
        .evidence_group_id = source->evidence_group_id
    };

    paf.ops = grow(NULL, paf.n_ops, sizeof(*paf.ops));
    reverse_operations(source->ops, source->n_ops, source->strand, paf.ops);
    if (paf.n_identity_ops > 0) {
        paf.identity_ops = grow(NULL, paf.n_identity_ops,
                                sizeof(*paf.identity_ops));
        reverse_operations(source->identity_ops, source->n_identity_ops,
                           source->strand, paf.identity_ops);
    }
    push_paf(run, &paf);
}

static bool duplicate_native_group(const TvRun *run, uint64_t group_id)
{
    for (size_t index = 0; index < run->n_pafs; index++) {
        if (run->pafs[index].origin == TV_EVIDENCE_NATIVE
            && run->pafs[index].evidence_group_id == group_id) {
            return true;
        }
    }
    return false;
}

int tv_add_paf_file(TvRun *run, int query_genome, int target_genome,
                    const char *path)
{
    FILE *stream;
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length;
    int line_number = 0;
    int records = 0;

    if (query_genome < 0 || target_genome < 0
        || query_genome == target_genome
        || query_genome >= (int)run->n_genomes
        || target_genome >= (int)run->n_genomes) {
        return -1;
    }
    stream = fopen(path, "r");
    if (stream == NULL) {
        tv_print_error("cannot open PAF '%s'", path);
        return -1;
    }
    while ((length = getline(&line, &capacity, stream)) >= 0) {
        char *row;
        char *record_text;
        char *fields[128] = {0};
        int count;
        const char *cg = NULL;
        const char *cs = NULL;
        TvPaf paf = {
            .query_genome = query_genome,
            .target_genome = target_genome,
            .origin = TV_EVIDENCE_NATIVE,
            .identity_method = TV_IDENTITY_MISSING
        };
        int parsed_mapq;
        TvContig *query_contig;
        TvContig *target_contig;

        (void)length;
        line_number++;
        row = strip(line);
        if (*row == '\0' || *row == '#') {
            continue;
        }
        record_text = dupstr(row);
        count = split_tabs(row, fields, 128);
        if (count < 12) {
            free(record_text);
            tv_print_error("%s:%d PAF requires 12 fields", path, line_number);
            goto fail;
        }
        paf.qname = dupstr(fields[0]);
        paf.strand = *fields[4];
        paf.tname = dupstr(fields[5]);
        if (!parse_i64(fields[1], &paf.qlen)
            || !parse_i64(fields[2], &paf.qstart)
            || !parse_i64(fields[3], &paf.qend)
            || !parse_i64(fields[6], &paf.tlen)
            || !parse_i64(fields[7], &paf.tstart)
            || !parse_i64(fields[8], &paf.tend)
            || !parse_i64(fields[9], &paf.matches)
            || !parse_i64(fields[10], &paf.block_len)
            || !parse_int(fields[11], &parsed_mapq)
            || (paf.strand != '+' && paf.strand != '-')
            || paf.qstart < 0 || paf.qend <= paf.qstart
            || paf.tstart < 0 || paf.tend <= paf.tstart
            || paf.qend > paf.qlen || paf.tend > paf.tlen
            || paf.matches < 0 || paf.block_len <= 0
            || paf.matches > paf.block_len
            || parsed_mapq < 0 || parsed_mapq > 255) {
            free(record_text);
            goto record_bad;
        }
        query_contig = find_contig(&run->genomes[query_genome], paf.qname);
        target_contig = find_contig(&run->genomes[target_genome], paf.tname);
        if (query_contig == NULL || target_contig == NULL
            || query_contig->length != paf.qlen
            || target_contig->length != paf.tlen) {
            free(record_text);
            goto record_bad;
        }
        for (int index = 12; index < count; index++) {
            if (strncmp(fields[index], "cg:Z:", 5) == 0) {
                cg = fields[index] + 5;
            } else if (strncmp(fields[index], "cs:Z:", 5) == 0) {
                cs = fields[index] + 5;
            }
        }
        if (cg == NULL
            || parse_cigar(cg, &paf.ops, &paf.n_ops,
                           paf.qend - paf.qstart, paf.tend - paf.tstart) != 0) {
            free(record_text);
            goto record_bad;
        }
        if (cs != NULL) {
            if (parse_cs(cs, &paf.identity_ops, &paf.n_identity_ops,
                         paf.qend - paf.qstart,
                         paf.tend - paf.tstart) != 0
                || !cigar_cs_compatible(paf.ops, paf.n_ops,
                                        paf.identity_ops,
                                        paf.n_identity_ops)) {
                free(record_text);
                goto record_bad;
            }
            paf.identity_method = TV_IDENTITY_CS;
        } else if (cigar_has_exact_identity(paf.ops, paf.n_ops)) {
            paf.identity_ops = grow(NULL, paf.n_ops, sizeof(*paf.identity_ops));
            memcpy(paf.identity_ops, paf.ops, paf.n_ops * sizeof(*paf.ops));
            paf.n_identity_ops = paf.n_ops;
            paf.identity_method = TV_IDENTITY_EQX;
        }
        paf.mapq_observed = parsed_mapq != 255;
        paf.mapq = paf.mapq_observed ? parsed_mapq : -1;
        paf.source_path = dupstr(path);
        paf.source_line = line_number;
        paf.evidence_group_id = evidence_group_hash(
            run, query_genome, target_genome, record_text);
        free(record_text);
        if (duplicate_native_group(run, paf.evidence_group_id)) {
            free(paf.qname);
            free(paf.tname);
            free(paf.ops);
            free(paf.identity_ops);
            free(paf.source_path);
            if (run->cfg.verbose) {
                fprintf(stderr, "tevox: ignored duplicate PAF evidence at %s:%d\n",
                        path, line_number);
            }
            continue;
        }
        push_paf(run, &paf);
        reverse_paf(run, &paf);
        records++;
        continue;

record_bad:
        free(paf.qname);
        free(paf.tname);
        free(paf.ops);
        free(paf.identity_ops);
        free(paf.source_path);
        tv_print_error("%s:%d invalid PAF, cg:Z/cs:Z, or declared direction",
                       path, line_number);
        goto fail;
    }
    free(line);
    fclose(stream);
    if (records == 0) {
        tv_print_error("PAF '%s' contains no new records", path);
        return -1;
    }
    return 0;

fail:
    free(line);
    fclose(stream);
    return -1;
}

int tv_load_alignments(TvRun *run, const char *path)
{
    FILE *stream = fopen(path, "r");
    char *base;
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length;
    int line_number = 0;
    int pairs = 0;

    if (stream == NULL) {
        tv_print_error("cannot open alignments table '%s'", path);
        return -1;
    }
    base = directory(path);
    while ((length = getline(&line, &capacity, stream)) >= 0) {
        char *row;
        char *fields[4] = {0};
        int count;
        int query;
        int target;
        char *paf_path;
        int loaded;

        (void)length;
        line_number++;
        row = strip(line);
        if (*row == '\0' || *row == '#') {
            continue;
        }
        count = split_tabs(row, fields, 4);
        if (strcmp(fields[0], "query_id") == 0) {
            continue;
        }
        if (count != 3) {
            tv_print_error("%s:%d expected three columns", path, line_number);
            goto fail;
        }
        query = genome_id(run, fields[0]);
        target = genome_id(run, fields[1]);
        if (query < 0 || target < 0 || query == target) {
            tv_print_error("%s:%d unknown or identical query/target genome",
                           path, line_number);
            goto fail;
        }
        paf_path = resolve(base, fields[2]);
        loaded = tv_add_paf_file(run, query, target, paf_path);
        free(paf_path);
        if (loaded != 0) {
            goto fail;
        }
        pairs++;
    }
    free(base);
    free(line);
    fclose(stream);
    if (pairs == 0) {
        tv_print_error("alignment table is empty");
        return -1;
    }
    return 0;

fail:
    free(base);
    free(line);
    fclose(stream);
    return -1;
}
