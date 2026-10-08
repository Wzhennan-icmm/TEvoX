#include "tevox.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

typedef enum {
    DELTA_EXPECT_FILE_HEADER,
    DELTA_EXPECT_DATA_TYPE,
    DELTA_EXPECT_RECORD,
    DELTA_READ_VALUES
} DeltaParserState;

typedef struct {
    char *data;
    size_t length;
    size_t capacity;
} TextBuffer;

static int split_whitespace(char *text, char **fields, size_t capacity)
{
    size_t count = 0;

    while (*text != '\0') {
        while (isspace((unsigned char)*text)) {
            text++;
        }
        if (*text == '\0') {
            break;
        }
        if (count == capacity) {
            return -1;
        }
        fields[count++] = text;
        while (*text != '\0' && !isspace((unsigned char)*text)) {
            text++;
        }
        if (*text != '\0') {
            *text++ = '\0';
        }
    }
    return (int)count;
}

static int append_operation(TvCigarOp **operations, size_t *count,
                            size_t *capacity, char code, int64_t length)
{
    if (length < 0 || strchr("MID", code) == NULL) {
        return -1;
    }
    if (length == 0) {
        return 0;
    }
    if (*count > 0 && (*operations)[*count - 1].code == code) {
        if ((*operations)[*count - 1].length > INT64_MAX - length) {
            return -1;
        }
        (*operations)[*count - 1].length += length;
        return 0;
    }
    if (*count == *capacity) {
        if (*capacity > SIZE_MAX / 2) {
            return -1;
        }
        *capacity = *capacity == 0 ? 16 : *capacity * 2;
        *operations = tv_grow(*operations, *capacity, sizeof(**operations));
    }
    (*operations)[(*count)++] = (TvCigarOp){code, length};
    return 0;
}

static void reverse_operation_order(TvCigarOp *operations, size_t count)
{
    for (size_t left = 0, right = count == 0 ? 0 : count - 1;
         left < right; left++, right--) {
        TvCigarOp temporary = operations[left];

        operations[left] = operations[right];
        operations[right] = temporary;
    }
}

static void free_alignment(TvPaf *alignment)
{
    free(alignment->qname);
    free(alignment->tname);
    free(alignment->ops);
    free(alignment->identity_ops);
    free(alignment->source_path);
    memset(alignment, 0, sizeof(*alignment));
}

static int text_buffer_appendf(TextBuffer *buffer, const char *format, ...)
{
    va_list arguments;
    va_list copy;
    int needed;
    size_t required;

    va_start(arguments, format);
    va_copy(copy, arguments);
    needed = vsnprintf(NULL, 0, format, copy);
    va_end(copy);
    if (needed < 0 || (size_t)needed > SIZE_MAX - buffer->length - 1) {
        va_end(arguments);
        return -1;
    }
    required = buffer->length + (size_t)needed + 1;
    if (required > buffer->capacity) {
        size_t new_capacity = buffer->capacity == 0 ? 256 : buffer->capacity;

        while (new_capacity < required) {
            if (new_capacity > SIZE_MAX / 2) {
                new_capacity = required;
                break;
            }
            new_capacity *= 2;
        }
        buffer->data = tv_grow(buffer->data, new_capacity, 1);
        buffer->capacity = new_capacity;
    }
    (void)vsnprintf(buffer->data + buffer->length,
                    buffer->capacity - buffer->length, format, arguments);
    va_end(arguments);
    buffer->length += (size_t)needed;
    return 0;
}

static char *canonical_record(const TvPaf *alignment)
{
    TextBuffer buffer = {0};

    if (text_buffer_appendf(
            &buffer,
            "MUMMER_DELTA_V1;Q%zu:%s;%lld;%lld;%lld;S%c;T%zu:%s;"
            "%lld;%lld;%lld;E%lld;%lld;%lld;C",
            strlen(alignment->qname), alignment->qname,
            (long long)alignment->qlen, (long long)alignment->qstart,
            (long long)alignment->qend, alignment->strand,
            strlen(alignment->tname), alignment->tname,
            (long long)alignment->tlen, (long long)alignment->tstart,
            (long long)alignment->tend, (long long)alignment->error_count,
            (long long)alignment->similarity_error_count,
            (long long)alignment->nonalpha_count) != 0) {
        free(buffer.data);
        return NULL;
    }
    for (size_t index = 0; index < alignment->n_ops; index++) {
        if (text_buffer_appendf(&buffer, "%lld%c;",
                                (long long)alignment->ops[index].length,
                                alignment->ops[index].code) != 0) {
            free(buffer.data);
            return NULL;
        }
    }
    return buffer.data;
}

static int invalid_delta(const char *path, int line_number,
                         const char *message)
{
    tv_print_error("%s:%d invalid MUMmer delta: %s", path, line_number,
                   message);
    return -1;
}

static const char *path_basename(const char *path)
{
    const char *slash = strrchr(path, '/');

    return slash == NULL ? path : slash + 1;
}

static int validate_file_header(const TvRun *run, int query_genome,
                                int target_genome, const char *reference_path,
                                const char *query_path, const char *delta_path,
                                int line_number)
{
    const char *expected_reference =
        run->genomes[target_genome].fasta_path;
    const char *expected_query = run->genomes[query_genome].fasta_path;
    const char *expected_reference_base = path_basename(expected_reference);
    const char *expected_query_base = path_basename(expected_query);
    const char *reference_base = path_basename(reference_path);
    const char *query_base = path_basename(query_path);
    char *canonical_expected_reference = realpath(expected_reference, NULL);
    char *canonical_expected_query = realpath(expected_query, NULL);
    char *canonical_reference = realpath(reference_path, NULL);
    char *canonical_query = realpath(query_path, NULL);
    bool reference_comparable = canonical_expected_reference != NULL
        && canonical_reference != NULL;
    bool query_comparable = canonical_expected_query != NULL
        && canonical_query != NULL;
    int result = -1;

    if ((reference_comparable
         && strcmp(canonical_expected_reference, canonical_reference) != 0)
        || (query_comparable
            && strcmp(canonical_expected_query, canonical_query) != 0)) {
        (void)invalid_delta(
            delta_path, line_number,
            "file header direction disagrees with target/reference and query FASTA paths");
        goto done;
    }
    if (reference_comparable && query_comparable) {
        result = 0;
        goto done;
    }

    if (*expected_reference_base == '\0' || *expected_query_base == '\0'
        || *reference_base == '\0' || *query_base == '\0'
        || strcmp(expected_reference_base, expected_query_base) == 0) {
        (void)invalid_delta(
            delta_path, line_number,
            "file header direction is ambiguous because canonical paths are unavailable and expected FASTA basenames are not unique");
        goto done;
    }
    if (strcmp(reference_base, expected_reference_base) != 0
        || strcmp(query_base, expected_query_base) != 0) {
        (void)invalid_delta(
            delta_path, line_number,
            "file header must list target/reference FASTA first and query FASTA second");
        goto done;
    }
    result = 0;

done:
    free(canonical_expected_reference);
    free(canonical_expected_query);
    free(canonical_reference);
    free(canonical_query);
    return result;
}

int tv_add_mummer_delta_file(TvRun *run, int query_genome, int target_genome,
                             const char *path)
{
    FILE *stream = NULL;
    char *line = NULL;
    size_t line_capacity = 0;
    ssize_t line_length;
    int line_number = 0;
    int alignment_ordinal = 0;
    int parsed_records = 0;
    int added_records = 0;
    DeltaParserState state = DELTA_EXPECT_FILE_HEADER;
    bool have_sequence_header = false;
    int records_in_sequence = 0;
    char *reference_name = NULL;
    char *query_name = NULL;
    int64_t reference_length = 0;
    int64_t query_length = 0;
    int reference_contig_index = -1;
    int query_contig_index = -1;
    TvPaf alignment;
    size_t operation_capacity = 0;
    int64_t query_span = 0;
    int64_t target_span = 0;
    int64_t query_consumed = 0;
    int64_t target_consumed = 0;
    int64_t positive_deltas = 0;
    int64_t negative_deltas = 0;
    int64_t delta_count = 0;

    memset(&alignment, 0, sizeof(alignment));
    if (query_genome < 0 || target_genome < 0
        || query_genome == target_genome
        || query_genome >= (int)run->n_genomes
        || target_genome >= (int)run->n_genomes) {
        return -1;
    }
    if (tv_register_input(run, "alignment", path) != 0) {
        return -1;
    }
    stream = fopen(path, "r");
    if (stream == NULL) {
        tv_print_error("cannot open MUMmer delta '%s': %s", path,
                       strerror(errno));
        return -1;
    }

    while ((line_length = getline(&line, &line_capacity, stream)) >= 0) {
        char *row;
        char *fields[8] = {0};
        int field_count;

        (void)line_length;
        if (line_number == INT_MAX) {
            (void)invalid_delta(path, line_number, "too many input lines");
            goto fail;
        }
        line_number++;
        row = tv_strip(line);
        if (*row == '\0') {
            continue;
        }

        if (state == DELTA_EXPECT_FILE_HEADER) {
            field_count = split_whitespace(row, fields, 3);
            if (field_count != 2) {
                (void)invalid_delta(path, line_number,
                                    "expected reference and query paths");
                goto fail;
            }
            if (validate_file_header(run, query_genome, target_genome,
                                     fields[0], fields[1], path,
                                     line_number) != 0) {
                goto fail;
            }
            state = DELTA_EXPECT_DATA_TYPE;
            continue;
        }

        if (state == DELTA_EXPECT_DATA_TYPE) {
            field_count = split_whitespace(row, fields, 2);
            if (field_count != 1 || strcmp(fields[0], "NUCMER") != 0) {
                if (field_count == 1 && strcmp(fields[0], "PROMER") == 0) {
                    (void)invalid_delta(path, line_number,
                                        "PROMER input is not DNA NUCMER data");
                } else {
                    (void)invalid_delta(path, line_number,
                                        "expected NUCMER data type");
                }
                goto fail;
            }
            state = DELTA_EXPECT_RECORD;
            continue;
        }

        if (state == DELTA_EXPECT_RECORD && *row == '>') {
            TvContig *reference_contig;
            TvContig *query_contig;

            if (have_sequence_header && records_in_sequence == 0) {
                (void)invalid_delta(path, line_number,
                                    "sequence header has no alignments");
                goto fail;
            }
            field_count = split_whitespace(tv_strip(row + 1), fields, 5);
            if (field_count != 4
                || *fields[0] == '\0' || *fields[1] == '\0'
                || !tv_parse_i64(fields[2], &reference_length)
                || !tv_parse_i64(fields[3], &query_length)
                || reference_length <= 0 || query_length <= 0) {
                (void)invalid_delta(path, line_number,
                                    "malformed sequence header");
                goto fail;
            }
            reference_contig = tv_find_contig(
                &run->genomes[target_genome], fields[0]);
            query_contig = tv_find_contig(
                &run->genomes[query_genome], fields[1]);
            if (reference_contig == NULL || query_contig == NULL
                || reference_contig->length != reference_length
                || query_contig->length != query_length) {
                (void)invalid_delta(
                    path, line_number,
                    "reference=target/query=query contigs or lengths disagree with FASTA");
                goto fail;
            }
            free(reference_name);
            free(query_name);
            reference_name = tv_dupstr(fields[0]);
            query_name = tv_dupstr(fields[1]);
            reference_contig_index = reference_contig->index;
            query_contig_index = query_contig->index;
            have_sequence_header = true;
            records_in_sequence = 0;
            continue;
        }

        if (state == DELTA_EXPECT_RECORD) {
            int64_t r_start;
            int64_t r_end;
            int64_t q_first;
            int64_t q_last;
            int64_t errors;
            int64_t similarity_errors;
            int64_t nonalpha;

            if (!have_sequence_header) {
                (void)invalid_delta(path, line_number,
                                    "alignment before sequence header");
                goto fail;
            }
            field_count = split_whitespace(row, fields, 8);
            if (field_count != 7
                || !tv_parse_i64(fields[0], &r_start)
                || !tv_parse_i64(fields[1], &r_end)
                || !tv_parse_i64(fields[2], &q_first)
                || !tv_parse_i64(fields[3], &q_last)
                || !tv_parse_i64(fields[4], &errors)
                || !tv_parse_i64(fields[5], &similarity_errors)
                || !tv_parse_i64(fields[6], &nonalpha)
                || r_start <= 0 || r_end < r_start
                || q_first <= 0 || q_last <= 0
                || r_end > reference_length
                || q_first > query_length || q_last > query_length
                || errors < 0 || similarity_errors < 0 || nonalpha < 0) {
                (void)invalid_delta(path, line_number,
                                    "malformed alignment header");
                goto fail;
            }
            if (alignment_ordinal == INT_MAX) {
                (void)invalid_delta(path, line_number,
                                    "too many alignment records");
                goto fail;
            }
            alignment_ordinal++;
            memset(&alignment, 0, sizeof(alignment));
            alignment.provider = TV_ALIGNMENT_MUMMER_DELTA;
            alignment.query_genome = query_genome;
            alignment.target_genome = target_genome;
            alignment.qname = tv_dupstr(query_name);
            alignment.qcontig_index = query_contig_index;
            alignment.qlen = query_length;
            alignment.qstart = (q_first < q_last ? q_first : q_last) - 1;
            alignment.qend = q_first > q_last ? q_first : q_last;
            alignment.strand = q_first <= q_last ? '+' : '-';
            alignment.tname = tv_dupstr(reference_name);
            alignment.tcontig_index = reference_contig_index;
            alignment.tlen = reference_length;
            alignment.tstart = r_start - 1;
            alignment.tend = r_end;
            alignment.mapq = -1;
            alignment.mapq_observed = false;
            alignment.mapq_status = TV_MAPQ_NOT_PROVIDED;
            alignment.aggregate_identity = NAN;
            alignment.aggregate_identity_method =
                TV_AGG_IDENTITY_DELTA_ERRORS;
            alignment.error_count = errors;
            alignment.similarity_error_count = similarity_errors;
            alignment.nonalpha_count = nonalpha;
            alignment.identity_method = TV_IDENTITY_MISSING;
            alignment.origin = TV_EVIDENCE_NATIVE;
            alignment.source_path = tv_dupstr(path);
            alignment.source_line = line_number;
            alignment.source_record = alignment_ordinal;
            operation_capacity = 0;
            query_span = alignment.qend - alignment.qstart;
            target_span = alignment.tend - alignment.tstart;
            query_consumed = 0;
            target_consumed = 0;
            positive_deltas = 0;
            negative_deltas = 0;
            delta_count = 0;
            state = DELTA_READ_VALUES;
            continue;
        }

        if (state == DELTA_READ_VALUES) {
            int64_t delta;
            int64_t absolute;
            int64_t match_length;

            field_count = split_whitespace(row, fields, 2);
            if (field_count != 1 || !tv_parse_i64(fields[0], &delta)) {
                (void)invalid_delta(path, line_number,
                                    "expected one signed delta value");
                goto fail;
            }
            if (delta == 0) {
                int64_t remaining_query = query_span - query_consumed;
                int64_t remaining_target = target_span - target_consumed;
                int64_t columns_by_target;
                int64_t columns_by_query;
                char *canonical;
                int pushed;

                if (remaining_query < 0 || remaining_target < 0
                    || remaining_query != remaining_target
                    || append_operation(&alignment.ops, &alignment.n_ops,
                                        &operation_capacity, 'M',
                                        remaining_query) != 0
                    || alignment.n_ops == 0
                    || target_span > INT64_MAX - negative_deltas
                    || query_span > INT64_MAX - positive_deltas) {
                    (void)invalid_delta(path, line_number,
                                        "delta path does not consume declared spans");
                    goto fail;
                }
                columns_by_target = target_span + negative_deltas;
                columns_by_query = query_span + positive_deltas;
                if (columns_by_target <= 0
                    || columns_by_target != columns_by_query
                    || alignment.error_count < delta_count
                    || alignment.similarity_error_count < delta_count
                    || alignment.error_count > columns_by_target
                    || alignment.similarity_error_count > columns_by_target
                    || (columns_by_target <= INT64_MAX / 2
                        && alignment.nonalpha_count > 2 * columns_by_target)) {
                    (void)invalid_delta(path, line_number,
                                        "inconsistent error counts or alignment length");
                    goto fail;
                }
                if (alignment.strand == '-') {
                    reverse_operation_order(alignment.ops, alignment.n_ops);
                }
                alignment.block_len = columns_by_target;
                alignment.matches = columns_by_target - alignment.error_count;
                alignment.aggregate_identity =
                    (double)alignment.matches / (double)alignment.block_len;
                canonical = canonical_record(&alignment);
                if (canonical == NULL) {
                    (void)invalid_delta(path, line_number,
                                        "cannot canonicalize alignment");
                    goto fail;
                }
                pushed = tv_push_native_alignment(run, &alignment, canonical);
                free(canonical);
                if (pushed < 0) {
                    (void)invalid_delta(path, line_number,
                                        "cannot store alignment evidence");
                    goto fail;
                }
                if (pushed > 0) {
                    if (run->cfg.verbose) {
                        fprintf(stderr,
                                "tevox: ignored duplicate MUMmer delta evidence at %s:%d\n",
                                path, line_number);
                    }
                } else {
                    added_records++;
                }
                parsed_records++;
                records_in_sequence++;
                operation_capacity = 0;
                state = DELTA_EXPECT_RECORD;
                continue;
            }
            if (delta == INT64_MIN) {
                (void)invalid_delta(path, line_number,
                                    "delta magnitude overflows int64");
                goto fail;
            }
            absolute = delta < 0 ? -delta : delta;
            match_length = absolute - 1;
            if (match_length > query_span - query_consumed
                || match_length > target_span - target_consumed
                || append_operation(&alignment.ops, &alignment.n_ops,
                                    &operation_capacity, 'M',
                                    match_length) != 0) {
                (void)invalid_delta(path, line_number,
                                    "delta path exceeds declared spans");
                goto fail;
            }
            query_consumed += match_length;
            target_consumed += match_length;
            if (delta > 0) {
                if (target_consumed >= target_span
                    || positive_deltas == INT64_MAX
                    || append_operation(&alignment.ops, &alignment.n_ops,
                                        &operation_capacity, 'D', 1) != 0) {
                    (void)invalid_delta(path, line_number,
                                        "reference insertion exceeds span");
                    goto fail;
                }
                target_consumed++;
                positive_deltas++;
            } else {
                if (query_consumed >= query_span
                    || negative_deltas == INT64_MAX
                    || append_operation(&alignment.ops, &alignment.n_ops,
                                        &operation_capacity, 'I', 1) != 0) {
                    (void)invalid_delta(path, line_number,
                                        "query insertion exceeds span");
                    goto fail;
                }
                query_consumed++;
                negative_deltas++;
            }
            if (delta_count == INT64_MAX) {
                (void)invalid_delta(path, line_number,
                                    "too many delta values");
                goto fail;
            }
            delta_count++;
            continue;
        }
    }

    if (ferror(stream)) {
        tv_print_error("error reading MUMmer delta '%s'", path);
        goto fail;
    }
    if (state == DELTA_READ_VALUES) {
        (void)invalid_delta(path, line_number,
                            "unterminated alignment (missing zero)");
        goto fail;
    }
    if (state != DELTA_EXPECT_RECORD || !have_sequence_header
        || records_in_sequence == 0 || parsed_records == 0) {
        (void)invalid_delta(path, line_number, "incomplete or empty file");
        goto fail;
    }
    free(reference_name);
    free(query_name);
    free(line);
    fclose(stream);
    if (added_records == 0) {
        tv_print_error("MUMmer delta '%s' contains no new records", path);
        return -1;
    }
    return 0;

fail:
    free_alignment(&alignment);
    free(reference_name);
    free(query_name);
    free(line);
    fclose(stream);
    return -1;
}
