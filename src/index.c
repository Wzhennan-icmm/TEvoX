#include "tevox.h"

#include <limits.h>
#include <stdlib.h>

static int compare_int(int left, int right)
{
    if (left < right) {
        return -1;
    }
    if (left > right) {
        return 1;
    }
    return 0;
}

static int compare_i64(int64_t left, int64_t right)
{
    if (left < right) {
        return -1;
    }
    if (left > right) {
        return 1;
    }
    return 0;
}

static int compare_interval_entry(const void *left, const void *right)
{
    const TvIntervalEntry *a = left;
    const TvIntervalEntry *b = right;
    int comparison;

    comparison = compare_int(a->key_a, b->key_a);
    if (comparison != 0) {
        return comparison;
    }
    comparison = compare_int(a->key_b, b->key_b);
    if (comparison != 0) {
        return comparison;
    }
    comparison = compare_int(a->key_c, b->key_c);
    if (comparison != 0) {
        return comparison;
    }
    comparison = compare_i64(a->start, b->start);
    if (comparison != 0) {
        return comparison;
    }
    comparison = compare_i64(a->end, b->end);
    if (comparison != 0) {
        return comparison;
    }
    return compare_int(a->value, b->value);
}

static int compare_entry_key(const TvIntervalEntry *entry, int key_a,
                             int key_b, int key_c)
{
    int comparison = compare_int(entry->key_a, key_a);

    if (comparison != 0) {
        return comparison;
    }
    comparison = compare_int(entry->key_b, key_b);
    if (comparison != 0) {
        return comparison;
    }
    return compare_int(entry->key_c, key_c);
}

static void build_prefix_max_end(TvIntervalEntry *entries, size_t count)
{
    int64_t maximum = 0;
    int previous_a = 0;
    int previous_b = 0;
    int previous_c = 0;

    for (size_t index = 0; index < count; index++) {
        TvIntervalEntry *entry = &entries[index];
        bool new_key = index == 0 || entry->key_a != previous_a
            || entry->key_b != previous_b || entry->key_c != previous_c;

        if (new_key || entry->end > maximum) {
            maximum = entry->end;
        }
        entry->prefix_max_end = maximum;
        previous_a = entry->key_a;
        previous_b = entry->key_b;
        previous_c = entry->key_c;
    }
}

static bool flattened_node_count(const TvRun *run, size_t *count)
{
    size_t total = 0;

    for (size_t genome = 0; genome < run->n_genomes; genome++) {
        size_t genome_count = run->genomes[genome].n_tes;

        if (genome_count > SIZE_MAX - total) {
            return false;
        }
        total += genome_count;
    }
    *count = total;
    return true;
}

static int build_paf_index(const TvRun *run, TvIntervalEntry **output,
                           size_t *output_count)
{
    TvIntervalEntry *entries = NULL;

    *output = NULL;
    *output_count = 0;
    if (run->n_pafs == 0) {
        return 0;
    }
    if (run->n_pafs > (size_t)INT_MAX) {
        tv_print_error("too many alignment views for the interval index");
        return -1;
    }
    entries = tv_grow(NULL, run->n_pafs, sizeof(*entries));
    for (size_t index = 0; index < run->n_pafs; index++) {
        const TvPaf *paf = &run->pafs[index];

        entries[index] = (TvIntervalEntry){
            .key_a = paf->query_genome,
            .key_b = paf->target_genome,
            .key_c = paf->qcontig_index,
            .value = (int)index,
            .start = paf->qstart,
            .end = paf->qend
        };
    }
    qsort(entries, run->n_pafs, sizeof(*entries), compare_interval_entry);
    build_prefix_max_end(entries, run->n_pafs);
    *output = entries;
    *output_count = run->n_pafs;
    return 0;
}

static int build_te_index(const TvRun *run, TvIntervalEntry **output,
                          size_t *output_count)
{
    TvIntervalEntry *entries = NULL;
    size_t count;

    *output = NULL;
    *output_count = 0;
    if (run->n_genomes > (size_t)INT_MAX) {
        tv_print_error("too many genomes for the interval index");
        return -1;
    }
    if (run->nodes != NULL) {
        count = run->n_nodes;
    } else if (!flattened_node_count(run, &count)) {
        tv_print_error("too many TE annotations for the interval index");
        return -1;
    }
    if (count == 0) {
        return 0;
    }
    if (count > (size_t)INT_MAX) {
        tv_print_error("too many TE annotations for the interval index");
        return -1;
    }
    entries = tv_grow(NULL, count, sizeof(*entries));
    if (run->nodes != NULL) {
        for (size_t index = 0; index < count; index++) {
            const TvTE *te = &run->nodes[index];

            entries[index] = (TvIntervalEntry){
                .key_a = te->genome,
                .key_b = -1,
                .key_c = te->contig_index,
                .value = (int)index,
                .start = te->start,
                .end = te->end
            };
        }
    } else {
        size_t node = 0;

        for (size_t genome = 0; genome < run->n_genomes; genome++) {
            const TvGenome *source = &run->genomes[genome];

            for (size_t local = 0; local < source->n_tes; local++) {
                const TvTE *te = &source->tes[local];

                entries[node] = (TvIntervalEntry){
                    .key_a = (int)genome,
                    .key_b = -1,
                    .key_c = te->contig_index,
                    .value = (int)node,
                    .start = te->start,
                    .end = te->end
                };
                node++;
            }
        }
    }
    qsort(entries, count, sizeof(*entries), compare_interval_entry);
    build_prefix_max_end(entries, count);
    *output = entries;
    *output_count = count;
    return 0;
}

int tv_build_interval_indexes(TvRun *run)
{
    TvIntervalEntry *paf_entries = NULL;
    TvIntervalEntry *te_entries = NULL;
    size_t paf_count = 0;
    size_t te_count = 0;

    if (run == NULL) {
        return -1;
    }
    if (build_paf_index(run, &paf_entries, &paf_count) != 0
        || build_te_index(run, &te_entries, &te_count) != 0) {
        free(paf_entries);
        free(te_entries);
        return -1;
    }
    free(run->paf_interval_index);
    free(run->te_interval_index);
    run->paf_interval_index = paf_entries;
    run->n_paf_interval_index = paf_count;
    run->te_interval_index = te_entries;
    run->n_te_interval_index = te_count;
    return 0;
}

static size_t key_lower_bound(const TvIntervalEntry *entries, size_t count,
                              int key_a, int key_b, int key_c)
{
    size_t low = 0;
    size_t high = count;

    while (low < high) {
        size_t middle = low + (high - low) / 2;

        if (compare_entry_key(&entries[middle], key_a, key_b, key_c) < 0) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    return low;
}

static size_t key_upper_bound(const TvIntervalEntry *entries, size_t count,
                              int key_a, int key_b, int key_c)
{
    size_t low = 0;
    size_t high = count;

    while (low < high) {
        size_t middle = low + (high - low) / 2;

        if (compare_entry_key(&entries[middle], key_a, key_b, key_c) <= 0) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    return low;
}

static TvIndexRange interval_range(const TvIntervalEntry *entries, size_t count,
                                   int key_a, int key_b, int key_c,
                                   int64_t start, int64_t end)
{
    TvIndexRange empty = {0, 0};
    size_t key_begin;
    size_t key_end;
    size_t candidate_begin;
    size_t candidate_end;
    size_t low;
    size_t high;

    if (entries == NULL || count == 0 || start > end) {
        return empty;
    }
    key_begin = key_lower_bound(entries, count, key_a, key_b, key_c);
    if (key_begin == count
        || compare_entry_key(&entries[key_begin], key_a, key_b, key_c) != 0) {
        return empty;
    }
    key_end = key_upper_bound(entries, count, key_a, key_b, key_c);

    /* Keep entries whose start may touch the inclusive query envelope. */
    low = key_begin;
    high = key_end;
    while (low < high) {
        size_t middle = low + (high - low) / 2;

        if (entries[middle].start <= end) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    candidate_end = low;

    /* Prefix maxima are monotone inside a key group. Entries before this
       point all end before the inclusive query envelope and cannot match. */
    low = key_begin;
    high = candidate_end;
    while (low < high) {
        size_t middle = low + (high - low) / 2;

        if (entries[middle].prefix_max_end < start) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    candidate_begin = low;
    if (candidate_begin == candidate_end) {
        return (TvIndexRange){candidate_end, candidate_end};
    }
    return (TvIndexRange){candidate_begin, candidate_end};
}

static void increment_counter(uint64_t *counter)
{
    if (*counter != UINT64_MAX) {
        (*counter)++;
    }
}

TvIndexRange tv_paf_index_range(TvRun *run, int query_genome,
                                int target_genome, int contig_index,
                                int64_t start, int64_t end)
{
    TvIndexRange range = {0, 0};

    if (run == NULL) {
        return range;
    }
    increment_counter(&run->performance.paf_interval_queries);
    if (query_genome < 0 || target_genome < 0
        || (size_t)query_genome >= run->n_genomes
        || (size_t)target_genome >= run->n_genomes
        || query_genome == target_genome || contig_index < 0
        || (size_t)contig_index
           >= run->genomes[query_genome].n_contigs) {
        return range;
    }
    range = interval_range(run->paf_interval_index,
                           run->n_paf_interval_index, query_genome,
                           target_genome, contig_index, start, end);
    return range;
}

TvIndexRange tv_te_index_range(TvRun *run, int genome, int contig_index,
                               int64_t start, int64_t end)
{
    TvIndexRange range = {0, 0};

    if (run == NULL) {
        return range;
    }
    increment_counter(&run->performance.te_interval_queries);
    if (genome < 0 || (size_t)genome >= run->n_genomes || contig_index < 0
        || (size_t)contig_index >= run->genomes[genome].n_contigs) {
        return range;
    }
    range = interval_range(run->te_interval_index,
                           run->n_te_interval_index, genome, -1,
                           contig_index, start, end);
    return range;
}
