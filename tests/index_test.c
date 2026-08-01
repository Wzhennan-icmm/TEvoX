#include "tevox.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void *tv_grow(void *pointer, size_t count, size_t width)
{
    void *result = realloc(pointer, (count == 0 ? 1 : count) * width);

    assert(result != NULL);
    return result;
}

void tv_print_error(const char *format, ...)
{
    va_list arguments;

    va_start(arguments, format);
    (void)vfprintf(stderr, format, arguments);
    fputc('\n', stderr);
    va_end(arguments);
}

static void check_paf_index(TvRun *run)
{
    TvIndexRange range;

    assert(run->n_paf_interval_index == 6);
    for (size_t index = 1; index < 4; index++) {
        assert(run->paf_interval_index[index - 1].start
               <= run->paf_interval_index[index].start);
    }

    range = tv_paf_index_range(run, 0, 1, 0, 60, 100);
    assert(range.end - range.begin == 3);
    assert(run->paf_interval_index[range.begin].value == 3);
    assert(run->paf_interval_index[range.end - 1].value == 0);

    range = tv_paf_index_range(run, 0, 1, 0, 21, 49);
    assert(range.end - range.begin == 1);
    assert(run->paf_interval_index[range.begin].value == 3);

    range = tv_paf_index_range(run, 0, 1, 0, 201, 300);
    assert(range.begin == range.end);
    range = tv_paf_index_range(run, 0, 1, 0, -10, 5);
    assert(range.begin == range.end);

    range = tv_paf_index_range(run, 0, 1, 0, 60, 60);
    assert(range.end - range.begin == 2);
    range = tv_paf_index_range(run, 0, 1, 0, 100, 60);
    assert(range.begin == 0 && range.end == 0);
    range = tv_paf_index_range(run, 0, 0, 0, 0, 100);
    assert(range.begin == 0 && range.end == 0);

    range = tv_paf_index_range(run, 0, 1, 1, 0, 5);
    assert(range.end - range.begin == 1);
    assert(run->paf_interval_index[range.begin].value == 4);
    assert(run->paf_interval_index[range.begin].prefix_max_end == 5);
    assert(run->performance.paf_interval_queries == 8);
    assert(run->performance.paf_records_examined == 0);
}

static void check_te_index(TvRun *run)
{
    TvIndexRange range;

    assert(run->n_te_interval_index == 4);
    range = tv_te_index_range(run, 0, 0, 20, 20);
    assert(range.end - range.begin == 3);
    assert(run->te_interval_index[range.begin].value == 1);
    assert(run->te_interval_index[range.end - 1].value == 0);

    range = tv_te_index_range(run, 1, 0, 10, 11);
    assert(range.end - range.begin == 1);
    assert(run->te_interval_index[range.begin].value == 3);
    assert(run->te_interval_index[range.begin].prefix_max_end == 12);

    range = tv_te_index_range(run, -1, 0, 0, 10);
    assert(range.begin == 0 && range.end == 0);
    assert(run->performance.te_interval_queries == 3);
    assert(run->performance.te_records_examined == 0);
}

int main(void)
{
    TvRun run;
    TvGenome genomes[2];
    TvContig contigs_a[2];
    TvContig contigs_b[1];
    TvPaf pafs[6];
    TvTE nodes[4];

    memset(&run, 0, sizeof(run));
    memset(genomes, 0, sizeof(genomes));
    memset(contigs_a, 0, sizeof(contigs_a));
    memset(contigs_b, 0, sizeof(contigs_b));
    memset(pafs, 0, sizeof(pafs));
    memset(nodes, 0, sizeof(nodes));
    genomes[0].contigs = contigs_a;
    genomes[0].n_contigs = 2;
    genomes[1].contigs = contigs_b;
    genomes[1].n_contigs = 1;
    run.genomes = genomes;
    run.n_genomes = 2;

    pafs[0] = (TvPaf){.query_genome = 0, .target_genome = 1,
                      .qcontig_index = 0, .qstart = 100, .qend = 200};
    pafs[1] = (TvPaf){.query_genome = 0, .target_genome = 1,
                      .qcontig_index = 0, .qstart = 10, .qend = 20};
    pafs[2] = (TvPaf){.query_genome = 0, .target_genome = 1,
                      .qcontig_index = 0, .qstart = 50, .qend = 60};
    pafs[3] = (TvPaf){.query_genome = 0, .target_genome = 1,
                      .qcontig_index = 0, .qstart = 15, .qend = 150};
    pafs[4] = (TvPaf){.query_genome = 0, .target_genome = 1,
                      .qcontig_index = 1, .qstart = 0, .qend = 5};
    pafs[5] = (TvPaf){.query_genome = 1, .target_genome = 0,
                      .qcontig_index = 0, .qstart = 2, .qend = 3};
    run.pafs = pafs;
    run.n_pafs = 6;

    nodes[0] = (TvTE){.genome = 0, .contig_index = 0,
                      .start = 20, .end = 30};
    nodes[1] = (TvTE){.genome = 0, .contig_index = 0,
                      .start = 0, .end = 20};
    nodes[2] = (TvTE){.genome = 0, .contig_index = 0,
                      .start = 10, .end = 40};
    nodes[3] = (TvTE){.genome = 1, .contig_index = 0,
                      .start = 8, .end = 12};
    run.nodes = nodes;
    run.n_nodes = 4;

    assert(tv_build_interval_indexes(&run) == 0);
    check_paf_index(&run);
    check_te_index(&run);
    assert(tv_build_interval_indexes(&run) == 0);
    assert(run.n_paf_interval_index == 6);
    assert(run.n_te_interval_index == 4);
    free(run.paf_interval_index);
    free(run.te_interval_index);
    return 0;
}
