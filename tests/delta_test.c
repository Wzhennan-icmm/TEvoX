#include "tevox.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

static bool require(bool condition, const char *message)
{
    if (!condition) {
        (void)fprintf(stderr, "delta_test: %s\n", message);
    }
    return condition;
}

static bool operations_equal(const TvPaf *alignment,
                             const TvCigarOp *expected, size_t count)
{
    if (alignment->n_ops != count) {
        return false;
    }
    for (size_t index = 0; index < count; index++) {
        if (alignment->ops[index].code != expected[index].code
            || alignment->ops[index].length != expected[index].length) {
            return false;
        }
    }
    return true;
}

static int check_delta(const char *query_fasta, const char *query_annotation,
                       const char *target_fasta,
                       const char *target_annotation, const char *delta_path,
                       char expected_strand, const TvCigarOp *native_ops,
                       size_t native_op_count)
{
    static const TvCigarOp derived_ops[] = {
        {'M', 1}, {'I', 1}, {'M', 2}, {'D', 1}, {'M', 4}
    };
    TvRun run;
    int result = 1;

    tv_run_init(&run);
    if (!require(tv_add_genome(&run, "A", query_fasta,
                               query_annotation, 1) == 0,
                 "cannot load query fixture")
        || !require(tv_add_genome(&run, "B", target_fasta,
                                  target_annotation, 1) == 1,
                    "cannot load target fixture")
        || !require(tv_add_mummer_delta_file(&run, 0, 1, delta_path) == 0,
                    "cannot import valid delta fixture")
        || !require(run.n_pafs == 2,
                    "one native record must create two directed traversals")) {
        goto done;
    }

    TvPaf *native = &run.pafs[0];
    TvPaf *derived = &run.pafs[1];
    if (!require(native->provider == TV_ALIGNMENT_MUMMER_DELTA,
                 "native provider is not MUMMER_DELTA")
        || !require(native->origin == TV_EVIDENCE_NATIVE,
                    "first traversal is not native")
        || !require(derived->origin == TV_EVIDENCE_DERIVED_REVERSE,
                    "second traversal is not derived reverse")
        || !require(native->evidence_group_id == derived->evidence_group_id,
                    "native and derived traversals have different evidence groups")
        || !require(native->source_record == 1 && native->source_line == 4,
                    "native provider record/line is incorrect")
        || !require(derived->source_record == 1 && derived->source_line == 4,
                    "derived provider record/line is incorrect")
        || !require(native->qstart == 20 && native->qend == 28
                    && native->tstart == 10 && native->tend == 18,
                    "native 0-based half-open coordinates are incorrect")
        || !require(derived->qstart == 10 && derived->qend == 18
                    && derived->tstart == 20 && derived->tend == 28,
                    "derived 0-based half-open coordinates are incorrect")
        || !require(native->strand == expected_strand
                    && derived->strand == expected_strand,
                    "normalized strand is incorrect")
        || !require(operations_equal(native, native_ops, native_op_count),
                    "native normalized operation path is incorrect")
        || !require(operations_equal(
                        derived, derived_ops,
                        sizeof(derived_ops) / sizeof(derived_ops[0])),
                    "derived normalized operation path is incorrect")
        || !require(native->block_len == 9 && native->matches == 7
                    && native->error_count == 2,
                    "delta alignment column/error accounting is incorrect")
        || !require(fabs(native->aggregate_identity - 7.0 / 9.0) < 1e-12,
                    "native aggregate identity is not 7/9")
        || !require(fabs(derived->aggregate_identity - 7.0 / 9.0) < 1e-12,
                    "derived aggregate identity is not 7/9")
        || !require(native->aggregate_identity_method
                        == TV_AGG_IDENTITY_DELTA_ERRORS
                    && native->identity_method == TV_IDENTITY_MISSING,
                    "delta identity provenance is incorrect")
        || !require(!native->mapq_observed
                    && native->mapq_status == TV_MAPQ_NOT_PROVIDED,
                    "delta MAPQ must be unobserved")) {
        goto done;
    }
    result = 0;

done:
    tv_run_free(&run);
    return result;
}

int main(int argc, char **argv)
{
    static const TvCigarOp plus_ops[] = {
        {'M', 1}, {'D', 1}, {'M', 2}, {'I', 1}, {'M', 4}
    };
    static const TvCigarOp reverse_ops[] = {
        {'M', 4}, {'I', 1}, {'M', 2}, {'D', 1}, {'M', 1}
    };

    if (argc != 7) {
        (void)fprintf(stderr,
                      "usage: delta_test QUERY_FA QUERY_GFF TARGET_FA TARGET_GFF PLUS_DELTA REVERSE_DELTA\n");
        return 2;
    }
    if (check_delta(argv[1], argv[2], argv[3], argv[4], argv[5], '+',
                    plus_ops, sizeof(plus_ops) / sizeof(plus_ops[0])) != 0) {
        return 1;
    }
    return check_delta(argv[1], argv[2], argv[3], argv[4], argv[6], '-',
                       reverse_ops,
                       sizeof(reverse_ops) / sizeof(reverse_ops[0]));
}
