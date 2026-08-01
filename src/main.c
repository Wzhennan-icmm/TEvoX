#include "tevox.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char *manifest;
    const char *alignments;
    const char *output;
    const char *genome_a;
    const char *fasta_a;
    const char *te_a;
    const char *genome_b;
    const char *fasta_b;
    const char *te_b;
    const char *paf;
    int copies_a;
    int copies_b;
} Arguments;

static void usage(FILE *stream)
{
    (void)fprintf(
        stream,
        "TEvoX %s - uncertainty-aware TE locus reconstruction\n\n"
        "Usage:\n"
        "  tevox graph --manifest genomes.tsv --alignments alignments.tsv [options]\n"
        "  tevox pair --genome-a A --fasta-a A.fa --te-a A.gff3 "
        "--genome-b B --fasta-b B.fa --te-b B.gff3 "
        "--paf A_query_B_target.paf [options]\n\n"
        "PAF direction:\n"
        "  PAF query must be genome A and PAF target must be genome B. Because\n"
        "  minimap2 takes target first, generate it with: minimap2 [opts] B.fa A.fa\n\n"
        "Options:\n"
        "  -o, --output PREFIX       output prefix (default: tevox)\n"
        "  --flank INT               flank length (100)\n"
        "  --candidate-window INT    candidate window (100)\n"
        "  --min-mapq INT            minimum observed MAPQ (20; 255 is missing)\n"
        "  --min-flank FLOAT         paired-flank fraction (0.60)\n"
        "  --min-recip-overlap FLOAT minimum reciprocal TE overlap (0.50)\n"
        "  --near-best-delta FLOAT   ambiguity score delta (5.0)\n"
        "  --min-edge FLOAT          graph edge threshold (45)\n"
        "  --max-n FLOAT             maximum N fraction (0.25)\n"
        "  --max-copies-a/b INT      legacy pair-mode component quota\n"
        "  -v, --verbose\n",
        TEVOX_VERSION);
}

static const char *option_value(int argc, char **argv, int *index)
{
    if (*index + 1 >= argc) {
        tv_print_error("option '%s' requires a value", argv[*index]);
        return NULL;
    }
    (*index)++;
    return argv[*index];
}

static bool parse_long(const char *text, long *value)
{
    char *end = NULL;

    errno = 0;
    *value = strtol(text, &end, 10);
    return errno == 0 && end != text && *end == '\0';
}

static bool parse_double(const char *text, double *value)
{
    char *end = NULL;

    errno = 0;
    *value = strtod(text, &end);
    return errno == 0 && end != text && *end == '\0';
}

static int invalid_value(const char *option, const char *value)
{
    tv_print_error("invalid value '%s' for option '%s'", value, option);
    return -1;
}

static int parse_options(int argc, char **argv, Arguments *arguments,
                         TvConfig *config)
{
    for (int index = 2; index < argc; index++) {
        const char *option = argv[index];
        const char *value;
        long integer;
        double decimal;

        if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
            return 1;
        }
        if (strcmp(option, "-v") == 0 || strcmp(option, "--verbose") == 0) {
            config->verbose = true;
            continue;
        }
        value = option_value(argc, argv, &index);
        if (value == NULL) {
            return -1;
        }
        if (strcmp(option, "--manifest") == 0) {
            arguments->manifest = value;
        } else if (strcmp(option, "--alignments") == 0) {
            arguments->alignments = value;
        } else if (strcmp(option, "-o") == 0
                   || strcmp(option, "--output") == 0) {
            arguments->output = value;
        } else if (strcmp(option, "--genome-a") == 0) {
            arguments->genome_a = value;
        } else if (strcmp(option, "--fasta-a") == 0) {
            arguments->fasta_a = value;
        } else if (strcmp(option, "--te-a") == 0) {
            arguments->te_a = value;
        } else if (strcmp(option, "--genome-b") == 0) {
            arguments->genome_b = value;
        } else if (strcmp(option, "--fasta-b") == 0) {
            arguments->fasta_b = value;
        } else if (strcmp(option, "--te-b") == 0) {
            arguments->te_b = value;
        } else if (strcmp(option, "--paf") == 0) {
            arguments->paf = value;
        } else if (strcmp(option, "--flank") == 0) {
            if (!parse_long(value, &integer) || integer < 1
                || integer > 1000000) {
                return invalid_value(option, value);
            }
            config->flank = (int)integer;
        } else if (strcmp(option, "--candidate-window") == 0) {
            if (!parse_long(value, &integer) || integer < 0
                || integer > 1000000) {
                return invalid_value(option, value);
            }
            config->candidate_window = (int)integer;
        } else if (strcmp(option, "--min-mapq") == 0) {
            if (!parse_long(value, &integer) || integer < 0 || integer > 254) {
                return invalid_value(option, value);
            }
            config->min_mapq = (int)integer;
        } else if (strcmp(option, "--max-copies-a") == 0) {
            if (!parse_long(value, &integer) || integer < 1 || integer > 1000) {
                return invalid_value(option, value);
            }
            arguments->copies_a = (int)integer;
        } else if (strcmp(option, "--max-copies-b") == 0) {
            if (!parse_long(value, &integer) || integer < 1 || integer > 1000) {
                return invalid_value(option, value);
            }
            arguments->copies_b = (int)integer;
        } else if (strcmp(option, "--min-flank") == 0) {
            if (!parse_double(value, &decimal) || decimal < 0.0
                || decimal > 1.0) {
                return invalid_value(option, value);
            }
            config->min_flank_fraction = decimal;
        } else if (strcmp(option, "--min-recip-overlap") == 0) {
            if (!parse_double(value, &decimal) || decimal < 0.0
                || decimal > 1.0) {
                return invalid_value(option, value);
            }
            config->min_reciprocal_overlap = decimal;
        } else if (strcmp(option, "--near-best-delta") == 0) {
            if (!parse_double(value, &decimal) || decimal < 0.0
                || decimal > 100.0) {
                return invalid_value(option, value);
            }
            config->near_best_delta = decimal;
        } else if (strcmp(option, "--min-edge") == 0) {
            if (!parse_double(value, &decimal) || decimal < 0.0
                || decimal > 100.0) {
                return invalid_value(option, value);
            }
            config->min_edge_score = decimal;
        } else if (strcmp(option, "--max-n") == 0) {
            if (!parse_double(value, &decimal) || decimal < 0.0
                || decimal > 1.0) {
                return invalid_value(option, value);
            }
            config->max_n_fraction = decimal;
        } else {
            tv_print_error("unknown option '%s'", option);
            return -1;
        }
    }
    return 0;
}

static bool pair_arguments_complete(const Arguments *arguments)
{
    return arguments->genome_a != NULL && arguments->fasta_a != NULL
        && arguments->te_a != NULL && arguments->genome_b != NULL
        && arguments->fasta_b != NULL && arguments->te_b != NULL
        && arguments->paf != NULL;
}

int main(int argc, char **argv)
{
    bool graph;
    bool pair;
    TvRun run;
    Arguments arguments = {
        .output = "tevox",
        .copies_a = 1,
        .copies_b = 1
    };
    int parsed;
    int status = 0;

    if (argc == 2 && strcmp(argv[1], "--version") == 0) {
        (void)printf("tevox %s (schema %s)\n", TEVOX_VERSION,
                     TEVOX_SCHEMA_VERSION);
        return 0;
    }
    if (argc < 2 || strcmp(argv[1], "-h") == 0
        || strcmp(argv[1], "--help") == 0) {
        usage(argc < 2 ? stderr : stdout);
        return argc < 2 ? 2 : 0;
    }
    graph = strcmp(argv[1], "graph") == 0;
    pair = strcmp(argv[1], "pair") == 0;
    if (!graph && !pair) {
        tv_print_error("expected subcommand 'graph' or 'pair'");
        return 2;
    }

    tv_run_init(&run);
    parsed = parse_options(argc, argv, &arguments, &run.cfg);
    if (parsed == 1) {
        usage(stdout);
        tv_run_free(&run);
        return 0;
    }
    if (parsed < 0) {
        status = -1;
    } else if (graph) {
        if (arguments.manifest == NULL || arguments.alignments == NULL) {
            tv_print_error("graph requires --manifest and --alignments");
            status = -1;
        } else if (tv_load_manifest(&run, arguments.manifest) != 0
                   || tv_load_alignments(&run, arguments.alignments) != 0) {
            status = -1;
        }
    } else if (!pair_arguments_complete(&arguments)) {
        tv_print_error("pair requires both genomes, FASTAs, annotations and --paf");
        status = -1;
    } else {
        int genome_a = tv_add_genome(&run, arguments.genome_a,
                                     arguments.fasta_a, arguments.te_a,
                                     arguments.copies_a);
        int genome_b = genome_a < 0 ? -1
            : tv_add_genome(&run, arguments.genome_b, arguments.fasta_b,
                            arguments.te_b, arguments.copies_b);

        if (genome_a < 0 || genome_b < 0
            || tv_add_paf_file(&run, genome_a, genome_b, arguments.paf) != 0) {
            status = -1;
        }
    }
    if (status == 0
        && (tv_analyze(&run) != 0
            || tv_write_outputs(&run, arguments.output) != 0)) {
        status = -1;
    }
    if (status == 0) {
        (void)printf("TEvoX reconstructed %d locus/loci across %zu genomes.\n",
                     run.n_loci, run.n_genomes);
        (void)printf(
            "Evidence schema: %s; outputs: "
            "%s.{evidence,candidates,decisions,edges,loci,instances,states,summary}.tsv "
            "and %s.run.json\n",
            TEVOX_SCHEMA_VERSION, arguments.output, arguments.output);
    }
    tv_run_free(&run);
    return status == 0 ? 0 : 1;
}
