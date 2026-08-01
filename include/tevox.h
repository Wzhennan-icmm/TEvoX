#ifndef TEVOX_H
#define TEVOX_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define TEVOX_VERSION "0.3.0-alpha.1"
#define TEVOX_SCHEMA_VERSION "1.0.0"
#define TEVOX_UNKNOWN "."

typedef struct {
    char *name;
    int64_t length;
    int64_t seq_offset;
    int line_bases;
    int line_bytes;
} TvContig;

typedef struct {
    int genome;
    char *id;
    char *contig;
    int64_t start;
    int64_t end;
    char strand;
    char *type;
    char *family;
} TvTE;

typedef struct {
    char *id;
    char *fasta_path;
    char *te_path;
    int max_locus_copies;
    TvContig *contigs;
    size_t n_contigs;
    size_t cap_contigs;
    TvTE *tes;
    size_t n_tes;
    size_t cap_tes;
} TvGenome;

typedef struct {
    char code;
    int64_t length;
} TvCigarOp;

typedef enum {
    TV_IDENTITY_MISSING,
    TV_IDENTITY_EQX,
    TV_IDENTITY_CS
} TvIdentityMethod;

typedef enum {
    TV_EVIDENCE_NATIVE,
    TV_EVIDENCE_DERIVED_REVERSE
} TvEvidenceOrigin;

typedef struct {
    int query_genome;
    int target_genome;
    char *qname;
    int64_t qlen;
    int64_t qstart;
    int64_t qend;
    char strand;
    char *tname;
    int64_t tlen;
    int64_t tstart;
    int64_t tend;
    int64_t matches;
    int64_t block_len;
    int mapq;
    bool mapq_observed;
    TvCigarOp *ops;
    size_t n_ops;
    TvCigarOp *identity_ops;
    size_t n_identity_ops;
    TvIdentityMethod identity_method;
    TvEvidenceOrigin origin;
    char *source_path;
    int source_line;
    uint64_t evidence_group_id;
} TvPaf;

typedef enum {
    TV_PRESENT_ANNOTATED,
    TV_PRESENT_UNANNOTATED,
    TV_EMPTY_SITE_CONFIRMED,
    TV_STRUCTURAL_ALTERNATIVE,
    TV_FAMILY_OR_BOUNDARY_DISCORDANCE,
    TV_ASSEMBLY_GAP,
    TV_PROJECTION_AMBIGUOUS,
    TV_UNCALLABLE
} TvState;

typedef enum {
    TV_TECH_CALLABLE,
    TV_TECH_GAP,
    TV_TECH_AMBIGUOUS,
    TV_TECH_UNCALLABLE
} TvTechnicalState;

typedef enum {
    TV_BIO_PRESENT,
    TV_BIO_EMPTY,
    TV_BIO_STRUCTURAL_ALTERNATIVE,
    TV_BIO_UNKNOWN
} TvBiologicalState;

typedef enum {
    TV_ANN_MATCHED,
    TV_ANN_MISSING,
    TV_ANN_FAMILY_CONFLICT,
    TV_ANN_NOT_APPLICABLE,
    TV_ANN_UNKNOWN
} TvAnnotationState;

typedef enum {
    TV_CLAIM_PRESENCE,
    TV_CLAIM_EMPTY_SITE,
    TV_CLAIM_STRUCTURAL_ALTERNATIVE,
    TV_CLAIM_ANNOTATION_DISCORDANCE,
    TV_CLAIM_NONE
} TvClaimType;

typedef enum {
    TV_FAMILY_MATCH,
    TV_FAMILY_CONFLICT,
    TV_FAMILY_UNKNOWN
} TvFamilyRelation;

typedef struct {
    int source_te;
    int target_genome;
    int paf_index;
    int decision_index;
    TvState state;
    TvTechnicalState technical_state;
    TvBiologicalState biological_state;
    TvAnnotationState annotation_state;
    TvClaimType claim_type;
    bool claimable;
    bool primary;
    bool near_best;
    bool decision_ambiguous;
    int rank;
    double quality;
    double evidence_completeness;
    double identity;
    double left_flank;
    double right_flank;
    double n_fraction;
    double te_aligned_fraction;
    double insertion_fraction;
    int mapq;
    bool mapq_observed;
    bool left_flank_observed;
    bool right_flank_observed;
    char *contig;
    int64_t start;
    int64_t end;
    int target_te;
    int nearby_candidate_count;
    int eligible_candidate_count;
    uint64_t evidence_id;
    char decision_code[96];
    char claimability_reason[96];
} TvProjection;

typedef struct {
    int projection_index;
    int target_te;
    double score;
    double reciprocal_overlap;
    double boundary_score;
    int breakpoint_distance;
    TvFamilyRelation family_relation;
    bool eligible;
    bool selected;
    uint64_t candidate_id;
    char decision_code[64];
} TvCandidate;

typedef struct {
    int source_te;
    int target_genome;
    size_t projection_start;
    size_t projection_count;
    int winner_projection;
    int near_best_count;
    TvState state;
    TvTechnicalState technical_state;
    TvBiologicalState biological_state;
    TvAnnotationState annotation_state;
    TvClaimType claim_type;
    bool claimable;
    bool ambiguous;
    double quality;
    double evidence_completeness;
    uint64_t decision_id;
    char decision_code[96];
    char claimability_reason[96];
} TvDecision;

typedef struct {
    int a;
    int b;
    double score;
    bool independent_reciprocal;
    bool family_compatible;
    bool selected;
    int breakpoint_distance;
    size_t support_count;
    char selection_reason[64];
} TvEdge;

typedef struct {
    int flank;
    int candidate_window;
    int min_mapq;
    double min_flank_fraction;
    double min_edge_score;
    double max_n_fraction;
    double min_reciprocal_overlap;
    double near_best_delta;
    bool verbose;
} TvConfig;

typedef struct {
    TvGenome *genomes;
    size_t n_genomes;
    size_t cap_genomes;
    TvTE *nodes;
    size_t n_nodes;
    TvPaf *pafs;
    size_t n_pafs;
    size_t cap_pafs;
    TvProjection *projections;
    size_t n_projections;
    size_t cap_projections;
    TvCandidate *candidates;
    size_t n_candidates;
    size_t cap_candidates;
    TvDecision *decisions;
    size_t n_decisions;
    size_t cap_decisions;
    TvEdge *edges;
    size_t n_edges;
    size_t cap_edges;
    int *components;
    int n_loci;
    TvConfig cfg;
} TvRun;

void tv_run_init(TvRun *run);
void tv_run_free(TvRun *run);
int tv_load_manifest(TvRun *run, const char *path);
int tv_load_alignments(TvRun *run, const char *path);
int tv_add_genome(TvRun *run, const char *id, const char *fasta,
                  const char *annotation, int copies);
int tv_add_paf_file(TvRun *run, int query_genome, int target_genome,
                    const char *path);
int tv_analyze(TvRun *run);
int tv_write_outputs(TvRun *run, const char *prefix);

const char *tv_state_name(TvState state);
const char *tv_technical_state_name(TvTechnicalState state);
const char *tv_biological_state_name(TvBiologicalState state);
const char *tv_annotation_state_name(TvAnnotationState state);
const char *tv_claim_type_name(TvClaimType type);
const char *tv_identity_method_name(TvIdentityMethod method);
const char *tv_evidence_origin_name(TvEvidenceOrigin origin);
const char *tv_family_relation_name(TvFamilyRelation relation);
void tv_print_error(const char *format, ...);

#endif
