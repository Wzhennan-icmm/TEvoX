#ifndef TEVOX_H
#define TEVOX_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "sha256.h"

#define TEVOX_VERSION "0.5.0-alpha.2"
#define TEVOX_SCHEMA_VERSION "1.2.0"
#define TEVOX_MODEL_ID "BUILTIN_UNCALIBRATED_V1"
#define TEVOX_CALIBRATION_STATUS "UNCALIBRATED"
#define TEVOX_UNKNOWN "."

typedef struct {
    int64_t start;
    int64_t end;
} TvGapRun;

typedef struct {
    char *role;
    char *path;
    char sha256[TV_SHA256_HEX_LENGTH + 1];
} TvInputDigest;

typedef struct {
    char *name;
    int index;
    int64_t length;
    int64_t seq_offset;
    int line_bases;
    int line_bytes;
    TvGapRun *gap_runs;
    size_t n_gap_runs;
    size_t cap_gap_runs;
} TvContig;

typedef struct {
    int genome;
    char *id;
    char *contig;
    int contig_index;
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
    size_t node_offset;
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
    TV_ALIGNMENT_PAF,
    TV_ALIGNMENT_MUMMER_DELTA
} TvAlignmentProvider;

typedef enum {
    TV_MAPQ_OBSERVED,
    TV_MAPQ_MISSING_255,
    TV_MAPQ_NOT_PROVIDED
} TvMapqStatus;

typedef enum {
    TV_AGG_IDENTITY_MISSING,
    TV_AGG_IDENTITY_PAF_CORE,
    TV_AGG_IDENTITY_DELTA_ERRORS
} TvAggregateIdentityMethod;

typedef enum {
    TV_MAPPING_NOT_ESTABLISHED,
    TV_MAPPING_MAPQ_PASS,
    TV_MAPPING_UNIQUE_ALIGNMENT,
    TV_MAPPING_UNIQUE_COPY_CONTEXT,
    TV_MAPPING_AMBIGUOUS,
    TV_MAPPING_FAILED
} TvMappingConfidence;

typedef enum {
    TV_EVIDENCE_NATIVE,
    TV_EVIDENCE_DERIVED_REVERSE
} TvEvidenceOrigin;

typedef struct {
    TvAlignmentProvider provider;
    int query_genome;
    int target_genome;
    char *qname;
    int qcontig_index;
    int64_t qlen;
    int64_t qstart;
    int64_t qend;
    char strand;
    char *tname;
    int tcontig_index;
    int64_t tlen;
    int64_t tstart;
    int64_t tend;
    int64_t matches;
    int64_t block_len;
    int mapq;
    bool mapq_observed;
    TvMapqStatus mapq_status;
    double aggregate_identity;
    TvAggregateIdentityMethod aggregate_identity_method;
    int64_t error_count;
    int64_t similarity_error_count;
    int64_t nonalpha_count;
    TvCigarOp *ops;
    size_t n_ops;
    TvCigarOp *identity_ops;
    size_t n_identity_ops;
    TvIdentityMethod identity_method;
    TvEvidenceOrigin origin;
    char *source_path;
    int source_line;
    int source_record;
    uint64_t evidence_group_id;
} TvPaf;

typedef struct {
    int genome;
    char *id;
    char *contig;
    int contig_index;
    int64_t start;
    int64_t end;
    char strand;
    char *subgenome_id;
    char *haplotype_id;
} TvGene;

typedef struct {
    uint64_t block_id;
    uint64_t evidence_group_id;
    uint64_t homology_group_id;
    char *source_id;
    char *source_path;
    int source_record;
    int source_line;
    int genome_a;
    int genome_b;
    char *contig_a;
    char *contig_b;
    int contig_a_index;
    int contig_b_index;
    int64_t start_a;
    int64_t end_a;
    int64_t start_b;
    int64_t end_b;
    char orientation;
    int anchor_count;
    size_t anchor_start;
    double reported_score;
    double reported_evalue;
    char *wgd_node;
    int context_a;
    int context_b;
    char status[32];
    char reason[64];
} TvSyntenyBlock;

typedef struct {
    uint64_t anchor_id;
    int block_index;
    int rank;
    int provider_rank;
    int gene_a;
    int gene_b;
    double reported_evalue;
} TvSyntenyAnchor;

typedef struct {
    uint64_t context_id;
    uint64_t homology_group_id;
    int genome;
    int contig_index;
    char *contig;
    int64_t start;
    int64_t end;
    char *subgenome_id;
    char *haplotype_id;
    char *syntenic_copy_id;
    char *wgd_node;
    int anchor_count;
    char status[32];
} TvCopyContext;

typedef enum {
    TV_CONTEXT_BRACKETED,
    TV_CONTEXT_BLOCK_INTERIOR,
    TV_CONTEXT_AMBIGUOUS
} TvContextAssignment;

typedef struct {
    uint64_t te_context_id;
    int te_node;
    int context_index;
    TvContextAssignment assignment;
    double overlap_fraction;
    int left_anchor;
    int right_anchor;
} TvTEContext;

typedef enum {
    TV_CONTEXT_RELATION_UNKNOWN,
    TV_CONTEXT_RELATION_SUPPORTED,
    TV_CONTEXT_RELATION_CONFLICT,
    TV_CONTEXT_RELATION_AMBIGUOUS
} TvContextRelation;

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

typedef enum {
    TV_RELATION_ORTHOLOG,
    TV_RELATION_WGD_HOMEOLOG,
    TV_RELATION_ALLELIC,
    TV_RELATION_TANDEM_PARALOG,
    TV_RELATION_SEGMENTAL_PARALOG,
    TV_RELATION_TRANSPOSED_PARALOG,
    TV_RELATION_UNKNOWN,
    TV_RELATION_COUNT
} TvRelationClass;

typedef enum {
    TV_MATCH_NOT_APPLICABLE,
    TV_MATCH_OPTIMAL,
    TV_MATCH_HEURISTIC
} TvMatchingMethod;

typedef enum {
    TV_SOLVER_TRIVIAL,
    TV_SOLVER_EXACT_ENUMERATION,
    TV_SOLVER_DETERMINISTIC_GREEDY
} TvSolverMethod;

typedef enum {
    TV_SOLVER_OPTIMAL,
    TV_SOLVER_HEURISTIC
} TvSolverStatus;

enum {
    TV_FEATURE_FLANK = UINT32_C(1) << 0,
    TV_FEATURE_LOCAL_IDENTITY = UINT32_C(1) << 1,
    TV_FEATURE_AGGREGATE_IDENTITY = UINT32_C(1) << 2,
    TV_FEATURE_MAPQ = UINT32_C(1) << 3,
    TV_FEATURE_N_FRACTION = UINT32_C(1) << 4,
    TV_FEATURE_TE_ALIGNMENT = UINT32_C(1) << 5,
    TV_FEATURE_INSERTION = UINT32_C(1) << 6,
    TV_FEATURE_RECIPROCAL_OVERLAP = UINT32_C(1) << 7,
    TV_FEATURE_BOUNDARY = UINT32_C(1) << 8,
    TV_FEATURE_FAMILY = UINT32_C(1) << 9,
    TV_FEATURE_CONTEXT = UINT32_C(1) << 10
};

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
    TvMappingConfidence mapping_confidence;
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
    int target_contig_index;
    int nearby_candidate_count;
    int eligible_candidate_count;
    int retained_candidate_count;
    int graph_candidate_count;
    size_t candidate_start;
    size_t candidate_count;
    int selected_candidate_index;
    uint64_t evidence_id;
    double technical_scores[4];
    double biological_scores[4];
    double annotation_scores[5];
    double technical_entropy;
    double biological_entropy;
    double annotation_entropy;
    bool inference_out_of_domain;
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
    TvContextRelation context_relation;
    uint64_t shared_homology_group_id;
    bool context_compatible;
    bool eligible;
    bool selected;
    bool output_retained;
    bool graph_retained;
    uint32_t observed_feature_mask;
    double membership_logit;
    double membership_score;
    double membership_entropy;
    bool inference_out_of_domain;
    int rank;
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
    uint64_t edge_id;
    int a;
    int b;
    double score;
    double membership_logit;
    double membership_score;
    double membership_entropy;
    bool independent_reciprocal;
    bool family_compatible;
    TvContextRelation context_relation;
    uint64_t shared_homology_group_id;
    bool inference_out_of_domain;
    bool matching_selected;
    uint64_t matching_group_id;
    TvMatchingMethod matching_method;
    bool selected;
    int solver_component;
    int breakpoint_distance;
    size_t support_count;
    size_t support_start;
    size_t support_record_count;
    char selection_reason[64];
} TvEdge;

typedef struct {
    uint64_t solver_id;
    int node_count;
    int edge_count;
    TvSolverMethod method;
    TvSolverStatus status;
    double objective;
    double upper_bound;
    double relative_gap;
    uint64_t states_explored;
} TvSolverComponent;

typedef struct {
    uint64_t relation_id;
    int a;
    int b;
    int edge_index;
    int locus;
    TvRelationClass predicted;
    double scores[TV_RELATION_COUNT];
    double entropy;
    bool out_of_domain;
} TvRelation;

typedef struct {
    int key_a;
    int key_b;
    int key_c;
    int value;
    int64_t start;
    int64_t end;
    int64_t prefix_max_end;
} TvIntervalEntry;

typedef struct {
    size_t begin;
    size_t end;
} TvIndexRange;

typedef struct {
    int a;
    int b;
    uint64_t evidence_group_id;
    uint64_t evidence_id;
    bool direction_ab;
    bool native;
} TvEdgeSupport;

typedef struct {
    uint64_t paf_interval_queries;
    uint64_t paf_records_examined;
    uint64_t te_interval_queries;
    uint64_t te_records_examined;
    uint64_t gap_queries;
    uint64_t fasta_reopens_after_index;
    uint64_t edge_support_records;
} TvPerformanceCounters;

typedef struct {
    int flank;
    int candidate_window;
    int min_mapq;
    double min_flank_fraction;
    double min_edge_score;
    double max_n_fraction;
    double min_reciprocal_overlap;
    double near_best_delta;
    double min_delta_identity;
    double min_membership_score;
    double prediction_set_mass;
    int max_candidates;
    int max_graph_candidates;
    int exact_max_edges;
    int exact_matching_max_nodes;
    int tandem_distance;
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
    TvEdgeSupport *edge_support;
    size_t n_edge_support;
    size_t cap_edge_support;
    TvSolverComponent *solver_components;
    size_t n_solver_components;
    size_t cap_solver_components;
    TvRelation *relations;
    size_t n_relations;
    size_t cap_relations;
    TvGene *genes;
    size_t n_genes;
    size_t cap_genes;
    TvSyntenyBlock *synteny_blocks;
    size_t n_synteny_blocks;
    size_t cap_synteny_blocks;
    TvSyntenyAnchor *synteny_anchors;
    size_t n_synteny_anchors;
    size_t cap_synteny_anchors;
    TvCopyContext *contexts;
    size_t n_contexts;
    size_t cap_contexts;
    TvTEContext *te_contexts;
    size_t n_te_contexts;
    size_t cap_te_contexts;
    size_t *te_context_offsets;
    int *context_scratch_left;
    int *context_scratch_right;
    size_t context_scratch_capacity;
    char **synteny_source_paths;
    char **synteny_gene_paths;
    size_t n_synteny_source_paths;
    size_t cap_synteny_source_paths;
    TvIntervalEntry *paf_interval_index;
    size_t n_paf_interval_index;
    TvIntervalEntry *te_interval_index;
    size_t n_te_interval_index;
    int *components;
    int n_loci;
    TvInputDigest *input_digests;
    size_t n_input_digests;
    size_t cap_input_digests;
    TvPerformanceCounters performance;
    TvConfig cfg;
} TvRun;

void tv_run_init(TvRun *run);
void tv_run_free(TvRun *run);
int tv_load_manifest(TvRun *run, const char *path);
int tv_load_alignments(TvRun *run, const char *path);
int tv_load_synteny_sources(TvRun *run, const char *path);
int tv_register_input(TvRun *run, const char *role, const char *path);
int tv_verify_inputs(const TvRun *run);
const char *tv_input_path(const TvRun *run, const char *path);
const char *tv_input_sha256(const TvRun *run, const char *path);
int tv_add_genome(TvRun *run, const char *id, const char *fasta,
                  const char *annotation, int copies);
int tv_add_paf_file(TvRun *run, int query_genome, int target_genome,
                    const char *path);
int tv_add_mummer_delta_file(TvRun *run, int query_genome, int target_genome,
                             const char *path);
int tv_push_native_alignment(TvRun *run, TvPaf *alignment,
                             const char *canonical_record);
int tv_genome_id(const TvRun *run, const char *id);
TvContig *tv_find_contig(TvGenome *genome, const char *name);
void *tv_grow(void *pointer, size_t count, size_t width);
char *tv_dupstr(const char *text);
char *tv_directory(const char *path);
char *tv_resolve(const char *base, const char *path);
char *tv_strip(char *text);
int tv_split_tabs(char *text, char **fields, int capacity);
int tv_parse_i64(const char *text, int64_t *value);
int tv_parse_int(const char *text, int *value);
uint64_t tv_hash_bytes(uint64_t hash, const void *data, size_t length);
uint64_t tv_hash_text(uint64_t hash, const char *text);
int tv_build_synteny_contexts(TvRun *run);
TvContextRelation tv_candidate_context_relation(
    const TvRun *run, int source_te, int target_te,
    uint64_t *shared_homology_group_id);
const char *tv_context_merge_reason(const TvRun *run, int *parents,
                                    int left_root, int right_root);
int tv_build_interval_indexes(TvRun *run);
TvIndexRange tv_paf_index_range(TvRun *run, int query_genome,
                                int target_genome, int contig_index,
                                int64_t start, int64_t end);
TvIndexRange tv_te_index_range(TvRun *run, int genome, int contig_index,
                               int64_t start, int64_t end);
int tv_analyze(TvRun *run);
int tv_validate_output_prefix(const TvRun *run, const char *prefix);
int tv_write_outputs(TvRun *run, const char *prefix);
void tv_score_inference(TvRun *run);
int tv_infer_loci(TvRun *run);

const char *tv_state_name(TvState state);
const char *tv_technical_state_name(TvTechnicalState state);
const char *tv_biological_state_name(TvBiologicalState state);
const char *tv_annotation_state_name(TvAnnotationState state);
const char *tv_claim_type_name(TvClaimType type);
const char *tv_identity_method_name(TvIdentityMethod method);
const char *tv_alignment_provider_name(TvAlignmentProvider provider);
const char *tv_mapq_status_name(TvMapqStatus status);
const char *tv_aggregate_identity_method_name(
    TvAggregateIdentityMethod method);
const char *tv_mapping_confidence_name(TvMappingConfidence confidence);
const char *tv_evidence_origin_name(TvEvidenceOrigin origin);
const char *tv_family_relation_name(TvFamilyRelation relation);
const char *tv_context_assignment_name(TvContextAssignment assignment);
const char *tv_context_relation_name(TvContextRelation relation);
const char *tv_relation_class_name(TvRelationClass relation);
const char *tv_matching_method_name(TvMatchingMethod method);
const char *tv_solver_method_name(TvSolverMethod method);
const char *tv_solver_status_name(TvSolverStatus status);
void tv_print_error(const char *format, ...);

#endif
