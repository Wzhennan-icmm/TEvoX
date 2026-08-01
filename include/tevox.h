#ifndef TEVOX_H
#define TEVOX_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#define TEVOX_VERSION "0.2.0-alpha.1"
#define TEVOX_UNKNOWN "."
typedef struct { char *name; int64_t length, seq_offset; int line_bases, line_bytes; } TvContig;
typedef struct { int genome; char *id,*contig; int64_t start,end; char strand; char *type,*family; } TvTE;
typedef struct { char *id,*fasta_path,*te_path; int max_locus_copies; TvContig *contigs; size_t n_contigs,cap_contigs; TvTE *tes; size_t n_tes,cap_tes; } TvGenome;
typedef struct { char code; int64_t length; } TvCigarOp;
typedef struct { int query_genome,target_genome; char *qname; int64_t qlen,qstart,qend; char strand; char *tname; int64_t tlen,tstart,tend,matches,block_len; int mapq; TvCigarOp *ops; size_t n_ops; } TvPaf;
typedef enum { TV_PRESENT_ANNOTATED,TV_PRESENT_UNANNOTATED,TV_EMPTY_SITE_CONFIRMED,TV_STRUCTURAL_ALTERNATIVE,TV_FAMILY_OR_BOUNDARY_DISCORDANCE,TV_ASSEMBLY_GAP,TV_PROJECTION_AMBIGUOUS,TV_UNCALLABLE } TvState;
typedef struct { int source_te,target_genome,paf_index; TvState state; double quality,identity,left_flank,right_flank,n_fraction; int mapq; char *contig; int64_t start,end; int target_te,candidate_count; char evidence[96]; } TvProjection;
typedef struct { int a,b; double score; bool reciprocal,family_match,selected; int breakpoint_distance; } TvEdge;
typedef struct { int flank,candidate_window,min_mapq; double min_flank_fraction,min_edge_score,max_n_fraction; bool verbose; } TvConfig;
typedef struct { TvGenome *genomes; size_t n_genomes,cap_genomes; TvTE *nodes; size_t n_nodes; TvPaf *pafs; size_t n_pafs,cap_pafs; TvProjection *projections; size_t n_projections,cap_projections; TvEdge *edges; size_t n_edges,cap_edges; int *components,n_loci; TvConfig cfg; } TvRun;
void tv_run_init(TvRun*); void tv_run_free(TvRun*);
int tv_load_manifest(TvRun*,const char*); int tv_load_alignments(TvRun*,const char*);
int tv_add_genome(TvRun*,const char*,const char*,const char*,int);
int tv_add_paf_file(TvRun*,int,int,const char*);
int tv_analyze(TvRun*); int tv_write_outputs(TvRun*,const char*);
const char *tv_state_name(TvState); void tv_print_error(const char*,...);
#endif
