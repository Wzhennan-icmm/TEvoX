#ifndef TE_COMPARATOR_H
#define TE_COMPARATOR_H

// 为了兼容性定义
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>

// 数据结构定义
typedef struct {
    char* chr;
    int start;
    int end;
    char* strand;
} GenomicRegion;

typedef struct {
    char* id;
    char* chr;
    int start;
    int end;
    char* strand;
    char* type;
    char* family;
    char* name;
} Transposon;

typedef struct {
    char* chr1;
    int start1;
    int end1;
    char* chr2;
    int start2;
    int end2;
    double score;
} SyntenyBlock;

typedef struct {
    Transposon* transposons;
    int count;
    int capacity;
} TEList;

typedef struct {
    SyntenyBlock* blocks;
    int count;
    int capacity;
} SyntenyList;

typedef struct {
    char* type;
    int count;
} TypeCount;

typedef struct {
    char* family;
    int count;
} FamilyCount;

// 文件类型枚举
typedef enum {
    FILE_GFF3,
    FILE_BED,
    FILE_SYNTENY,
    FILE_UNKNOWN
} FileType;

// 主要函数声明
FileType detect_file_type(const char* filename);
int parse_gff3(const char* filename, TEList* te_list);
int parse_bed(const char* filename, TEList* te_list);
int parse_synteny(const char* filename, SyntenyList* synteny_list);
int compare_te_differences(TEList* te1, TEList* te2, SyntenyList* synteny, 
                           TEList* unique_te1, TEList* unique_te2);
void print_te_list(TEList* te_list, const char* title);
void print_synteny_list(SyntenyList* synteny_list, const char* title);
void free_te_list(TEList* te_list);
void free_synteny_list(SyntenyList* synteny_list);
bool is_in_synteny_region(Transposon* te, SyntenyList* synteny, int genome_id);
void analyze_te_types(TEList* unique_te1, TEList* unique_te2);
void analyze_te_families(TEList* unique_te1, TEList* unique_te2);
TypeCount* count_te_types(TEList* te_list);
FamilyCount* count_te_families(TEList* te_list);
void write_results_to_file(TEList* unique_te1, TEList* unique_te2, const char* output_prefix);
void init_te_list(TEList* te_list);
void init_synteny_list(SyntenyList* synteny_list);
void add_transposon(TEList* te_list, Transposon* te);
void add_synteny_block(SyntenyList* synteny_list, SyntenyBlock* block);

// 工具函数
char* strdup_safe(const char* str);
void* safe_malloc(size_t size);
void* safe_realloc(void* ptr, size_t size);

#endif // TE_COMPARATOR_H