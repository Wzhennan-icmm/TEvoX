#include "te_comparator.h"
#include <ctype.h>

// 简单的strcasecmp替代函数
int strcasecmp_safe(const char* s1, const char* s2) {
    while (*s1 && *s2) {
        int c1 = tolower((unsigned char)*s1);
        int c2 = tolower((unsigned char)*s2);
        if (c1 != c2) {
            return c1 - c2;
        }
        s1++;
        s2++;
    }
    return tolower((unsigned char)*s1) - tolower((unsigned char)*s2);
}

// 安全内存分配函数
void* safe_malloc(size_t size) {
    void* ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "Error: Memory allocation failed for %zu bytes\n", size);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

// 安全内存重新分配函数
void* safe_realloc(void* ptr, size_t size) {
    void* new_ptr = realloc(ptr, size);
    if (!new_ptr) {
        fprintf(stderr, "Error: Memory reallocation failed for %zu bytes\n", size);
        exit(EXIT_FAILURE);
    }
    return new_ptr;
}

// 安全字符串复制函数
char* strdup_safe(const char* str) {
    if (!str) return NULL;
    char* new_str = strdup(str);
    if (!new_str) {
        fprintf(stderr, "Error: String duplication failed\n");
        exit(EXIT_FAILURE);
    }
    return new_str;
}

// 检测文件类型
FileType detect_file_type(const char* filename) {
    if (!filename) return FILE_UNKNOWN;
    
    FILE* file = fopen(filename, "r");
    if (!file) return FILE_UNKNOWN;
    
    char line[1024];
    if (fgets(line, sizeof(line), file) == NULL) {
        fclose(file);
        return FILE_UNKNOWN;
    }
    
    fclose(file);
    
    // 检查文件扩展名
    const char* ext = strrchr(filename, '.');
    if (!ext) return FILE_UNKNOWN;
    
    if (strcasecmp_safe(ext, ".gff3") == 0 || strcasecmp_safe(ext, ".gff") == 0) {
        return FILE_GFF3;
    } else if (strcasecmp_safe(ext, ".bed") == 0) {
        return FILE_BED;
    }
    
    // 检查文件内容
    file = fopen(filename, "r");
    if (!file) return FILE_UNKNOWN;
    
    FileType type = FILE_UNKNOWN;
    int line_count = 0;
    
    while (fgets(line, sizeof(line), file) && line_count < 10) {
        line_count++;
        
        // 跳过注释行
        if (line[0] == '#') continue;
        
        // 计算制表符数量
        int tabs = 0;
        for (int i = 0; line[i] != '\0'; i++) {
            if (line[i] == '\t') tabs++;
        }
        
        if (tabs >= 8) {
            type = FILE_GFF3;
            break;
        } else if (tabs >= 2) {
            type = FILE_BED;
            break;
        }
    }
    
    fclose(file);
    return type;
}

// 初始化TE列表
void init_te_list(TEList* te_list) {
    if (!te_list) return;
    te_list->transposons = NULL;
    te_list->count = 0;
    te_list->capacity = 0;
}

// 初始化共线性列表
void init_synteny_list(SyntenyList* synteny_list) {
    if (!synteny_list) return;
    synteny_list->blocks = NULL;
    synteny_list->count = 0;
    synteny_list->capacity = 0;
}

// 添加转座子到列表
void add_transposon(TEList* te_list, Transposon* te) {
    if (!te_list || !te) return;
    
    if (te_list->count >= te_list->capacity) {
        int new_capacity = te_list->capacity == 0 ? 100 : te_list->capacity * 2;
        te_list->transposons = (Transposon*)safe_realloc(te_list->transposons, 
                                                         new_capacity * sizeof(Transposon));
        te_list->capacity = new_capacity;
    }
    
    // 深拷贝转座子数据，避免浅拷贝导致的内存问题
    Transposon* new_te = &te_list->transposons[te_list->count];
    memset(new_te, 0, sizeof(Transposon));
    
    new_te->id = te->id ? strdup_safe(te->id) : NULL;
    new_te->chr = te->chr ? strdup_safe(te->chr) : NULL;
    new_te->start = te->start;
    new_te->end = te->end;
    new_te->strand = te->strand ? strdup_safe(te->strand) : NULL;
    new_te->type = te->type ? strdup_safe(te->type) : NULL;
    new_te->family = te->family ? strdup_safe(te->family) : NULL;
    new_te->name = te->name ? strdup_safe(te->name) : NULL;
    
    te_list->count++;
}

// 添加共线性区块到列表
void add_synteny_block(SyntenyList* synteny_list, SyntenyBlock* block) {
    if (!synteny_list || !block) return;
    
    if (synteny_list->count >= synteny_list->capacity) {
        int new_capacity = synteny_list->capacity == 0 ? 100 : synteny_list->capacity * 2;
        synteny_list->blocks = (SyntenyBlock*)safe_realloc(synteny_list->blocks, 
                                                           new_capacity * sizeof(SyntenyBlock));
        synteny_list->capacity = new_capacity;
    }
    
    // 深拷贝共线性区块数据，避免浅拷贝导致的内存问题
    SyntenyBlock* new_block = &synteny_list->blocks[synteny_list->count];
    memset(new_block, 0, sizeof(SyntenyBlock));
    
    new_block->chr1 = block->chr1 ? strdup_safe(block->chr1) : NULL;
    new_block->chr2 = block->chr2 ? strdup_safe(block->chr2) : NULL;
    new_block->start1 = block->start1;
    new_block->end1 = block->end1;
    new_block->start2 = block->start2;
    new_block->end2 = block->end2;
    new_block->score = block->score;
    
    synteny_list->count++;
}

// 释放TE列表内存
void free_te_list(TEList* te_list) {
    if (!te_list) return;
    
    for (int i = 0; i < te_list->count; i++) {
        Transposon* te = &te_list->transposons[i];
        free(te->id);
        free(te->chr);
        free(te->strand);
        free(te->type);
        free(te->family);
        free(te->name);
    }
    
    free(te_list->transposons);
    te_list->transposons = NULL;
    te_list->count = 0;
    te_list->capacity = 0;
}

// 释放共线性列表内存
void free_synteny_list(SyntenyList* synteny_list) {
    if (!synteny_list) return;
    
    for (int i = 0; i < synteny_list->count; i++) {
        SyntenyBlock* block = &synteny_list->blocks[i];
        free(block->chr1);
        free(block->chr2);
    }
    
    free(synteny_list->blocks);
    synteny_list->blocks = NULL;
    synteny_list->count = 0;
    synteny_list->capacity = 0;
}
