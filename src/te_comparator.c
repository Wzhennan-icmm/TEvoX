#include "te_comparator.h"

// 比较两个基因组间的TE差异
int compare_te_differences(TEList* te1, TEList* te2, SyntenyList* synteny, 
                           TEList* unique_te1, TEList* unique_te2) {
    if (!te1 || !te2 || !unique_te1 || !unique_te2) {
        fprintf(stderr, "Error: Invalid parameters for compare_te_differences\n");
        return -1;
    }
    
    init_te_list(unique_te1);
    init_te_list(unique_te2);
    
    printf("Comparing TE differences between two genomes...\n");
    printf("Genome 1: %d transposons\n", te1->count);
    printf("Genome 2: %d transposons\n", te2->count);
    printf("Synteny blocks: %d\n", synteny ? synteny->count : 0);
    
    // 查找基因组1中独有的转座子
    int unique_count_1 = 0;
    for (int i = 0; i < te1->count; i++) {
        Transposon* te = &te1->transposons[i];
        
        // 如果没有共线性信息，或者转座子不在共线性区域内，则认为是独有的
        bool in_synteny = false;
        if (synteny && synteny->count > 0) {
            in_synteny = is_in_synteny_region(te, synteny, 1);
        }
        
        if (!in_synteny) {
            add_transposon(unique_te1, te);
            unique_count_1++;
        }
    }
    
    // 查找基因组2中独有的转座子
    int unique_count_2 = 0;
    for (int i = 0; i < te2->count; i++) {
        Transposon* te = &te2->transposons[i];
        
        // 如果没有共线性信息，或者转座子不在共线性区域内，则认为是独有的
        bool in_synteny = false;
        if (synteny && synteny->count > 0) {
            in_synteny = is_in_synteny_region(te, synteny, 2);
        }
        
        if (!in_synteny) {
            add_transposon(unique_te2, te);
            unique_count_2++;
        }
    }
    
    printf("\n=== TE Difference Analysis Results ===\n");
    printf("Genome 1 unique transposons: %d (%.1f%% of total)\n", 
           unique_count_1, te1->count > 0 ? (100.0 * unique_count_1 / te1->count) : 0.0);
    printf("Genome 2 unique transposons: %d (%.1f%% of total)\n", 
           unique_count_2, te2->count > 0 ? (100.0 * unique_count_2 / te2->count) : 0.0);
    
    // 按类型统计独有转座子
    analyze_te_types(unique_te1, unique_te2);
    
    // 按家族统计独有转座子
    analyze_te_families(unique_te1, unique_te2);
    
    return unique_count_1 + unique_count_2;
}

// 按类型分析转座子
void analyze_te_types(TEList* unique_te1, TEList* unique_te2) {
    printf("\n=== Transposon Type Analysis ===\n");
    
    // 统计基因组1的转座子类型
    printf("Genome 1 unique TE types:\n");
    TypeCount* types1 = count_te_types(unique_te1);
    if (types1) {
        for (int i = 0; types1[i].type != NULL; i++) {
            printf("  %s: %d\n", types1[i].type, types1[i].count);
        }
        free(types1);
    }
    
    // 统计基因组2的转座子类型
    printf("Genome 2 unique TE types:\n");
    TypeCount* types2 = count_te_types(unique_te2);
    if (types2) {
        for (int i = 0; types2[i].type != NULL; i++) {
            printf("  %s: %d\n", types2[i].type, types2[i].count);
        }
        free(types2);
    }
}

// 按家族分析转座子
void analyze_te_families(TEList* unique_te1, TEList* unique_te2) {
    printf("\n=== Transposon Family Analysis ===\n");
    
    // 统计基因组1的转座子家族
    printf("Genome 1 unique TE families:\n");
    FamilyCount* families1 = count_te_families(unique_te1);
    if (families1) {
        for (int i = 0; families1[i].family != NULL; i++) {
            printf("  %s: %d\n", families1[i].family, families1[i].count);
        }
        free(families1);
    }
    
    // 统计基因组2的转座子家族
    printf("Genome 2 unique TE families:\n");
    FamilyCount* families2 = count_te_families(unique_te2);
    if (families2) {
        for (int i = 0; families2[i].family != NULL; i++) {
            printf("  %s: %d\n", families2[i].family, families2[i].count);
        }
        free(families2);
    }
}

// 统计转座子类型
TypeCount* count_te_types(TEList* te_list) {
    if (!te_list || te_list->count == 0) {
        return NULL;
    }
    
    // 首先统计有多少种不同的类型
    int type_capacity = 10;
    TypeCount* types = (TypeCount*)safe_malloc((type_capacity + 1) * sizeof(TypeCount));
    int type_count = 0;
    
    for (int i = 0; i < te_list->count; i++) {
        Transposon* te = &te_list->transposons[i];
        char* type = te->type ? te->type : "unknown";
        
        // 查找是否已经存在这个类型
        int found = 0;
        for (int j = 0; j < type_count; j++) {
            if (strcmp(types[j].type, type) == 0) {
                types[j].count++;
                found = 1;
                break;
            }
        }
        
        // 如果没找到，添加新类型
        if (!found) {
            if (type_count >= type_capacity) {
                type_capacity *= 2;
                types = (TypeCount*)safe_realloc(types, (type_capacity + 1) * sizeof(TypeCount));
            }
            
            types[type_count].type = strdup_safe(type);
            types[type_count].count = 1;
            type_count++;
        }
    }
    
    // 添加结束标记
    types[type_count].type = NULL;
    types[type_count].count = 0;
    
    return types;
}

// 统计转座子家族
FamilyCount* count_te_families(TEList* te_list) {
    if (!te_list || te_list->count == 0) {
        return NULL;
    }
    
    // 首先统计有多少种不同的家族
    int family_capacity = 10;
    FamilyCount* families = (FamilyCount*)safe_malloc((family_capacity + 1) * sizeof(FamilyCount));
    int family_count = 0;
    
    for (int i = 0; i < te_list->count; i++) {
        Transposon* te = &te_list->transposons[i];
        char* family = te->family ? te->family : "unknown";
        
        // 查找是否已经存在这个家族
        int found = 0;
        for (int j = 0; j < family_count; j++) {
            if (strcmp(families[j].family, family) == 0) {
                families[j].count++;
                found = 1;
                break;
            }
        }
        
        // 如果没找到，添加新家族
        if (!found) {
            if (family_count >= family_capacity) {
                family_capacity *= 2;
                families = (FamilyCount*)safe_realloc(families, (family_capacity + 1) * sizeof(FamilyCount));
            }
            
            families[family_count].family = strdup_safe(family);
            families[family_count].count = 1;
            family_count++;
        }
    }
    
    // 添加结束标记
    families[family_count].family = NULL;
    families[family_count].count = 0;
    
    return families;
}

// 将结果写入文件
void write_results_to_file(TEList* unique_te1, TEList* unique_te2, 
                           const char* output_prefix) {
    if (!unique_te1 || !unique_te2 || !output_prefix) {
        return;
    }
    
    // 写入基因组1的独有转座子
    char filename1[512];
    snprintf(filename1, sizeof(filename1), "%s_genome1_unique.txt", output_prefix);
    
    FILE* file1 = fopen(filename1, "w");
    if (file1) {
        fprintf(file1, "# Unique transposons in Genome 1\n");
        fprintf(file1, "# ID\tChr\tStart\tEnd\tStrand\tType\tFamily\tName\n");
        
        for (int i = 0; i < unique_te1->count; i++) {
            Transposon* te = &unique_te1->transposons[i];
            fprintf(file1, "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
                   te->id ? te->id : "N/A",
                   te->chr ? te->chr : "N/A",
                   te->start,
                   te->end,
                   te->strand ? te->strand : ".",
                   te->type ? te->type : "N/A",
                   te->family ? te->family : "N/A",
                   te->name ? te->name : "N/A");
        }
        fclose(file1);
        printf("Genome 1 unique TEs written to: %s\n", filename1);
    }
    
    // 写入基因组2的独有转座子
    char filename2[512];
    snprintf(filename2, sizeof(filename2), "%s_genome2_unique.txt", output_prefix);
    
    FILE* file2 = fopen(filename2, "w");
    if (file2) {
        fprintf(file2, "# Unique transposons in Genome 2\n");
        fprintf(file2, "# ID\tChr\tStart\tEnd\tStrand\tType\tFamily\tName\n");
        
        for (int i = 0; i < unique_te2->count; i++) {
            Transposon* te = &unique_te2->transposons[i];
            fprintf(file2, "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
                   te->id ? te->id : "N/A",
                   te->chr ? te->chr : "N/A",
                   te->start,
                   te->end,
                   te->strand ? te->strand : ".",
                   te->type ? te->type : "N/A",
                   te->family ? te->family : "N/A",
                   te->name ? te->name : "N/A");
        }
        fclose(file2);
        printf("Genome 2 unique TEs written to: %s\n", filename2);
    }
}