#include "te_comparator.h"

// 解析GFF3文件的属性字段
void parse_gff3_attributes(const char* attributes_str, Transposon* te) {
    if (!attributes_str || !te) return;
    
    char* attrs = strdup_safe(attributes_str);
    char* attr = strtok(attrs, ";");
    
    while (attr != NULL) {
        // 跳过空格
        while (*attr == ' ') attr++;
        
        char* equal = strchr(attr, '=');
        if (equal != NULL) {
            *equal = '\0';
            char* key = attr;
            char* value = equal + 1;
            
            if (strcmp(key, "ID") == 0) {
                te->id = strdup_safe(value);
            } else if (strcmp(key, "Name") == 0) {
                te->name = strdup_safe(value);
            } else if (strcmp(key, "family") == 0) {
                te->family = strdup_safe(value);
            } else if (strcmp(key, "TE_family") == 0) {
                te->family = strdup_safe(value);
            }
        }
        
        attr = strtok(NULL, ";");
    }
    
    free(attrs);
}

// 解析GFF3文件
int parse_gff3(const char* filename, TEList* te_list) {
    if (!filename || !te_list) {
        fprintf(stderr, "Error: Invalid parameters for parse_gff3\n");
        return -1;
    }
    
    init_te_list(te_list);
    
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open GFF3 file %s\n", filename);
        return -1;
    }
    
    char line[2048];
    int line_num = 0;
    
    while (fgets(line, sizeof(line), file)) {
        line_num++;
        
        // 跳过注释行和空行
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
            continue;
        }
        
        // 移除行尾的换行符
        line[strcspn(line, "\r\n")] = '\0';
        
        // 解析GFF3行（9列）
        char* token;
        char* tokens[10];
        int token_count = 0;
        
        token = strtok(line, "\t");
        while (token != NULL && token_count < 10) {
            tokens[token_count++] = token;
            token = strtok(NULL, "\t");
        }
        
        if (token_count < 9) {
            fprintf(stderr, "Warning: Line %d has insufficient columns (%d), skipping\n", 
                    line_num, token_count);
            continue;
        }
        
        // 检查是否为转座子相关特征
        char* feature_type = tokens[2];
        bool is_te = false;
        
        // 常见的转座子特征类型
        const char* te_types[] = {
            "transposable_element", "TE", "retrotransposon", "DNA_transposon",
            "LINE", "SINE", "LTR", "TIR", "MITE", "helitron"
        };
        
        for (int i = 0; i < 10; i++) {
            if (strstr(feature_type, te_types[i]) != NULL) {
                is_te = true;
                break;
            }
        }
        
        // 如果不是转座子特征，跳过
        if (!is_te) {
            continue;
        }
        
        Transposon te;
        memset(&te, 0, sizeof(te));
        
        te.chr = strdup_safe(tokens[0]);
        te.start = atoi(tokens[3]);
        te.end = atoi(tokens[4]);
        te.strand = strdup_safe(tokens[6]);
        te.type = strdup_safe(tokens[2]);
        
        // 解析属性字段
        parse_gff3_attributes(tokens[8], &te);
        
        // 如果没有ID，生成一个
        if (!te.id) {
            char temp_id[100];
            snprintf(temp_id, sizeof(temp_id), "TE_%d_%d_%d", line_num, te.start, te.end);
            te.id = strdup_safe(temp_id);
        }
        
        // 如果没有name，使用ID
        if (!te.name) {
            te.name = strdup_safe(te.id);
        }
        
        add_transposon(te_list, &te);
    }
    
    fclose(file);
    
    printf("Parsed %d transposons from GFF3 file %s\n", te_list->count, filename);
    return te_list->count;
}

// 解析BED文件
int parse_bed(const char* filename, TEList* te_list) {
    if (!filename || !te_list) {
        fprintf(stderr, "Error: Invalid parameters for parse_bed\n");
        return -1;
    }
    
    init_te_list(te_list);
    
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open BED file %s\n", filename);
        return -1;
    }
    
    char line[2048];
    int line_num = 0;
    
    while (fgets(line, sizeof(line), file)) {
        line_num++;
        
        // 跳过注释行和空行
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r' || 
            (line[0] == 't' && line[1] == 'r' && line[2] == 'a' && line[3] == 'c' && line[4] == 'k')) {
            continue;
        }
        
        // 移除行尾的换行符
        line[strcspn(line, "\r\n")] = '\0';
        
        // 解析BED行（至少3列）
        char* token;
        char* tokens[15];
        int token_count = 0;
        
        token = strtok(line, "\t");
        while (token != NULL && token_count < 15) {
            tokens[token_count++] = token;
            token = strtok(NULL, "\t");
        }
        
        if (token_count < 3) {
            fprintf(stderr, "Warning: Line %d has insufficient columns (%d), skipping\n", 
                    line_num, token_count);
            continue;
        }
        
        Transposon te;
        memset(&te, 0, sizeof(te));
        
        te.chr = strdup_safe(tokens[0]);
        te.start = atoi(tokens[1]) + 1; // BED是0-based，转换为1-based
        te.end = atoi(tokens[2]);
        
        // 可选字段
        if (token_count > 3) {
            te.name = strdup_safe(tokens[3]);
        } else {
            // 生成默认名称
            char temp_name[100];
            snprintf(temp_name, sizeof(temp_name), "TE_%d_%d_%d", line_num, te.start, te.end);
            te.name = strdup_safe(temp_name);
        }
        
        if (token_count > 5) {
            te.strand = strdup_safe(tokens[5]);
        } else {
            te.strand = strdup_safe(".");
        }
        
        if (token_count > 6) {
            te.type = strdup_safe(tokens[6]);
        } else {
            te.type = strdup_safe("transposable_element");
        }
        
        // 生成ID
        char temp_id[100];
        snprintf(temp_id, sizeof(temp_id), "TE_%d_%d_%d", line_num, te.start, te.end);
        te.id = strdup_safe(temp_id);
        
        add_transposon(te_list, &te);
    }
    
    fclose(file);
    
    printf("Parsed %d transposons from BED file %s\n", te_list->count, filename);
    return te_list->count;
}

// 打印转座子列表（用于调试）
void print_te_list(TEList* te_list, const char* title) {
    if (!te_list || !title) return;
    
    printf("\n=== %s ===\n", title);
    printf("Total transposons: %d\n", te_list->count);
    
    if (te_list->count > 0) {
        printf("First 10 transposons:\n");
        printf("ID\tChr\tStart\tEnd\tStrand\tType\tFamily\n");
        int max_print = te_list->count < 10 ? te_list->count : 10;
        
        for (int i = 0; i < max_print; i++) {
            Transposon* te = &te_list->transposons[i];
            printf("%s\t%s\t%d\t%d\t%s\t%s\t%s\n",
                   te->id ? te->id : "N/A",
                   te->chr ? te->chr : "N/A",
                   te->start,
                   te->end,
                   te->strand ? te->strand : ".",
                   te->type ? te->type : "N/A",
                   te->family ? te->family : "N/A");
        }
    }
    printf("\n");
}