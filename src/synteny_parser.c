#include "te_comparator.h"

// 解析共线性文件
int parse_synteny(const char* filename, SyntenyList* synteny_list) {
    if (!filename || !synteny_list) {
        fprintf(stderr, "Error: Invalid parameters for parse_synteny\n");
        return -1;
    }
    
    init_synteny_list(synteny_list);
    
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open synteny file %s\n", filename);
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
        
        // 解析共线性区块
        // 假设格式为: chr1    start1    end1    chr2    start2    end2    [score]
        char* token;
        char* tokens[10];
        int token_count = 0;
        
        token = strtok(line, "\t");
        while (token != NULL && token_count < 10) {
            tokens[token_count++] = token;
            token = strtok(NULL, "\t");
        }
        
        if (token_count < 6) {
            fprintf(stderr, "Warning: Line %d has insufficient columns (%d), skipping\n", 
                    line_num, token_count);
            continue;
        }
        
        SyntenyBlock block;
        block.chr1 = strdup_safe(tokens[0]);
        block.start1 = atoi(tokens[1]);
        block.end1 = atoi(tokens[2]);
        block.chr2 = strdup_safe(tokens[3]);
        block.start2 = atoi(tokens[4]);
        block.end2 = atoi(tokens[5]);
        
        // 可选的score字段
        if (token_count > 6) {
            block.score = atof(tokens[6]);
        } else {
            block.score = 1.0;
        }
        
        // 验证坐标的合理性
        if (block.start1 > block.end1) {
            int temp = block.start1;
            block.start1 = block.end1;
            block.end1 = temp;
        }
        
        if (block.start2 > block.end2) {
            int temp = block.start2;
            block.start2 = block.end2;
            block.end2 = temp;
        }
        
        add_synteny_block(synteny_list, &block);
    }
    
    fclose(file);
    
    printf("Parsed %d synteny blocks from %s\n", synteny_list->count, filename);
    return synteny_list->count;
}

// 检查转座子是否在共线性区域内
bool is_in_synteny_region(Transposon* te, SyntenyList* synteny, int genome_id) {
    if (!te || !synteny || synteny->count == 0) {
        return false;
    }
    
    for (int i = 0; i < synteny->count; i++) {
        SyntenyBlock* block = &synteny->blocks[i];
        
        if (genome_id == 1) {
            // 检查基因组1的转座子
            if (strcmp(te->chr, block->chr1) == 0) {
                // 检查是否有重叠
                if (!(te->end < block->start1 || te->start > block->end1)) {
                    return true;
                }
            }
        } else if (genome_id == 2) {
            // 检查基因组2的转座子
            if (strcmp(te->chr, block->chr2) == 0) {
                // 检查是否有重叠
                if (!(te->end < block->start2 || te->start > block->end2)) {
                    return true;
                }
            }
        }
    }
    
    return false;
}

// 打印共线性信息（用于调试）
void print_synteny_list(SyntenyList* synteny_list, const char* title) {
    if (!synteny_list || !title) return;
    
    printf("\n=== %s ===\n", title);
    printf("Total synteny blocks: %d\n", synteny_list->count);
    
    if (synteny_list->count > 0) {
        printf("First 5 blocks:\n");
        printf("Chr1\tStart1\tEnd1\tChr2\tStart2\tEnd2\tScore\n");
        int max_print = synteny_list->count < 5 ? synteny_list->count : 5;
        
        for (int i = 0; i < max_print; i++) {
            SyntenyBlock* block = &synteny_list->blocks[i];
            printf("%s\t%d\t%d\t%s\t%d\t%d\t%.2f\n",
                   block->chr1, block->start1, block->end1,
                   block->chr2, block->start2, block->end2,
                   block->score);
        }
    }
    printf("\n");
}