#include "te_comparator.h"

// 程序使用说明
void print_usage(const char* program_name) {
    printf("TE Comparator - Compare transposon differences between two genomes\n\n");
    printf("Usage: %s <synteny_file> <te_file1> <te_file2> [genome1_file] [genome2_file] [options]\n\n", program_name);
    printf("Required arguments:\n");
    printf("  synteny_file    File containing synteny blocks between two genomes\n");
    printf("  te_file1        Transposon annotation file for genome 1 (GFF3 or BED format)\n");
    printf("  te_file2        Transposon annotation file for genome 2 (GFF3 or BED format)\n\n");
    printf("Optional arguments:\n");
    printf("  genome1_file    Genome sequence file for genome 1 (for future use)\n");
    printf("  genome2_file    Genome sequence file for genome 2 (for future use)\n\n");
    printf("Options:\n");
    printf("  -o, --output PREFIX    Output file prefix (default: te_comparison)\n");
    printf("  -v, --verbose          Enable verbose output\n");
    printf("  -h, --help             Show this help message\n\n");
    printf("Examples:\n");
    printf("  %s synteny.txt genome1.te.gff3 genome2.te.bed\n", program_name);
    printf("  %s synteny.txt genome1.te.gff3 genome2.te.bed -o my_comparison\n", program_name);
    printf("\n");
}

// 解析命令行参数
typedef struct {
    char* synteny_file;
    char* te_file1;
    char* te_file2;
    char* genome1_file;
    char* genome2_file;
    char* output_prefix;
    bool verbose;
    bool show_help;
} ProgramArgs;

void init_args(ProgramArgs* args) {
    args->synteny_file = NULL;
    args->te_file1 = NULL;
    args->te_file2 = NULL;
    args->genome1_file = NULL;
    args->genome2_file = NULL;
    args->output_prefix = strdup_safe("te_comparison");
    args->verbose = false;
    args->show_help = false;
}

void free_args(ProgramArgs* args) {
    free(args->output_prefix);
}
int parse_arguments(int argc, char* argv[], ProgramArgs* args) {
    // 首先检查帮助选项
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            args->show_help = true;
            return 0;
        }
    }
    
    if (argc < 4) {
        return -1; // 参数不足
    }
    
    // 必需参数
    args->synteny_file = strdup_safe(argv[1]);
    args->te_file1 = strdup_safe(argv[2]);
    args->te_file2 = strdup_safe(argv[3]);
    
    // 可选的位置参数
    int positional_args = 4;
    if (argc > positional_args && argv[positional_args][0] != '-') {
        args->genome1_file = strdup_safe(argv[positional_args]);
        positional_args++;
    }
    
    if (argc > positional_args && argv[positional_args][0] != '-') {
        args->genome2_file = strdup_safe(argv[positional_args]);
        positional_args++;
    }
    
    // 解析选项参数
    for (int i = positional_args; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            args->show_help = true;
            return 0;
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            args->verbose = true;
        } else if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            free(args->output_prefix);
            args->output_prefix = strdup_safe(argv[++i]);
        } else {
            fprintf(stderr, "Error: Unknown option %s\n", argv[i]);
            return -1;
        }
    }
    
    return 0;
}

// 检查文件是否存在
bool file_exists(const char* filename) {
    if (!filename) return false;
    FILE* file = fopen(filename, "r");
    if (file) {
        fclose(file);
        return true;
    }
    return false;
}

// 验证输入参数
int validate_arguments(ProgramArgs* args) {
    if (!args->synteny_file || !args->te_file1 || !args->te_file2) {
        fprintf(stderr, "Error: Missing required arguments\n");
        return -1;
    }
    
    // 检查必需文件是否存在
    if (!file_exists(args->synteny_file)) {
        fprintf(stderr, "Error: Synteny file not found: %s\n", args->synteny_file);
        return -1;
    }
    
    if (!file_exists(args->te_file1)) {
        fprintf(stderr, "Error: TE file 1 not found: %s\n", args->te_file1);
        return -1;
    }
    
    if (!file_exists(args->te_file2)) {
        fprintf(stderr, "Error: TE file 2 not found: %s\n", args->te_file2);
        return -1;
    }
    
    // 检查可选文件
    if (args->genome1_file && !file_exists(args->genome1_file)) {
        fprintf(stderr, "Warning: Genome 1 file not found: %s\n", args->genome1_file);
    }
    
    if (args->genome2_file && !file_exists(args->genome2_file)) {
        fprintf(stderr, "Warning: Genome 2 file not found: %s\n", args->genome2_file);
    }
    
    return 0;
}

int main(int argc, char* argv[]) {
    ProgramArgs args;
    init_args(&args);
    
    // 解析命令行参数
    int parse_result = parse_arguments(argc, argv, &args);
    if (parse_result != 0 || args.show_help) {
        free_args(&args);
        return args.show_help ? 0 : 1;
    }
    
    // 验证参数
    if (validate_arguments(&args) != 0) {
        free_args(&args);
        return 1;
    }
    
    printf("=== TE Comparator ===\n");
    printf("Synteny file: %s\n", args.synteny_file);
    printf("TE file 1: %s\n", args.te_file1);
    printf("TE file 2: %s\n", args.te_file2);
    if (args.genome1_file) printf("Genome 1 file: %s\n", args.genome1_file);
    if (args.genome2_file) printf("Genome 2 file: %s\n", args.genome2_file);
    printf("Output prefix: %s\n", args.output_prefix);
    printf("Verbose mode: %s\n", args.verbose ? "ON" : "OFF");
    printf("\n");
    
    // 解析共线性文件
    SyntenyList synteny_list;
    if (parse_synteny(args.synteny_file, &synteny_list) < 0) {
        fprintf(stderr, "Error: Failed to parse synteny file\n");
        free_args(&args);
        return 1;
    }
    
    if (args.verbose) {
        print_synteny_list(&synteny_list, "Synteny Blocks");
    }
    
    // 解析TE文件1
    TEList te_list1;
    FileType type1 = detect_file_type(args.te_file1);
    int result1 = 0;
    
    if (type1 == FILE_GFF3) {
        result1 = parse_gff3(args.te_file1, &te_list1);
    } else if (type1 == FILE_BED) {
        result1 = parse_bed(args.te_file1, &te_list1);
    } else {
        fprintf(stderr, "Error: Unsupported file format for TE file 1: %s\n", args.te_file1);
        free_args(&args);
        free_synteny_list(&synteny_list);
        return 1;
    }
    
    if (result1 < 0) {
        fprintf(stderr, "Error: Failed to parse TE file 1\n");
        free_args(&args);
        free_synteny_list(&synteny_list);
        return 1;
    }
    
    if (args.verbose) {
        print_te_list(&te_list1, "Genome 1 Transposons");
    }
    
    // 解析TE文件2
    TEList te_list2;
    FileType type2 = detect_file_type(args.te_file2);
    int result2 = 0;
    
    if (type2 == FILE_GFF3) {
        result2 = parse_gff3(args.te_file2, &te_list2);
    } else if (type2 == FILE_BED) {
        result2 = parse_bed(args.te_file2, &te_list2);
    } else {
        fprintf(stderr, "Error: Unsupported file format for TE file 2: %s\n", args.te_file2);
        free_args(&args);
        free_synteny_list(&synteny_list);
        free_te_list(&te_list1);
        return 1;
    }
    
    if (result2 < 0) {
        fprintf(stderr, "Error: Failed to parse TE file 2\n");
        free_args(&args);
        free_synteny_list(&synteny_list);
        free_te_list(&te_list1);
        return 1;
    }
    
    if (args.verbose) {
        print_te_list(&te_list2, "Genome 2 Transposons");
    }
    
    // 比较TE差异
    TEList unique_te1, unique_te2;
    int total_unique = compare_te_differences(&te_list1, &te_list2, &synteny_list, &unique_te1, &unique_te2);
    
    if (args.verbose) {
        print_te_list(&unique_te1, "Genome 1 Unique Transposons");
        print_te_list(&unique_te2, "Genome 2 Unique Transposons");
    }
    
    // 写入结果文件
    write_results_to_file(&unique_te1, &unique_te2, args.output_prefix);
    
    printf("\n=== Analysis Complete ===\n");
    printf("Total unique transposons identified: %d\n", total_unique);
    printf("Results written to files with prefix: %s\n", args.output_prefix);
    
    // 清理内存
    free_args(&args);
    free_synteny_list(&synteny_list);
    free_te_list(&te_list1);
    free_te_list(&te_list2);
    free_te_list(&unique_te1);
    free_te_list(&unique_te2);
    
    return 0;
}