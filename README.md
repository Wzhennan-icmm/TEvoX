# TEvo - Transposable Element Evolution Analyzer

A C-language tool for comparing transposon (TE) differences between two genomes using synteny information.

## Features

- Parse synteny files to identify conserved regions between genomes
- Support for GFF3 and BED format TE annotation files
- Identify unique transposons in each genome (those outside synteny regions)
- Statistical analysis of TE types and families
- Generate detailed output reports

## Usage

### Basic Usage

```bash
./te_comparator <synteny_file> <te_file1> <te_file2> [options]
```

### Required Arguments

- `synteny_file`: File containing synteny blocks between two genomes
- `te_file1`: Transposon annotation file for genome 1 (GFF3 or BED format)
- `te_file2`: Transposon annotation file for genome 2 (GFF3 or BED format)

### Optional Arguments

- `genome1_file`: Genome sequence file for genome 1 (reserved for future use)
- `genome2_file`: Genome sequence file for genome 2 (reserved for future use)

### Options

- `-o, --output PREFIX`: Output file prefix (default: te_comparison)
- `-v, --verbose`: Enable verbose output
- `-h, --help`: Show help message

### Examples

```bash
# Basic comparison
./te_comparator synteny.txt genome1.te.gff3 genome2.te.bed

# With custom output prefix and verbose mode
./te_comparator synteny.txt genome1.te.gff3 genome2.te.bed -o my_comparison -v

# With genome files (for future use)
./te_comparator synteny.txt genome1.te.gff3 genome2.te.bed genome1.fa genome2.fa
```

## Input File Formats

### Synteny File Format

Tab-separated file with the following columns:
```
chr1    start1    end1    chr2    start2    end2    [score]
```

Example:
```
chr1    1000    5000    chr1    1000    5000    0.95
chr1    6000    12000   chr1    6200    12200   0.88
```

### GFF3 Format

Standard GFF3 format with TE-related features. The program looks for features with types containing:
- transposable_element
- TE
- retrotransposon
- DNA_transposon
- LINE
- SINE
- LTR
- TIR
- MITE
- helitron

### BED Format

Standard BED format with at least 3 columns:
```
chr    start    end    [name]    [score]    [strand]    [type]
```

## Output Files

The program generates two output files:

1. `{prefix}_genome1_unique.txt`: Unique transposons in genome 1
2. `{prefix}_genome2_unique.txt`: Unique transposons in genome 2

Each file contains:
- ID: Transposon identifier
- Chr: Chromosome name
- Start: Start position (1-based)
- End: End position
- Strand: Strand information (+, -, or .)
- Type: Transposon type
- Family: Transposon family
- Name: Transposon name

## Building

```bash
make
```

## Testing

```bash
make test
# or
./test/run_tests.sh
```

## Dependencies

- GCC compiler
- Standard C library

## Algorithm

1. Parse synteny blocks to identify conserved regions
2. Parse TE annotation files for both genomes
3. For each transposon, check if it overlaps with any synteny region
4. Transposons outside synteny regions are considered "unique"
5. Generate statistics and output files

## License

This project is open source. Please check the license file for details.

## Contributing

Contributions are welcome! Please ensure all code follows the project coding standards and includes appropriate tests.