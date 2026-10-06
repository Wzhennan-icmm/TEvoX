#!/usr/bin/env python3
"""Annotate native TEvoX candidate TSV using exact GFF3 gene feature intervals.

TEvoX Start/End are ONE-based inclusive. Gene distances are the smallest base
coordinate difference between the candidate and a gene body: overlap = 0,
bookended intervals = 1. All equidistant nearest genes are reported. Gene
proximity is descriptive and does not imply altered expression or causation.
"""
import argparse
from bisect import bisect_left, bisect_right
import collections
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import sys
from normalize_annotations import attributes, open_text, protect_inputs, read_lengths, validate_interval


@dataclass(frozen=True)
class Gene:
    start: int
    end: int
    identifier: str
    name: str
    strand: str


class GeneIndex:
    def __init__(self, genes):
        self.by_start = sorted(genes, key=lambda g: (g.start, g.end, g.identifier))
        self.starts = [g.start for g in self.by_start]
        self.prefix_ends = []
        maxend = -1
        for gene in self.by_start:
            maxend = max(maxend, gene.end)
            self.prefix_ends.append(maxend)
        self.by_end = sorted(genes, key=lambda g: (g.end, g.start, g.identifier))
        self.ends = [g.end for g in self.by_end]

    def query(self, start, end):
        overlaps = []
        pos = bisect_left(self.starts, end) - 1
        # Prefix maxima find long/nested genes missed by checking adjacent starts.
        while pos >= 0 and self.prefix_ends[pos] > start:
            gene = self.by_start[pos]
            if gene.end > start:
                overlaps.append(gene)
            pos -= 1
        if overlaps:
            return overlaps, overlaps, 0
        candidates = []
        left = bisect_right(self.ends, start) - 1
        if left >= 0:
            nearest_end = self.ends[left]
            while left >= 0 and self.ends[left] == nearest_end:
                candidates.append((start - nearest_end + 1, self.by_end[left]))
                left -= 1
        right = bisect_left(self.starts, end)
        if right < len(self.starts):
            nearest_start = self.starts[right]
            while right < len(self.starts) and self.starts[right] == nearest_start:
                candidates.append((nearest_start - end + 1, self.by_start[right]))
                right += 1
        if not candidates:
            return [], [], None
        distance = min(item[0] for item in candidates)
        return [], [gene for gap, gene in candidates if gap == distance], distance


def load_genes(path, lengths, reject_path):
    genes, counts = collections.defaultdict(list), collections.Counter()
    with open_text(path) as handle, open(reject_path, "w") as reject_file:
        rejects = csv.writer(reject_file, delimiter="\t", lineterminator="\n")
        rejects.writerow(["source_line", "reason", "source_record"])
        for lineno, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            counts["source_feature_rows"] += 1
            if len(fields) == 9 and fields[2] != "gene":
                counts["other_feature_rows_skipped"] += 1
                continue
            try:
                if len(fields) != 9:
                    raise ValueError("gff_requires_nine_columns")
                chrom, start, end = fields[0], int(fields[3]) - 1, int(fields[4])
                strand = "." if fields[6] == "?" else fields[6]
                validate_interval(chrom, start, end, strand, lengths)
                attr = attributes(fields[8])
                if not attr.get("ID") or attr["ID"] == ".":
                    raise ValueError("gene_missing_id")
                genes[chrom].append(Gene(start, end, attr["ID"], attr.get("Name", "."), strand))
                counts["retained_gene_feature_rows"] += 1
            except (ValueError, IndexError) as exc:
                reason = str(exc)
                if reason.startswith("invalid literal"):
                    reason = "malformed_integer"
                counts["rejected_" + reason] += 1
                rejects.writerow([lineno, reason, line.rstrip("\r\n")])
    return {chrom: GeneIndex(items) for chrom, items in genes.items()}, counts


def read_candidates(path, lengths):
    with open_text(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 8:
                raise ValueError("TEvoX output requires exactly 8 columns at line %d" % lineno)
            start, end = int(fields[2]) - 1, int(fields[3])
            validate_interval(fields[1], start, end, fields[4], lengths)
            yield fields, start, end


def gene_arrays(genes):
    records = sorted({(g.identifier, g.name) for g in genes})
    return [item[0] for item in records], [item[1] for item in records]


def annotate(args):
    protect_inputs([args.candidates, args.genes, args.metadata, args.lengths], [args.output, str(args.output) + ".qc.json", str(args.output) + ".genes_rejected.tsv"])
    lengths = read_lengths(args.lengths)
    keys = {fields[7] for fields, _, _ in read_candidates(args.candidates, lengths)}
    metadata = {}
    with open_text(args.metadata) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"fragment_id", "chrom", "start0", "end0", "repeat_name", "canonical_class", "family", "classification_note"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError("Not a normalize_annotations metadata TSV")
        for row in reader:
            key = row["fragment_id"]
            if key in keys:
                if key in metadata:
                    raise ValueError("Duplicate fragment_id in metadata: " + key)
                metadata[key] = row
    if keys - metadata.keys():
        raise ValueError("Candidate Name missing from metadata: " + next(iter(keys - metadata.keys())))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    indices, gene_counts = load_genes(args.genes, lengths, str(output) + ".genes_rejected.tsv")
    counts = collections.Counter()
    fields = ["tevox_id", "chrom", "start1", "end1", "strand", "tevox_type", "tevox_family", "fragment_id", "repeat_name", "canonical_class", "repeat_family", "classification_note", "overlap_gene_ids_json", "overlap_gene_names_json", "nearest_gene_ids_json", "nearest_gene_names_json", "nearest_gene_body_distance_bp", "gene_annotation_status"]
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for record, start, end in read_candidates(args.candidates, lengths):
            meta = metadata[record[7]]
            if (meta["chrom"], int(meta["start0"]), int(meta["end0"])) != (record[1], start, end):
                raise ValueError("Candidate/metadata interval mismatch for " + record[7])
            index = indices.get(record[1])
            overlaps, nearest, distance = index.query(start, end) if index else ([], [], None)
            overlap_ids, overlap_names = gene_arrays(overlaps)
            nearest_ids, nearest_names = gene_arrays(nearest)
            status = "overlap" if overlaps else ("nearest" if nearest else "no_gene_feature_on_sequence")
            counts["candidates"] += 1
            counts[status] += 1
            writer.writerow(record + [meta["repeat_name"], meta["canonical_class"], meta["family"], meta["classification_note"]] + [json.dumps(value, ensure_ascii=True) for value in [overlap_ids, overlap_names, nearest_ids, nearest_names]] + [distance if distance is not None else ".", status])
    warnings = []
    if not gene_counts["retained_gene_feature_rows"]:
        warnings.append("No exact GFF3 gene features retained; no transcript-to-gene fallback was invented.")
    if any(k.startswith("rejected_") for k in gene_counts):
        warnings.append("Some gene rows failed validation; inspect genes_rejected.tsv before interpreting gene neighborhoods.")
    qc = {"candidates": str(args.candidates), "genes": str(args.genes), "metadata": str(args.metadata), "counts": dict(counts), "gene_counts": dict(gene_counts), "warnings": warnings, "candidate_coordinate_system": "one-based-inclusive", "distance_definition": "Minimum gene-body base-coordinate difference; overlapping intervals=0; adjacent intervals=1; all nearest ties reported.", "interpretation": "Gene proximity is descriptive; annotation fragments are not confirmed TE insertions or expression/causal effects."}
    with open(str(output) + ".qc.json", "w") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return qc


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", required=True, help="Native TEvoX *_genomeN_unique.txt (Start/End one-based inclusive)")
    p.add_argument("--genes", required=True, help="Matching assembly GFF3 or GFF3.gz; only exact gene feature rows")
    p.add_argument("--metadata", required=True, help="Normalizer metadata.tsv, joined to TEvoX Name")
    p.add_argument("--lengths", required=True)
    p.add_argument("--output", required=True, help="TSV path; adds .qc.json and .genes_rejected.tsv")
    args = p.parse_args()
    try:
        qc = annotate(args)
    except (OSError, ValueError) as exc:
        p.exit(2, "Error: %s\n" % exc)
    print(json.dumps(qc, sort_keys=True))
    return 2 if any(k.startswith("rejected_") for k in qc["gene_counts"]) else 0


if __name__ == "__main__":
    sys.exit(main())
