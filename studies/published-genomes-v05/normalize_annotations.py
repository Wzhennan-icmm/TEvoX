#!/usr/bin/env python3
"""Normalize explicit repeat formats to TEvoX v0.5 BED8 and an identity-preserving sidecar.

BED col4 is a unique annotation-fragment key (TEvoX preserves it as Name); col7
is comparison family and col8 is TE class. Fragments are never merged into presumed biological insertions.
All output intervals are zero-based half-open. Original GFF3/RM coordinates are
one-based inclusive. Never infer coordinate convention from a filename.
"""
import argparse
import collections
import csv
import gzip
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote

NON_TE = {"simple_repeat", "low_complexity", "satellite", "rrna", "trna", "snrna", "scrna", "srprna", "rna"}
CANONICAL = {"line": "LINE", "sva": "SVA", "l1-dep": "LINE", "sine": "SINE", "ltr": "LTR", "dna": "DNA", "tir": "DNA", "mite": "DNA", "rc": "RC", "helitron": "RC", "unknown": "Unknown"}
META_FIELDS = ["fragment_id", "chrom", "start0", "end0", "strand", "repeat_name", "canonical_class", "family", "original_class", "source_id", "source_line", "classification_note"]


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def protect_inputs(inputs, outputs):
    source_paths = {Path(path).resolve() for path in inputs if path}
    for output in outputs:
        if Path(output).resolve() in source_paths:
            raise ValueError("Output would overwrite an input: " + str(output))


def attributes(text):
    """GFF3 splitting must precede percent decoding (%3B is not a delimiter)."""
    result = {}
    for item in text.split(";"):
        if not item or item == ".":
            continue
        if "=" not in item:
            raise ValueError("malformed_gff_attribute")
        key, value = item.split("=", 1)
        key = unquote(key.strip())
        if key in result:
            raise ValueError("duplicate_gff_attribute")
        result[key] = unquote(value)
    return result


def read_lengths(path):
    with open_text(path) as handle:
        if (str(path)[:-3] if str(path).endswith(".gz") else str(path)).endswith(".json"):
            lengths = json.load(handle)
            if not isinstance(lengths, dict):
                raise ValueError("Length JSON must be a {sequence: length} object")
        else:
            lengths = {}
            for lineno, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 2:
                    raise ValueError("Invalid lengths TSV line %d" % lineno)
                name, value = fields[:2]
                if name in lengths:
                    raise ValueError("Duplicate sequence in lengths: " + name)
                lengths[name] = int(value)
    if not lengths:
        raise ValueError("No FASTA sequence lengths supplied")
    for name, value in lengths.items():
        if not isinstance(name, str) or not name or any(c.isspace() for c in name):
            raise ValueError("Invalid FASTA sequence ID")
        if type(value) is not int or value <= 0:
            raise ValueError("Invalid sequence length for " + name)
    return lengths


def validate_interval(chrom, start, end, strand, lengths):
    if chrom not in lengths:
        raise ValueError("sequence_not_in_fasta")
    if start < 0 or end <= start or end > lengths[chrom]:
        raise ValueError("interval_out_of_bounds")
    if end > 2147483647:
        raise ValueError("interval_exceeds_tevox_signed_int")
    if strand not in {"+", "-", "."}:
        raise ValueError("invalid_strand")


def classification(raw):
    """Only explicit class values are interpreted, never repeat names/features."""
    raw = raw.strip()
    if not raw or raw == ".":
        raise ValueError("unresolved_class")
    base, _, family = raw.partition("/")
    key = base.lower()
    if key in NON_TE:
        raise ValueError("excluded_non_te")
    if key not in CANONICAL:
        raise ValueError("unresolved_class")
    canonical = CANONICAL[key]
    if key == "helitron" and not family:
        family = "Helitron"
    note = "unknown_class_not_confirmed_te" if canonical == "Unknown" else "explicit_source_class"
    return canonical, family or ".", note


def read_mapping(path):
    if not path:
        return {}
    result = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"repeat_name", "class"} <= set(reader.fieldnames):
            raise ValueError("Class mapping needs repeat_name and class TSV headers")
        for row in reader:
            key = row["repeat_name"]
            value = row["class"] + ("/" + row["family"] if row.get("family", ".") not in {"", "."} and "/" not in row["class"] else "")
            if key in result and result[key] != value:
                raise ValueError("Conflicting class mapping for " + key)
            result[key] = value
    return result


def parse_record(line, fmt, args):
    """Return chrom,start0,end0,strand,name,class,source_id,optional_family."""
    if fmt == "repeatmasker-out":
        fields = line.split()
        if not fields or not fields[0].isdigit():
            # RepeatMasker reports have a three-line header and optional no-repeat note.
            if not fields or fields[0] in {"SW", "score", "There"}:
                return None
            raise ValueError("malformed_repeatmasker_row")
        if len(fields) < 14:
            raise ValueError("insufficient_columns")
        return fields[4], int(fields[5]) - 1, int(fields[6]), "-" if fields[8] == "C" else fields[8], fields[9], fields[10], fields[14] if len(fields) > 14 and fields[14] != "*" else ".", "."
    fields = line.rstrip("\r\n").split("\t")
    if fmt == "ucsc-rmsk":
        # Standard rmsk table: 17 columns with bin, 16 without. bigRmsk is different.
        offset = 1 if len(fields) == 17 else 0
        if len(fields) not in {16, 17}:
            raise ValueError("rmsk_requires_16_or_17_columns_not_bigrmsk")
        f = fields[offset:]
        return f[4], int(f[5]), int(f[6]), f[8], f[9], f[10], f[15], f[11]
    if fmt == "repeat-gff3":
        if len(fields) != 9:
            raise ValueError("gff_requires_nine_columns")
        if fields[2] not in args.gff_feature:
            raise ValueError("excluded_feature_type")
        attr = attributes(fields[8])
        return fields[0], int(fields[3]) - 1, int(fields[4]), "." if fields[6] == "?" else fields[6], attr.get(args.name_attribute, attr.get("ID", ".")), attr.get(args.class_attribute, "."), attr.get("ID", "."), attr.get(args.family_attribute, ".") if args.family_attribute else "."
    if fmt == "te-bed":
        needed = max(3, args.bed_name_column or 0, args.bed_class_column or 0, args.bed_strand_column or 0, args.bed_family_column or 0, args.bed_id_column or 0)
        if len(fields) < needed:
            raise ValueError("insufficient_columns")
        start = int(fields[1]) - (1 if args.bed_coordinates == "one-based-inclusive" else 0)
        return fields[0], start, int(fields[2]), fields[args.bed_strand_column - 1] if args.bed_strand_column else ".", fields[args.bed_name_column - 1] if args.bed_name_column else ".", fields[args.bed_class_column - 1] if args.bed_class_column else ".", fields[args.bed_id_column - 1] if args.bed_id_column else ".", fields[args.bed_family_column - 1] if args.bed_family_column else "."
    raise ValueError("unsupported_format")


def input_records(line, args):
    if args.format != "ucsc-bigrmsk-bed":
        yield None, parse_record(line, args.format, args)
        return
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != 14:
        raise ValueError("bigrmsk_requires_fourteen_columns")
    aligned_start, aligned_end = int(fields[6]), int(fields[7])
    if aligned_start < 0 or aligned_end <= aligned_start:
        raise ValueError("invalid_bigrmsk_aligned_bounds")
    descriptions = [item.strip() for item in fields[13].split(",") if item.strip()]
    if not descriptions:
        raise ValueError("missing_bigrmsk_repeatmasker_descriptions")
    # Validate the entire joined record before yielding any retained fragments.
    records = []
    for description in descriptions:
        record = parse_record(description, "repeatmasker-out", args)
        if record is None:
            raise ValueError("invalid_bigrmsk_repeatmasker_description")
        chrom, start, end, strand = record[:4]
        if chrom != fields[0]:
            raise ValueError("bigrmsk_description_chromosome_mismatch")
        if strand != fields[5] or start < aligned_start or end > aligned_end:
            raise ValueError("bigrmsk_description_alignment_mismatch")
        records.append(record)
    for part, record in enumerate(records, 1):
        yield part, record


def normalize(args):
    lengths = read_lengths(args.lengths)
    mapping = read_mapping(args.class_map)
    counts, classes = collections.Counter(), collections.Counter()
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = {kind: str(prefix) + suffix for kind, suffix in [("bed", ".bed"), ("metadata", ".metadata.tsv"), ("rejected", ".rejected.tsv"), ("qc", ".qc.json")]}
    protect_inputs([args.input, args.lengths, args.class_map], paths.values())
    with open_text(args.input) as source, open(paths["bed"], "w") as bed, open(paths["metadata"], "w") as metadata, open(paths["rejected"], "w") as rejected:
        meta = csv.DictWriter(metadata, fieldnames=META_FIELDS, delimiter="\t", lineterminator="\n")
        meta.writeheader()
        rejects = csv.writer(rejected, delimiter="\t", lineterminator="\n")
        rejects.writerow(["source_line", "reason", "source_record"])
        for lineno, line in enumerate(source, 1):
            if line.startswith("##FASTA"):
                counts["embedded_fasta_not_processed"] += 1
                break
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            counts["source_records"] += 1
            try:
                records = list(input_records(line, args))
            except (ValueError, IndexError) as exc:
                reason = str(exc)
                if reason.startswith("invalid literal") or isinstance(exc, IndexError):
                    reason = "malformed_integer_or_columns"
                counts["rejected_" + reason] += 1
                counts["malformed_source_records"] += 1
                rejects.writerow([lineno, reason, line.rstrip("\r\n")])
                continue
            source_descriptions = [item.strip() for item in line.rstrip("\r\n").split("\t")[13].split(",") if item.strip()] if args.format == "ucsc-bigrmsk-bed" else None
            for part, record in records:
                if record is None:
                    counts["header_records"] += 1
                    continue
                counts["parsed_fragments"] += 1
                try:
                    chrom, start, end, strand, name, rawclass, source_id, source_family = record
                    validate_interval(chrom, start, end, strand, lengths)
                    original_class = rawclass
                    if name in mapping:
                        rawclass = mapping[name]
                    canonical, family, note = classification(rawclass)
                    if name in mapping:
                        note += ";explicit_class_map"
                    elif source_family not in {"", "."}:
                        family = source_family
                    fragment_id = "%s_F%012d" % (args.id_prefix, lineno)
                    if part is not None:
                        fragment_id += "_P%06d" % part
                    comparison_family = (name if name not in {"", "."} else ".") if args.family_policy == "repeat-name" else canonical + ("/" + family if family != "." else "")
                    if canonical == "Unknown":
                        comparison_family = "."
                    if any(c in comparison_family for c in "\t\r\n"):
                        raise ValueError("invalid_family_delimiter")
                    bed.write("\t".join([chrom, str(start), str(end), fragment_id, "0", strand, comparison_family, canonical]) + "\n")
                    meta.writerow(dict(zip(META_FIELDS, [fragment_id, chrom, start, end, strand, name, canonical, family, original_class, source_id, lineno, note])))
                    counts["retained_fragments"] += 1
                    classes[canonical] += 1
                    if canonical == "Unknown":
                        counts["unknown_class_fragments"] += 1
                except (ValueError, IndexError) as exc:
                    reason = str(exc)
                    if reason.startswith("invalid literal") or isinstance(exc, IndexError):
                        reason = "malformed_integer_or_columns"
                    counts["rejected_" + reason] += 1
                    source_line = str(lineno) + (".%d" % part if part is not None else "")
                    raw_record = source_descriptions[part - 1] if source_descriptions is not None else line.rstrip("\r\n")
                    rejects.writerow([source_line, reason, raw_record])
    qc = {
        "input": str(args.input), "format": args.format, "lengths": str(args.lengths),
        "input_coordinate_system": args.bed_coordinates if args.format == "te-bed" else ("zero-based-half-open" if args.format == "ucsc-rmsk" else ("original_one-based-inclusive_descriptions_in_bigrmsk_bed14" if args.format == "ucsc-bigrmsk-bed" else "one-based-inclusive")),
        "output_coordinate_system": "zero-based-half-open", "counts": dict(counts), "classes": dict(classes),
        "outputs": paths, "interpretation": "Annotation fragments only; not reconstructed elements or proven insertion/deletion events. Unknown class retained explicitly with uncertainty.",
        "configuration": {k: v for k, v in vars(args).items() if k != "input"},
    }
    with open(paths["qc"], "w") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return qc


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--family-policy", choices=["repeat-name", "classification"], default="repeat-name")
    p.add_argument("--input", required=True)
    p.add_argument("--format", required=True, choices=["repeatmasker-out", "ucsc-rmsk", "ucsc-bigrmsk-bed", "repeat-gff3", "te-bed"])
    p.add_argument("--lengths", required=True, help="FASTA-derived TSV (ID,length; FAI accepted), or JSON ID:length")
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--id-prefix", required=True, help="Assembly/material key, ASCII letters/digits/underscore/dot/hyphen")
    p.add_argument("--class-map", help="Explicit repeat_name,class[,family] headered TSV; replaces source class on exact name match")
    p.add_argument("--gff-feature", action="append", default=[], help="Exact feature type to include; repeat to include several (required for GFF)")
    p.add_argument("--name-attribute", default="Name")
    p.add_argument("--class-attribute", default="Classification")
    p.add_argument("--family-attribute")
    p.add_argument("--bed-coordinates", choices=["zero-based-half-open", "one-based-inclusive"], help="Required for te-bed; do not infer from .bed extension")
    p.add_argument("--bed-name-column", type=int, default=4)
    p.add_argument("--bed-class-column", type=int, help="One-based column, or supply --class-map keyed by name")
    p.add_argument("--bed-strand-column", type=int, help="One-based column; omitted means unknown strand")
    p.add_argument("--bed-family-column", type=int, help="One-based source family column")
    p.add_argument("--bed-id-column", type=int, help="One-based original annotation ID column")
    return p


def main():
    p = parser()
    args = p.parse_args()
    if len(args.id_prefix) > 128 or not re.fullmatch(r"[A-Za-z0-9_.-]+", args.id_prefix):
        p.error("--id-prefix must have <=128 ASCII letters/digits/underscore/dot/hyphen")
    if args.format == "repeat-gff3" and not args.gff_feature:
        p.error("repeat-gff3 requires explicit --gff-feature to avoid mixing parent and child annotations")
    if args.format == "te-bed" and (not args.bed_coordinates or not (args.bed_class_column or args.class_map)):
        p.error("te-bed requires --bed-coordinates and either --bed-class-column or --class-map")
    if any(value is not None and value < 1 for value in [args.bed_name_column, args.bed_class_column, args.bed_strand_column, args.bed_family_column, args.bed_id_column]):
        p.error("BED column numbers must be positive")
    try:
        qc = normalize(args)
    except (OSError, ValueError) as exc:
        p.exit(2, "Error: %s\n" % exc)
    print(json.dumps(qc, sort_keys=True))
    counts = qc["counts"]
    invalid = sum(v for k, v in counts.items() if k.startswith("rejected_") and k not in {"rejected_excluded_non_te", "rejected_excluded_feature_type", "rejected_unresolved_class"})
    return 2 if invalid or not counts.get("retained_fragments") else 0


if __name__ == "__main__":
    sys.exit(main())
