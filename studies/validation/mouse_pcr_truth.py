#!/usr/bin/env python3
"""Transfer published CAST PCR deletion labels using an exact reference allele.

Input candidates are the predeclared Table S2/S3 join from PMID 37228752:
DEL, nonempty TE_TYPE, CAST PCR 0/1. Only exact, globally unique GRCm39
REF plus 50 bp flanks in the selected C57BL/6NJ assembly can transfer.
C57BL/6J PCR is never substituted for a C57BL/6NJ experimental label.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

COMPLEMENT = bytes.maketrans(b'ACGT', b'TGCA')


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def fasta_records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rb') as handle:
        name, chunks = None, []
        for line in handle:
            if line.startswith(b'>'):
                if name is not None:
                    yield name, b''.join(chunks).upper()
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, b''.join(chunks).upper()


def prepare(args):
    source = json.loads(args.candidates.read_text())
    queries = {}
    for record in source['candidates']:
        sv, reference = record['sv'], record['reference']
        if sv['SVTYPE'] != 'DEL' or not sv['TE_TYPE'] or record['pcr']['CAST'] not in (0, 1):
            raise ValueError('ineligible PCR candidate')
        path = Path(reference['path'])
        if not path.is_absolute():
            path = args.candidates.parent / path
        elif not path.exists():
            path = args.candidates.parent / path.name
        if sha256(path) != reference['sha256']:
            raise ValueError('reference download changed')
        data = json.loads(path.read_text())
        sequence = data['dna'].upper().encode()
        if (reference['flank'] != 50 or reference['start'] != sv['POS'] - 50 or
                reference['end'] != sv['END'] + 50 or len(sequence) != sv['SVLEN'] + 100 or
                data['genome'] != 'mm39' or data['chrom'] != sv['CHROM'] or
                data['start'] != reference['start'] or data['end'] != reference['end'] or
                set(sequence) - set(b'ACGT')):
            raise ValueError('reference allele/coordinate mismatch')
        queries[sv['ID']] = (sequence, sequence.translate(COMPLEMENT)[::-1], record)
    hits = collections.defaultdict(list)
    for chrom, sequence in fasta_records(args.fasta):
        for key, (forward, reverse, record) in queries.items():
            for query, strand in [(forward, '+'), (reverse, '-')]:
                start = sequence.find(query)
                while start >= 0:
                    hits[key].append((chrom, start + 50, start + len(query) - 50, strand))
                    if len(hits[key]) > 1:
                        break
                    start = sequence.find(query, start + 1)
    mapped = collections.defaultdict(list)
    counts = collections.Counter(candidates=len(queries))
    for key in queries:
        count = len(hits[key])
        counts['unique_reference_transfer' if count == 1 else 'no_reference_transfer' if count == 0 else 'ambiguous_reference_transfer'] += 1
        if count == 1:
            chrom, start, end, strand = hits[key][0]
            mapped[chrom].append((start, end, strand, key))
    truth = []
    with open(args.te) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            chrom, start, end, key, score, strand, family, kind = line.rstrip().split('\t')
            start, end = int(start), int(end)
            matches = []
            for left, right, orientation, pcr_id in mapped.get(chrom, []):
                overlap = max(0, min(end, right) - max(start, left))
                if overlap / (end - start) >= .95 and overlap / (right - left) >= .8:
                    matches.append((left, right, orientation, pcr_id))
            if len(matches) != 1:
                continue
            left, right, orientation, pcr_id = matches[0]
            record = queries[pcr_id][2]
            truth.append({'pcr_id': pcr_id, 'source_genome_id': 'c57', 'source_te_id': key,
                          'target_genome_id': 'cast', 'chrom': chrom, 'te_start': start, 'te_end': end,
                          'deletion_start': left, 'deletion_end': right, 'mapping_strand': orientation,
                          'CAST_pcr': record['pcr']['CAST'],
                          'expected_biological_state': 'EMPTY' if record['pcr']['CAST'] == 1 else 'PRESENT',
                          'pcr_sheet': record['pcr']['sheet'], 'external_te_type': record['sv']['TE_TYPE']})
    counts['matched_te_fragments'] = len(truth)
    output = {'status': 'EXTERNAL_PCR_STATE_SUBSET', 'counts': dict(counts), 'truth': truth,
              'transfer_hits': dict(hits), 'protocol': source['protocol'],
              'input_sha256': {name: sha256(getattr(args, name)) for name in ['candidates', 'fasta', 'te']},
              'limitations': ['Small selected PCR SV subset; not a random genome-wide TE sample.',
                              'Reference allele in C57BL/6NJ is supported by exact assembly sequence, not C57BL/6NJ PCR.',
                              'CAST PCR establishes SV state, not fragment-to-fragment candidate membership.',
                              'No transfer to other strains, assembly releases or subgenomes is implied.']}
    args.output.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(dict(counts)))


def evaluate(args):
    truth = json.loads(args.truth.read_text())['truth']
    keys = {(r['source_genome_id'], r['source_te_id'], r['target_genome_id']): r for r in truth}
    if len(keys) != len(truth):
        raise ValueError('duplicate truth fragment')
    seen = {}
    opener = gzip.open if str(args.decisions).endswith('.gz') else open
    with opener(args.decisions, 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = (row['source_genome_id'], row['source_te_id'], row['target_genome_id'])
            if key in keys:
                if key in seen:
                    raise ValueError('duplicate decision')
                seen[key] = row
    confusion = collections.Counter()
    records = []
    for key, label in keys.items():
        row = seen.get(key)
        predicted = row['biological_state'] if row and row['claimable'] == 'true' else 'UNRESOLVED'
        confusion[label['expected_biological_state'] + '->' + predicted] += 1
        records.append(dict(label, predicted_biological_state=predicted,
                            legacy_state=row['legacy_state'] if row else 'NO_DECISION',
                            claimability_reason=row.get('claimability_reason') if row else None))
    report = {'truth_fragments': len(truth), 'matched_decisions': len(seen), 'confusion': dict(confusion),
              'records': records, 'truth_sha256': sha256(args.truth), 'decisions_sha256': sha256(args.decisions),
              'calibration_status': 'NOT_MEMBERSHIP_TRUTH', 'genome_wide_precision': None}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['confusion']))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    a = sub.add_parser('prepare')
    for name in ['candidates', 'fasta', 'te', 'output']:
        a.add_argument('--' + name, required=True, type=Path)
    a = sub.add_parser('evaluate')
    for name in ['truth', 'decisions', 'output']:
        a.add_argument('--' + name, required=True, type=Path)
    args = p.parse_args()
    (prepare if args.command == 'prepare' else evaluate)(args)
