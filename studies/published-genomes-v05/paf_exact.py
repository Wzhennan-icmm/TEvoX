"""Recover =/X from an existing cg alignment and the exact two FASTAs.

Alignment coordinates, strand, MAPQ and the indel path stay fixed. This is
sequence-based annotation of an existing alignment, not independent alignment
evidence. Non-ACGT comparisons are X, even for two equal ambiguous letters.
"""
import gzip
import hashlib
import json
import mmap
import re
from pathlib import Path

import numpy as np

COMPLEMENT = bytes.maketrans(b'ACGTRYMKBDHVNacgtrymkbdhvn', b'TGCAYRKMVHDBNtgcayrkmvhdbn')
CIGAR = re.compile(r'([1-9][0-9]*)([MID=X])')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


class Fasta:
    def __init__(self, path):
        self.path = Path(path)
        self.handle = self.path.open('rb')
        self.sequences = {}
        current = None
        short = False
        for line in self.handle:
            if line.startswith(b'>'):
                current = line[1:].split()[0].decode('ascii')
                if current in self.sequences:
                    raise ValueError('duplicate FASTA sequence: ' + current)
                self.sequences[current] = [0, self.handle.tell(), 0, 0]
                short = False
            else:
                bases = line.rstrip(b'\r\n')
                if current is None or not bases or short:
                    raise ValueError('invalid FASTA wrapping')
                row = self.sequences[current]
                if row[2] == 0:
                    row[2:] = [len(bases), len(line)]
                if len(bases) > row[2] or (len(bases) == row[2] and len(line) != row[3]):
                    raise ValueError('inconsistent FASTA line width')
                short = len(bases) < row[2]
                row[0] += len(bases)
        if not self.sequences or any(r[0] == 0 for r in self.sequences.values()):
            raise ValueError('empty FASTA')
        self.data = mmap.mmap(self.handle.fileno(), 0, access=mmap.ACCESS_READ)

    def fetch(self, name, start, end):
        length, offset, width, byte_width = self.sequences[name]
        if not 0 <= start <= end <= length:
            raise ValueError('FASTA interval outside sequence')
        if start == end:
            return b''
        begin = offset + start // width * byte_width + start % width
        stop = offset + (end - 1) // width * byte_width + (end - 1) % width + 1
        seq = self.data[begin:stop].replace(b'\n', b'').replace(b'\r', b'').upper()
        if len(seq) != end - start:
            raise ValueError('FASTA index mismatch')
        return seq

    def close(self):
        self.data.close()
        self.handle.close()


def exact_record(fields, query, target):
    if len(fields) < 12 or fields[4] not in {'+', '-'}:
        raise ValueError('invalid PAF')
    qname, tname = fields[0], fields[5]
    qlen, qs, qe = map(int, fields[1:4])
    tlen, ts, te = map(int, fields[6:9])
    if not (0 <= qs < qe <= qlen and 0 <= ts < te <= tlen):
        raise ValueError('invalid PAF bounds')
    if qlen != query.sequences[qname][0] or tlen != target.sequences[tname][0]:
        raise ValueError('PAF/FASTA length mismatch')
    tags = [v for v in fields[12:] if v.startswith('cg:Z:')]
    if len(tags) != 1:
        raise ValueError('exactly one cg tag required')
    cigar = tags[0][5:]
    operations = CIGAR.findall(cigar)
    if ''.join(n + op for n, op in operations) != cigar:
        raise ValueError('unsupported/malformed cg')
    qp, tp = (qs if fields[4] == '+' else qe), ts
    result = []
    matches = block_length = 0

    def append(n, op):
        if result and result[-1][1] == op:
            result[-1][0] += n
        else:
            result.append([n, op])

    for number, op in operations:
        n = int(number)
        block_length += n
        if op in 'M=X':
            qseq = query.fetch(qname, qp, qp+n) if fields[4] == '+' else query.fetch(qname, qp-n, qp).translate(COMPLEMENT)[::-1]
            tseq = target.fetch(tname, tp, tp+n)
            a = np.frombuffer(qseq, dtype=np.uint8)
            b = np.frombuffer(tseq, dtype=np.uint8)
            same = (a == b) & ((a == 65) | (a == 67) | (a == 71) | (a == 84))
            cuts = np.flatnonzero(same[1:] != same[:-1]) + 1
            begin = 0
            for stop in [*cuts, n]:
                size = int(stop) - begin
                code = '=' if same[begin] else 'X'
                append(size, code)
                if code == '=':
                    matches += size
                begin = int(stop)
            qp += n if fields[4] == '+' else -n
            tp += n
        elif op == 'I':
            append(n, op)
            qp += n if fields[4] == '+' else -n
        else:
            append(n, op)
            tp += n
    if qp != (qe if fields[4] == '+' else qs) or tp != te:
        raise ValueError('CIGAR consumption disagrees with PAF endpoints')
    # Drop alignment-score tags after re-annotation instead of leaving stale NM.
    # PAF core counts reflect the reconstructed exact path, including ambiguous X.
    core = fields[:12]
    core[9:11] = [str(matches), str(block_length)]
    return core + ['cg:Z:' + ''.join(str(n)+op for n, op in result)]


def convert(source, destination, query, target):
    source, destination = Path(source), Path(destination)
    if destination.exists() or destination.resolve() in {source.resolve(), query.path.resolve(), target.path.resolve()}:
        raise ValueError('output exists or aliases an input')
    opener = gzip.open if source.suffix == '.gz' else open
    temporary = destination.with_suffix(destination.suffix + '.part')
    count = 0
    with opener(source, 'rt') as inp, temporary.open('w') as out:
        for line in inp:
            if not line.strip() or line.startswith('#'):
                continue
            fields = exact_record(line.rstrip('\r\n').split('\t'), query, target)
            out.write('\t'.join(fields) + '\n')
            count += 1
    if count == 0:
        raise ValueError('empty alignment')
    temporary.replace(destination)
    return {'records': count, 'source_sha256': sha256(source), 'output_sha256': sha256(destination),
            'method': 'Exact =/X reconstructed from the unchanged PAF indel path and both FASTAs; no independent reciprocal support added',
            'ambiguous_bases': 'non-ACGT treated as mismatches', 'numpy_version': np.__version__}
