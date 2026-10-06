import contextlib
import io
import json
from pathlib import Path
import random
import tempfile
from types import SimpleNamespace
import unittest

from mouse_pcr_truth import COMPLEMENT, prepare, sha256


class PCRTransferTests(unittest.TestCase):
    def run_transfer(self, mode, label=1):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rng = random.Random(91)
            sequence = ''.join(rng.choice('ACGT') for _ in range(200))
            ref = root / 'reference.json'
            ref.write_text(json.dumps({'dna': sequence, 'genome': 'mm39', 'chrom': 'chr1', 'start': 950, 'end': 1150}))
            record = {'sv': {'ID': 'DEL1', 'SVTYPE': 'DEL', 'TE_TYPE': 'L1', 'POS': 1000, 'END': 1100, 'SVLEN': 100, 'CHROM': 'chr1'},
                      'pcr': {'CAST': label, 'sheet': 'PCR_validations'},
                      'reference': {'path': str(ref), 'sha256': sha256(ref), 'flank': 50, 'start': 950, 'end': 1150}}
            candidates = root / 'candidates.json'
            candidates.write_text(json.dumps({'protocol': 'test', 'candidates': [record]}))
            oriented = sequence.encode().translate(COMPLEMENT)[::-1].decode() if mode == 'reverse' else sequence
            if mode == 'mismatch':
                oriented = oriented[:90] + ('A' if oriented[90] != 'A' else 'C') + oriented[91:]
            fasta = root / 'input.fa'
            fasta.write_text('>contigA\n' + 'N' * 10 + oriented + '\n' +
                             ('>contigB\n' + sequence + '\n' if mode == 'duplicate' else ''))
            te = root / 'te.bed'
            te.write_text('contigA\t60\t160\tT1\t0\t+\tL1\tLINE\n')
            output = root / 'truth.json'
            with contextlib.redirect_stdout(io.StringIO()):
                prepare(SimpleNamespace(candidates=candidates, fasta=fasta, te=te, output=output))
            return json.loads(output.read_text())

    def test_reverse_transfer_preserves_deletion_and_pcr_label(self):
        result = self.run_transfer('reverse')
        self.assertEqual(len(result['truth']), 1)
        self.assertEqual(result['truth'][0]['mapping_strand'], '-')
        self.assertEqual(result['truth'][0]['expected_biological_state'], 'EMPTY')
        self.assertEqual(result['truth'][0]['deletion_start'], 60)

    def test_observed_zero_is_presence_not_missing_label(self):
        result = self.run_transfer('forward', 0)
        self.assertEqual(result['truth'][0]['expected_biological_state'], 'PRESENT')

    def test_duplicate_transfer_is_excluded(self):
        result = self.run_transfer('duplicate')
        self.assertEqual(result['truth'], [])
        self.assertEqual(result['counts']['ambiguous_reference_transfer'], 1)

    def test_changed_reference_allele_is_excluded(self):
        result = self.run_transfer('mismatch')
        self.assertEqual(result['truth'], [])
        self.assertEqual(result['counts']['no_reference_transfer'], 1)


if __name__ == '__main__':
    unittest.main()
