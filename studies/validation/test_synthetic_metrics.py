import gzip
from pathlib import Path
import tempfile
import unittest

from synthetic_suite import evaluate


class MembershipUniverseTests(unittest.TestCase):
    def evaluate_rows(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'run'
            with gzip.open(str(prefix) + '.instances.tsv.gz', 'wt') as handle:
                handle.write('genome_id\tlocus_id\tmember_ids\n' + rows)
            return evaluate(prefix, [('A', 'a', 'L'), ('B', 'b', 'L')])

    def test_duplicate_membership_cannot_be_overwritten(self):
        with self.assertRaisesRegex(ValueError, 'more than once'):
            self.evaluate_rows('A\twrong\ta\nA\tcorrect\ta\nB\tcorrect\tb\n')

    def test_missing_member_stays_in_denominator_contract(self):
        with self.assertRaisesRegex(ValueError, 'universe mismatch'):
            self.evaluate_rows('A\tL\ta\n')

    def test_exact_reconstruction(self):
        result = self.evaluate_rows('A\tpred\ta\nB\tpred\tb\n')
        self.assertEqual(result['b_cubed_f1'], 1)


if __name__ == '__main__':
    unittest.main()
