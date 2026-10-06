import gzip
from pathlib import Path
import tempfile
import unittest

from giab_empty_truth import Chains, Intervals


class ChainMappingTests(unittest.TestCase):
    def chain(self,text):
        directory=tempfile.TemporaryDirectory();self.addCleanup(directory.cleanup)
        path=Path(directory.name)/'chain.gz'
        with gzip.open(path,'wt') as handle:handle.write(text)
        return Chains(path)

    def test_forward_reverse_and_reciprocal_coordinates(self):
        forward=self.chain('chain 10 chr1 1000 + 0 100 chrX 1000 - 400 500 1\n100\n')
        backward=self.chain('chain 10 chrX 1000 + 500 600 chr1 1000 - 900 1000 2\n100\n')
        self.assertEqual(forward.map('chr1',10,20),[('chrX',580,590,'-')])
        self.assertEqual(backward.map('chrX',580,590),[('chr1',10,20,'-')])
        self.assertEqual(forward.map('chr1',95,105),[])

    def test_gapped_or_nonunique_windows_remain_ineligible(self):
        chain=self.chain('chain 10 chr1 1000 + 0 100 chrX 1000 + 10 110 1\n40 20 20\n40\n'
                         'chain 9 chr1 1000 + 0 40 chrY 1000 + 200 240 2\n40\n')
        self.assertEqual(len(chain.map('chr1',10,20)),2)
        self.assertEqual(chain.map('chr1',30,70),[])
        self.assertEqual(chain.map('chr1',60,100),[('chrX',70,110,'+')])

    def test_interval_index_retains_nested_long_intervals(self):
        index=Intervals([(0,100),(10,20),(30,40),(90,110)])
        self.assertEqual(set(index.overlap(50,60)),{(0,100)})
        self.assertEqual(set(index.overlap(100,101)),{(90,110)})


if __name__=='__main__':unittest.main()
