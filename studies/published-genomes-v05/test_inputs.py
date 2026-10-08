import gzip
from pathlib import Path
import tempfile
import unittest
import normalize_annotations as norm
from paf_exact import Fasta, exact_record


class InputContractTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
    def fasta(self,name,sequence):
        p=self.root/(name+'.fa');p.write_text('>chr1\n'+sequence+'\n')
        f=Fasta(p);self.addCleanup(f.close);return f
    def test_exact_forward_and_ambiguous_bases(self):
        q=self.fasta('A','AACNGT');t=self.fasta('B','AATNGT')
        fields='chr1 6 0 6 + chr1 6 0 6 5 6 60 cg:Z:6M'.split()
        got=exact_record(fields,q,t)
        self.assertEqual(got[-1],'cg:Z:2=2X2=')
        self.assertEqual(got[9:11],['4','6'])
    def test_negative_insertion_and_deletion_order(self):
        # Query forward AACCAGT -> traversal ACTGGTT, target ACATT.
        q=self.fasta('A','AACCAGT');t=self.fasta('B','ACATT')
        fields='chr1 7 0 7 - chr1 5 0 5 4 7 60 cg:Z:2M2I3M'.split()
        got=exact_record(fields,q,t)
        self.assertEqual(got[-1],'cg:Z:2=2I1X2=')
        self.assertEqual(got[9:11],['4','7'])
    def test_invalid_consumption_and_bounds_rejected(self):
        q=self.fasta('A','AAAA');t=self.fasta('B','AAAA')
        for text in ['chr1 4 0 4 + chr1 4 0 4 4 4 60 cg:Z:3M',
                     'chr1 4 -1 4 + chr1 4 0 4 4 4 60 cg:Z:5M']:
            with self.assertRaises(ValueError):exact_record(text.split(),q,t)
    def test_bed8_keeps_family_and_includes_sva(self):
        lengths=self.root/'lengths';lengths.write_text('chr1\t1000\n')
        source=self.root/'repeat.bed';source.write_text('chr1\t10\t40\tSVA_F\tSVA\t+\n')
        prefix=self.root/'out'
        args=norm.parser().parse_args(['--input',str(source),'--format','te-bed','--lengths',str(lengths),
            '--output-prefix',str(prefix),'--id-prefix','A','--bed-coordinates','zero-based-half-open',
            '--bed-class-column','5','--bed-strand-column','6'])
        qc=norm.normalize(args)
        fields=Path(str(prefix)+'.bed').read_text().strip().split('\t')
        self.assertEqual(len(fields),8);self.assertEqual(fields[6:],['SVA_F','SVA'])
        self.assertEqual(qc['counts']['retained_fragments'],1)


if __name__=='__main__':unittest.main()
