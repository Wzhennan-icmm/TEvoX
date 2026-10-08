import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
BINARY=Path(os.environ.get('TEVOX_BIN',str(HERE.parents[1]/'tevox'))).resolve()


class RunIntegrityTests(unittest.TestCase):
    def setUp(self):
        directory=tempfile.TemporaryDirectory();self.addCleanup(directory.cleanup)
        self.root=Path(directory.name);genomes=[]
        for name in ['A','B']:
            (self.root/(name+'.fa')).write_text('>chr1\n'+'A'*1000+'\n')
            (self.root/(name+'.bed')).write_text('chr1\t400\t440\tFAM\tDNA\t+\n')
            (self.root/(name+'.gff3')).write_text('##gff-version 3\nchr1\ttest\tgene\t701\t800\t.\t+\t.\tID='+name+'_gene\n')
            genomes.append({'id':name,'fasta':name+'.fa','repeats':name+'.bed','genes':name+'.gff3','repeat_format':'te-bed',
                'normalize_args':['--bed-coordinates','zero-based-half-open','--bed-class-column','5','--bed-strand-column','6']})
        (self.root/'pair.json').write_text(json.dumps({'genomes':genomes}))
        (self.root/'pair.paf').write_text('chr1\t1000\t0\t1000\t+\tchr1\t1000\t0\t1000\t1000\t1000\t60\tcg:Z:1000M\n')
        self.command=[sys.executable,str(HERE/'run_pair.py'),'--manifest',str(self.root/'pair.json'),'--paf',str(self.root/'pair.paf'),
                      '--tevox',str(BINARY),'--output',str(self.root/'run')]
        result=subprocess.run(self.command,cwd='/',capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        complete=json.loads((self.root/'run/complete.json').read_text())
        self.assertEqual(len(complete['signature']['source_inputs']),6)
        (self.root/'run/complete.json').unlink() # simulate interruption after analysis

    def test_changed_source_alignment_is_rejected_before_cache_reuse(self):
        with (self.root/'pair.paf').open('a') as handle:handle.write('# changed source\n')
        result=subprocess.run(self.command,capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('source/configuration changed',result.stderr)

    def test_changed_prepared_fasta_is_rejected(self):
        fasta=self.root/'run/A/genome.fa';fasta.write_text(fasta.read_text().replace('AAAA','CAAA',1))
        result=subprocess.run(self.command,capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('unverified or changed prepared FASTA',result.stderr)

    def test_changed_repeat_source_is_rejected(self):
        with (self.root/'A.bed').open('a') as handle:handle.write('chr1\t500\t520\tFAM\tDNA\t+\n')
        result=subprocess.run(self.command,capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('source/configuration changed',result.stderr)


if __name__=='__main__':unittest.main()
