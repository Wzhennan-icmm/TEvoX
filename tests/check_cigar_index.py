#!/usr/bin/env python3
"""Compare indexed projection to the retained linear traversal on dense CIGARs."""
import argparse, random, subprocess, tempfile
from pathlib import Path

def check(indexed,linear):
    with tempfile.TemporaryDirectory(prefix='tevox-cigar-index-') as directory:
        root=Path(directory);rng=random.Random(310521)
        for strand in ['+','-']:
            ops=[(rng.randint(1,11),rng.choice(['=','=','=','X','I','D'])) for _ in range(2000)]
            qlen=sum(n for n,op in ops if op!='D');tlen=sum(n for n,op in ops if op!='I')
            (root/'A.fa').write_text('>chr1\n'+'A'*qlen+'\n')
            (root/'B.fa').write_text('>chr1\n'+'A'*tlen+'\n')
            for sample,length in [('A',qlen),('B',tlen)]:
                # Include TE boundaries around checkpoint-adjacent indels as well
                # as long intervals spanning many checkpoint blocks.
                intervals=[(n,min(n+rng.randint(1,180),length)) for n in range(0,length,13)]
                (root/(sample+'.bed')).write_text(''.join('chr1\t%d\t%d\t%s%d\t0\t+\tFam1\tLINE\n'%(a,b,sample,i) for i,(a,b) in enumerate(intervals)))
            cigar=''.join(str(n)+op for n,op in ops)
            (root/'input.paf').write_text('\t'.join(map(str,['chr1',qlen,0,qlen,strand,'chr1',tlen,0,tlen,sum(n for n,op in ops if op=='='),sum(n for n,op in ops),60]))+'\tcg:Z:'+cigar+'\n')
            common=['pair','--genome-a','A','--fasta-a',str(root/'A.fa'),'--te-a',str(root/'A.bed'),
                    '--genome-b','B','--fasta-b',str(root/'B.fa'),'--te-b',str(root/'B.bed'),
                    '--paf',str(root/'input.paf'),'--flank','17','--candidate-window','5']
            for binary,label in [(indexed,'indexed'),(linear,'linear')]:
                subprocess.run([binary,*common,'--output',str(root/label)],check=True,capture_output=True)
            for output in root.glob('indexed.*.tsv'):
                other=root/('linear'+output.name[len('indexed'):])
                if output.read_bytes()!=other.read_bytes():raise AssertionError('different '+output.name+' strand '+strand)
        print('Dense CIGAR equivalence: both strands, all 17 TSV outputs identical')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('indexed');p.add_argument('linear');a=p.parse_args();check(a.indexed,a.linear)
