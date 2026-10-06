#!/usr/bin/env python3
"""Controlled copy-context ablations and sparse scaling; not biological truth.

The generator fixes ancestry labels before running TEvoX. Exact supplied
alignments isolate engine/optimizer behavior from alignment discovery.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def write_tsv(path, fields, rows):
    with path.open('w',newline='') as handle:
        w=csv.writer(handle,delimiter='\t');w.writerow(fields);w.writerows(rows)


def generate(root, loci, copies, seed, missing=False):
    root.mkdir(parents=True,exist_ok=False)
    length=loci*200+200
    rng=random.Random(seed);sequence=''.join(rng.choices('ACGT',k=length))
    truth=[]
    for genome,ncopies in [('O',1),('P',copies)]:
        with (root/(genome+'.fa')).open('w') as fa, (root/(genome+'.bed')).open('w') as bed:
            for copy in range(ncopies):
                chrom=genome+str(copy)
                fa.write('>'+chrom+'\n')
                for i in range(0,length,80):fa.write(sequence[i:i+80]+'\n')
                for locus in range(loci):
                    if missing and genome=='P' and copy==0 and locus%7==0:continue
                    start=100+locus*200;key='{}_L{}'.format(chrom,locus)
                    bed.write('\t'.join(map(str,[chrom,start,start+40,key,0,'+','FAM','DNA']))+'\n')
                    truth.append((genome,key,'ancestral_'+str(locus)))
    write_tsv(root/'truth.members.tsv',['genome_id','te_id','ancestral_locus'],truth)
    # Copy capacity one deliberately tests explicit WGD contexts overriding
    # fallback quota. The no-synteny ablation must retain that same manifest.
    write_tsv(root/'manifest.tsv',['genome_id','fasta','te_annotation','max_locus_copies'],
              [(g,g+'.fa',g+'.bed',1) for g in ['O','P']])
    with (root/'O_P.paf').open('w') as paf:
        for copy in range(copies):
            paf.write('\t'.join(map(str,['O0',length,0,length,'+','P'+str(copy),length,0,length,length,length,60,'cg:Z:'+str(length)+'=']))+'\n')
    write_tsv(root/'alignments.tsv',['query_id','target_id','paf'],[('O','P','O_P.paf')])
    genes=[];anchors=[20,length//2,length-40]
    for genome,ncopies in [('O',1),('P',copies)]:
        for copy in range(ncopies):
            for i,start in enumerate(anchors):
                genes.append([genome,f'{genome}{copy}_G{i}',f'{genome}{copy}',start,start+10,'+',f'{genome}{copy}','hap1'])
    write_tsv(root/'genes.tsv',['genome_id','gene_id','contig','start','end','strand','subgenome_id','haplotype_id'],genes)
    with (root/'blocks.collinearity').open('w') as blocks:
        for copy in range(copies):
            blocks.write(f'## Alignment {copy}: score=300.0 e_value=1e-30 N=3 O0&P{copy} plus\n')
            for i in range(3):blocks.write(f'  {copy}- {i}:\tO0_G{i}\tP{copy}_G{i}\t1e-30\n')
    write_tsv(root/'synteny.tsv',['source_id','format','collinearity','genes','wgd_node'],
              [('controlled','mcscanx','blocks.collinearity','genes.tsv','WGD_SIM')])
    provenance={'truth_origin':'SEEDED_CONTROLLED_GENERATOR_NOT_EXTERNAL_BIOLOGY','seed':seed,
                'loci':loci,'polyploid_subgenome_copies':copies,'missing_annotation':missing,
                'alignment':'exact supplied sequence-identical paths; no alignment search measured',
                'inputs':{p.name:digest(p) for p in sorted(root.iterdir())}}
    (root/'generation.json').write_text(json.dumps(provenance,indent=2)+'\n')
    return truth


def evaluate(prefix,truth):
    # Membership is reconstructed from selected node-to-locus instances.
    expected={(g,t):l for g,t,l in truth};observed={}
    with gzip.open(str(prefix)+'.instances.tsv.gz','rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            for key in row['member_ids'].split(','):
                if key not in {'','.'}:observed[(row['genome_id'],key)]=row['locus_id']
    if set(expected)!=set(observed):raise ValueError('truth/predicted membership universe mismatch')
    intersections=collections.Counter((expected[k],observed[k]) for k in expected)
    true_size=collections.Counter(expected.values());pred_size=collections.Counter(observed.values())
    n=len(expected)
    precision=sum(v*v/pred_size[p] for (t,p),v in intersections.items())/n
    recall=sum(v*v/true_size[t] for (t,p),v in intersections.items())/n
    same=lambda n:n*(n-1)//2
    tp=sum(same(v) for v in intersections.values());pp=sum(same(v) for v in pred_size.values());ap=sum(same(v) for v in true_size.values())
    return {'members':n,'truth_loci':len(true_size),'predicted_loci':len(pred_size),
            'b_cubed_precision':precision,'b_cubed_recall':recall,'b_cubed_f1':2*precision*recall/(precision+recall),
            'pairwise_precision':tp/pp if pp else None,'pairwise_recall':tp/ap if ap else None}


def run(binary,root,truth,condition,synteny):
    prefix=root/condition
    command=[str(binary),'graph','--manifest',str(root/'manifest.tsv'),'--alignments',str(root/'alignments.tsv'),
             '--flank','20','--candidate-window','10','--gzip-output','--output',str(prefix)]
    if synteny:command+=['--synteny',str(root/'synteny.tsv')]
    start=time.monotonic()
    with (root/(condition+'.log')).open('w') as handle:
        process=subprocess.Popen(command,stdout=handle,stderr=subprocess.STDOUT)
        _,status,usage=os.wait4(process.pid,0)
        process.returncode=os.waitstatus_to_exitcode(status) if hasattr(os,'waitstatus_to_exitcode') else (os.WEXITSTATUS(status) if os.WIFEXITED(status) else -os.WTERMSIG(status))
    report={'command':command,'binary_sha256':digest(binary),'condition':condition,'exit_code':process.returncode,
            'wall_seconds':time.monotonic()-start,'peak_rss_kib':usage.ru_maxrss if sys.platform!='darwin' else usage.ru_maxrss/1024,
            'user_seconds':usage.ru_utime,'system_seconds':usage.ru_stime}
    if process.returncode==0:
        report['metrics']=evaluate(prefix,truth)
        native=json.loads(Path(str(prefix)+'.run.json').read_text());report['counts']=native['counts']
        report['output_bytes']=sum(p.stat().st_size for p in root.glob(condition+'.*.gz'))
    (root/(condition+'.evaluation.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    if process.returncode:raise ValueError('simulation failed')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--tevox',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--mode',choices=['polyploid','scale'],required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    reports=[]
    cases=[(128,c,missing,seed) for c in [2,4,6] for missing in [False,True] for seed in [17,31,53]] if a.mode=='polyploid' else [(n,1,False,17) for n in [100,1000,10000,100000]]
    for n,c,missing,seed in cases:
        root=a.output/f'loci{n}_copies{c}_missing{int(missing)}_seed{seed}'
        truth=generate(root,n,c,seed,missing)
        for condition,synteny in ([('with_synteny',True),('without_synteny',False)] if a.mode=='polyploid' else [('sparse_pair',False)]):
            reports.append(run(a.tevox.resolve(),root.resolve(),truth,condition,synteny))
    (a.output/'summary.json').write_text(json.dumps({'scope':'CONTROLLED_SYNTHETIC_ONLY','runs':reports},indent=2)+'\n')
