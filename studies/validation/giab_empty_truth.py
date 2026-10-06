#!/usr/bin/env python3
"""Conservative external HG002 homozygous deletion subset for CHM13 TE absence.

Never treats absent VCF records as negatives. Requires PASS, 1/1, Tier1,
unique ungapped and reciprocal chain mapping including flanks, exact REF
sequence agreement and reciprocal TE/deletion overlap. Output is a positive
absence benchmark; it is not candidate-membership calibration truth.
"""
import argparse
import bisect
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'published-genomes-v05'))
from paf_exact import Fasta, COMPLEMENT, sha256


class Intervals:
    def __init__(self, rows):
        self.rows=sorted(rows, key=lambda r:(r[0],r[1]))
        self.starts=[r[0] for r in self.rows]
        self.maximum=[]
        end=-1
        for r in self.rows:
            end=max(end,r[1]); self.maximum.append(end)

    def overlap(self, start, end):
        i=bisect.bisect_left(self.starts,end)-1
        while i>=0 and self.maximum[i]>start:
            row=self.rows[i]
            if row[1]>start: yield row
            i-=1


class Chains:
    def __init__(self, path):
        blocks=collections.defaultdict(list)
        with gzip.open(path,'rt') as handle:
            for line in handle:
                f=line.split()
                if not f:continue
                if f[0]=='chain':
                    _,score,chrom,size,strand,start,end,target,tsize,tstrand,tstart,tend,key=f
                    if strand!='+':raise ValueError('unsupported chain reference orientation')
                    a,b=int(start),int(tstart); target_size=int(tsize)
                else:
                    n=int(f[0])
                    bs,be=(b,b+n) if tstrand=='+' else (target_size-b-n,target_size-b)
                    blocks[chrom].append((a,a+n,target,bs,be,tstrand,key))
                    if len(f)==3: a+=n+int(f[1]); b+=n+int(f[2])
        self.blocks={k:Intervals(v) for k,v in blocks.items()}

    def map(self, chrom, start, end):
        result=set()
        for a,b,target,bs,be,strand,key in self.blocks.get(chrom,Intervals([])).overlap(start,end):
            if a<=start and end<=b:
                x,y=(bs+start-a,bs+end-a) if strand=='+' else (be-(end-a),be-(start-a))
                result.add((target,x,y,strand))
        return list(result)


def prepare(args):
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    counters=collections.Counter()
    regions=collections.defaultdict(list)
    with open(args.regions) as handle:
        for line in handle:
            if line.startswith(('#','track')):continue
            chrom,a,b,*_=line.split(); regions[chrom.removeprefix('chr') if hasattr(chrom,'removeprefix') else chrom.replace('chr','',1)].append((int(a),int(b)))
    regions={k:Intervals(v) for k,v in regions.items()}
    forward,reverse=Chains(args.forward),Chains(args.reverse)
    fasta=Fasta(args.fasta); accepted=collections.defaultdict(list)
    try:
        with gzip.open(args.vcf,'rt') as handle:
            for line in handle:
                if line.startswith('#CHROM'):
                    if line.rstrip().split('\t')[9:]!=['HG002']:raise ValueError('expected one HG002 sample')
                if line.startswith('#'):continue
                f=line.rstrip().split('\t');chrom,pos,key,ref,alt,qual,filt,info,fmt,call=f
                counters['vcf_records']+=1
                tags=dict(s.split('=',1) for s in info.split(';') if '=' in s)
                gt=dict(zip(fmt.split(':'),call.split(':'))).get('GT')
                if filt!='PASS' or gt not in {'1/1','1|1'} or tags.get('SVTYPE')!='DEL':continue
                counters['pass_homozygous_deletion']+=1
                if len(alt)!=1 or not ref.startswith(alt) or len(ref)-1<50 or set(ref.upper())-set('ACGT'):continue
                start=int(pos)-1;end=start+len(ref);left=start-args.flank;right=end+args.flank
                if left<0 or not any(a<=left and right<=b for a,b in regions.get(chrom,Intervals([])).overlap(left,right)):
                    counters['outside_tier1_with_flanks']+=1;continue
                maps=forward.map('chr'+chrom,left,right)
                if len(maps)!=1:counters['nonunique_or_gapped_mapping']+=1;continue
                t,x,y,strand=maps[0]
                back=reverse.map(t,x,y)
                if back!=[('chr'+chrom,left,right,strand)]:counters['nonreciprocal_mapping']+=1;continue
                rs,re=x+args.flank,y-args.flank
                sequence=fasta.fetch(t,rs,re)
                if strand=='-':sequence=sequence.translate(COMPLEMENT)[::-1]
                if sequence!=ref.upper().encode():counters['chm13_ref_sequence_mismatch']+=1;continue
                ds,de=(rs+1,re) if strand=='+' else (rs,re-1)
                accepted[t].append((ds,de,key,chrom,start,end,strand,hashlib.sha256(sequence).hexdigest()))
                counters['compatible_deletions']+=1
    finally:fasta.close()
    index={k:Intervals(v) for k,v in accepted.items()}
    fields=['truth_record_id','source_genome_id','source_te_id','target_genome_id','chrom','te_start','te_end','sv_start','sv_end','giab_variant_id','grch37_chrom','grch37_ref_start','grch37_ref_end','mapping_strand','te_overlap_fraction','deletion_overlap_fraction','expected_state','validation_source','ref_sha256']
    with (out/'truth.tsv').open('w',newline='') as dest, open(args.te) as source:
        writer=csv.DictWriter(dest,fieldnames=fields,delimiter='\t');writer.writeheader()
        for line in source:
            if line.startswith('#') or not line.strip():continue
            c,a,b,key,score,strand,family,kind=line.rstrip().split('\t');a,b=int(a),int(b)
            matches=[]
            for row in index.get(c,Intervals([])).overlap(a,b):
                ds,de,sv,oldc,oldstart,oldend,orientation,digest=row
                overlap=max(0,min(b,de)-max(a,ds))
                if overlap/(b-a)>=args.te_coverage and overlap/(de-ds)>=args.deletion_coverage:matches.append(row)
            if len(matches)!=1:continue
            ds,de,sv,oldc,oldstart,oldend,orientation,digest=matches[0]
            overlap=min(b,de)-max(a,ds);counters['positive_te_fragments']+=1
            writer.writerow(dict(zip(fields,[f'GIAB_TE_{counters["positive_te_fragments"]:06}',args.source_genome,key,args.target_genome,c,a,b,ds,de,sv,oldc,oldstart,oldend,orientation,overlap/(b-a),overlap/(de-ds),'EMPTY_SITE_CONFIRMED','GIAB_HG002_SV_Tier1_v0.6',digest])))
    report={'status':'EXTERNAL_POSITIVE_ABSENCE_SUBSET','counts':dict(counters),
            'parameters':vars(args),'inputs':{name:{'path':str(getattr(args,name)),'sha256':sha256(getattr(args,name))} for name in ['vcf','regions','forward','reverse','fasta','te']},
            'limitations':['No negative labels: precision/FDR and membership calibration cannot be estimated.',
             'Read data may overlap assembly inputs; independence of all experimental reads is not established.',
             'Strict transfer/overlap filters select simple reference-like loci and do not represent genome-wide performance.',
             'Homozygous HG002 genotypes avoid inferring haplotype from unphased heterozygous calls.']}
    (out/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['counts']))


def evaluate(args):
    with open(args.truth) as handle:truth=list(csv.DictReader(handle,delimiter='\t'))
    keys={(r['source_genome_id'],r['source_te_id'],r['target_genome_id']):r for r in truth}
    if len(keys)!=len(truth):raise ValueError('duplicate truth TE')
    seen={};counts=collections.Counter()
    opener=gzip.open if str(args.decisions).endswith('.gz') else open
    with opener(args.decisions,'rt') as handle:
        for r in csv.DictReader(handle,delimiter='\t'):
            key=(r['source_genome_id'],r['source_te_id'],r['target_genome_id'])
            if key not in keys:continue
            if key in seen:raise ValueError('duplicate decision for truth TE')
            seen[key]=r;counts[r['legacy_state']]+=1
    success=sum(r['legacy_state']=='EMPTY_SITE_CONFIRMED' and r['claimable']=='true' for r in seen.values())
    report={'truth_fragments':len(truth),'matched_decisions':len(seen),'missing_decisions':len(truth)-len(seen),
            'decision_states':dict(counts),'claimable_empty':success,
            'positive_subset_recall':success/len(truth) if truth else None,
            'precision':None,'false_discovery_rate':None,'calibration_status':'NOT_ASSESSABLE_POSITIVE_ONLY',
            'truth_sha256':sha256(args.truth),'decisions_sha256':sha256(args.decisions)}
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for name in ['vcf','regions','forward','reverse','fasta','te','output']:a.add_argument('--'+name,required=True)
    a.add_argument('--source-genome',default='chm13');a.add_argument('--target-genome',default='hg002_mat')
    a.add_argument('--flank',type=int,default=50);a.add_argument('--te-coverage',type=float,default=.95);a.add_argument('--deletion-coverage',type=float,default=.8)
    b=sub.add_parser('evaluate')
    for name in ['truth','decisions','output']:b.add_argument('--'+name,required=True)
    args=p.parse_args();(prepare if args.command=='prepare' else evaluate)(args)
