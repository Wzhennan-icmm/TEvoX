#!/usr/bin/env python3
"""Run the v0.5 engine on a published pair using frozen legacy source manifests.

Reuses only the original FASTA, repeat/gene annotations and alignment path.
Does not reinterpret old unique-TE results as v0.5 states. Exact identity is
recomputed from both FASTAs. Each output directory belongs to one immutable run.
"""
import argparse
import collections
import csv
import datetime
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import annotate_candidates
import normalize_annotations
from paf_exact import Fasta, convert, sha256


def save(path, value):
    temporary = Path(str(path)+'.part')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n')
    temporary.replace(path)


def fingerprint(path):
    path = Path(path).resolve()
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha256(path)}


def rows(path):
    path=Path(path)
    if not path.exists(): path=Path(str(path)+'.gz')
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def run(args):
    source = Path(args.legacy_result).resolve() if args.legacy_result else None
    manifest_path = Path(args.manifest).resolve() if args.manifest else source/'manifest.json'
    alignment_source = Path(args.paf).resolve() if args.paf else source/'alignment.paf.gz'
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out/'complete.json').exists():
        raise ValueError('completed run exists; use another output directory')
    manifest = json.loads(manifest_path.read_text())
    if len(manifest['genomes']) != 2:
        raise ValueError('exactly two genomes are required, in PAF query/target order')
    if len({g['id'] for g in manifest['genomes']}) != 2:
        raise ValueError('genome identifiers must be unique')
    inputs = []
    for genome in manifest['genomes']:
        if not genome['id'] or genome['id'] in {'.','..'} or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-' for c in genome['id']):
            raise ValueError('genome identifier must be a safe directory name')
        for kind in ['fasta','repeats','genes']:
            path=Path(genome[kind])
            genome[kind]=str((path if path.is_absolute() else manifest_path.parent/path).resolve())
            inputs.append(dict(role=kind,genome=genome['id'],**fingerprint(genome[kind])))
    binary = Path(args.tevox).resolve()
    signature = {'manifest': fingerprint(manifest_path), 'binary': fingerprint(binary),
                 'alignment':fingerprint(alignment_source),'source_inputs':inputs,
                 'family_policy': args.family_policy,
                 'scripts': {p.name: sha256(p) for p in sorted(Path(__file__).parent.glob('*.py'))}}
    if (out/'signature.json').exists() and json.loads((out/'signature.json').read_text()) != signature:
        raise ValueError('source/configuration changed; use a new output directory')
    save(out/'signature.json', signature)
    save(out/'source-manifest.json', manifest)
    prepared = []
    fasta_objects = []
    def event(phase, **values):
        record = dict(phase=phase, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), **values)
        with (out/'events.jsonl').open('a') as handle: handle.write(json.dumps(record)+'\n')
        print(json.dumps(record), flush=True)

    try:
        for genome in manifest['genomes']:
            d = out/genome['id']; d.mkdir(exist_ok=True)
            fasta = d/'genome.fa'
            fasta_record = d/'fasta-cache.json'
            event('prepare_fasta', genome=genome['id'])
            if not fasta.exists():
                opener = gzip.open if str(genome['fasta']).endswith('.gz') else open
                with opener(genome['fasta'], 'rb') as inp, (d/'genome.fa.part').open('wb') as dest:
                    shutil.copyfileobj(inp, dest, 4*1024*1024)
                (d/'genome.fa.part').replace(fasta)
                save(fasta_record,fingerprint(fasta))
            elif not fasta_record.exists() or json.loads(fasta_record.read_text()) != fingerprint(fasta):
                raise ValueError('unverified or changed prepared FASTA: '+str(fasta))
            seq = Fasta(fasta); fasta_objects.append(seq)
            lengths = d/'lengths.tsv'
            lengths.write_text(''.join(name+'\t'+str(row[0])+'\n' for name,row in seq.sequences.items()))
            prefix = d/'repeats'
            event('prepare_annotations', genome=genome['id'])
            normalizer_args = normalize_annotations.parser().parse_args([
                '--input',genome['repeats'],'--format',genome['repeat_format'],
                '--lengths',str(lengths),'--output-prefix',str(prefix),
                '--id-prefix',genome['id'],'--family-policy',args.family_policy,
                *genome.get('normalize_args',[])])
            qc = normalize_annotations.normalize(normalizer_args)
            allowed = {'rejected_excluded_non_te', 'rejected_excluded_feature_type', 'rejected_unresolved_class'}
            invalid = {k:v for k,v in qc['counts'].items() if k.startswith('rejected_') and k not in allowed and v}
            if invalid or not qc['counts'].get('retained_fragments'):
                raise ValueError('invalid annotations: '+str(invalid))
            prepared.append(dict(id=genome['id'], fasta=str(fasta), te=str(prefix)+'.bed', genes=genome['genes']))
        paf = out/'alignment.eqx.paf'
        event('exact_identity')
        if not (out/'exact_identity.json').exists():
            if paf.exists(): raise ValueError('alignment without provenance marker')
            exact = convert(alignment_source, paf, *fasta_objects)
            save(out/'exact_identity.json', exact)
        elif sha256(paf) != json.loads((out/'exact_identity.json').read_text())['output_sha256']:
            raise ValueError('exact alignment hash mismatch')
    finally:
        for seq in fasta_objects: seq.close()
    prefix = out/'tevox'
    command = [str(binary),'pair']
    for suffix, genome in zip(['a','b'],prepared):
        command += ['--genome-'+suffix,genome['id'],'--fasta-'+suffix,genome['fasta'],'--te-'+suffix,genome['te']]
    command += ['--paf',str(paf),'--output',str(prefix),'--gzip-output']
    save(out/'inputs.json', inputs)
    if sha256(binary) != signature['binary']['sha256']:
        raise ValueError('TEvoX binary changed during preparation; use the recorded binary')
    event('tevox_start', argv=command)
    started = time.monotonic()
    with (out/'tevox.stdout.log').open('w') as stdout, (out/'tevox.stderr.log').open('w') as stderr:
        result = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        _, wait_status, usage = os.wait4(result.pid, 0)
        result.returncode = os.WEXITSTATUS(wait_status) if os.WIFEXITED(wait_status) else -os.WTERMSIG(wait_status)
        save(out/'resource_usage.json', {'wall_seconds':time.monotonic()-started,
             'user_seconds':usage.ru_utime, 'system_seconds':usage.ru_stime,
             'peak_rss_kib':usage.ru_maxrss if sys.platform != 'darwin' else usage.ru_maxrss/1024,
             'method':'wait4 child resource usage'})
    event('tevox_finished',exit_code=result.returncode,seconds=time.monotonic()-started)
    if result.returncode: raise ValueError('TEvoX failed; inspect stderr and time.txt')
    run_json = json.loads(Path(str(prefix)+'.run.json').read_text())
    if run_json['schema_version'] != '1.2.0': raise ValueError('unexpected schema')
    state_counts = collections.Counter()
    for row in rows(str(prefix)+'.states.tsv'):
        state_counts[(row['genome_id'],row['state'],row['claimable'])] += 1
        if row['claimable']=='true' and row['technical_state'] != 'CALLABLE':
            raise ValueError('claimable non-callable state')
    counts = collections.Counter()
    expected = {'evidence.tsv':'evidence_observations','observation_scores.tsv':'observation_score_rows',
                'candidates.tsv':'candidates','candidate_features.tsv':'candidate_feature_rows',
                'candidate_contexts.tsv':'candidates','decisions.tsv':'decisions','edges.tsv':'edges',
                'relations.tsv':'relations','solver.tsv':'solver_components','loci.tsv':'loci',
                'synteny.blocks.tsv':'synteny_blocks','synteny.anchors.tsv':'synteny_anchors',
                'contexts.tsv':'copy_contexts','te_contexts.tsv':'te_context_assignments'}
    for suffix in run_json['outputs']:
        path=Path(str(prefix)+'.'+suffix)
        opener=gzip.open if path.suffix=='.gz' else open
        with opener(path,'rb') as handle:
            first=handle.readline()
            if not first: raise ValueError('missing header '+suffix)
            counts[suffix]=sum(block.count(b'\n') for block in iter(lambda:handle.read(4*1024*1024),b''))
        contract_suffix=suffix[:-3] if suffix.endswith('.gz') else suffix
        count_expected = run_json['counts'][expected[contract_suffix]] if contract_suffix in expected else None
        if contract_suffix in {'instances.tsv','states.tsv'}: count_expected=run_json['counts']['loci']*run_json['counts']['genomes']
        if count_expected is not None and counts[suffix] != count_expected:
            raise ValueError('row count mismatch '+suffix)
    if run_json['performance']['fasta_reopens_after_index'] != 0: raise ValueError('FASTA reopen invariant')
    # Gene context for annotation fragments with a claimable empty-site decision.
    # These remain decision-level fragments, not independent insertion events.
    empty = collections.defaultdict(set)
    for row in rows(str(prefix)+'.decisions.tsv'):
        if row['legacy_state']=='EMPTY_SITE_CONFIRMED' and row['claimable']=='true':
            empty[row['source_genome_id']].add(row['source_te_id'])
    for genome in prepared:
        d=out/genome['id']; candidate=d/'empty_site_source_fragments.tsv'
        with candidate.open('w') as handle, open(genome['te']) as bed:
            handle.write('# ID\tChr\tStart\tEnd\tStrand\tType\tFamily\tName\n')
            for line in bed:
                chrom,start,end,key,score,strand,family,kind=line.rstrip('\n').split('\t')
                if key in empty[genome['id']]:
                    handle.write('\t'.join([key,chrom,str(int(start)+1),end,strand,kind,family,key])+'\n')
        gene_qc=annotate_candidates.annotate(SimpleNamespace(candidates=str(candidate), genes=genome['genes'],
            metadata=str(d/'repeats.metadata.tsv'),lengths=str(d/'lengths.tsv'),output=str(d/'empty_site_genes.tsv')))
        if any(k.startswith('rejected_') and v for k,v in gene_qc['gene_counts'].items()): raise ValueError('invalid gene annotation')
    summary={'source_pair':str(source or manifest_path),'counts':run_json['counts'],
             'state_counts':[dict(genome_id=k[0],state=k[1],claimable=k[2],count=v) for k,v in sorted(state_counts.items())],
             'empty_site_source_fragments':{k:len(v) for k,v in empty.items()},
             'row_counts':dict(counts),'validation':'streaming counts/state gates; full schema validator separately required on bounded fixtures',
             'family_policy':args.family_policy,'calibration_status':'UNCALIBRATED',
             'interpretation':'TE annotation fragments and reconstructed loci; no insertion/deletion evolutionary direction or independent accuracy claim'}
    save(out/'summary.json',summary)
    artifacts=[]
    for path in sorted(out.glob('tevox.*')):
        artifacts.append(fingerprint(path))
    save(out/'complete.json',dict(status='completed',signature=signature,artifacts=artifacts,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    event('complete',loci=run_json['counts']['loci'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    source_args=p.add_mutually_exclusive_group(required=True)
    source_args.add_argument('--legacy-result',help='Frozen directory containing manifest.json and alignment.paf.gz')
    source_args.add_argument('--manifest',help='Pair JSON; relative genome paths resolve beside this file')
    p.add_argument('--paf',help='PAF in manifest query/target order; required with --manifest')
    p.add_argument('--tevox',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--family-policy',choices=['repeat-name','classification'],default='repeat-name')
    args=p.parse_args()
    if args.manifest and not args.paf:p.error('--manifest requires --paf')
    run(args)
