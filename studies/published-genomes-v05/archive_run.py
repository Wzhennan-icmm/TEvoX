#!/usr/bin/env python3
"""Losslessly archive completed tables; record exact restoration of FASTAs."""
import argparse, gzip, hashlib, json, shutil
from pathlib import Path
from paf_exact import sha256

def checkpoint(root, records, status):
    temporary=root/'archive.json.part'
    temporary.write_text(json.dumps({'status':status,'files':records,
        'note':'run.json retains original paths/hashes; use restore_from and restore_method to recover original bytes'},indent=2)+'\n')
    temporary.replace(root/'archive.json')

def archive(directory):
    root=Path(directory).resolve()
    complete=root/'complete.json'
    if not complete.exists():raise ValueError('only completed runs may be archived')
    records=[]
    if (root/'archive.json').exists():
        prior=json.loads((root/'archive.json').read_text())
        if prior['status']=='lossless_archive':return
        records=prior['files']
    restored={r['plain']:r for r in records}
    manifest=json.loads((root/'source-manifest.json').read_text())
    run=json.loads((root/'tevox.run.json').read_text())
    for genome, native in zip(manifest['genomes'], run['genomes']):
        fasta=root/genome['id']/'genome.fa'
        if str(fasta) in restored:
            if fasta.exists():
                if sha256(fasta)!=restored[str(fasta)]['plain_sha256']:raise ValueError('resumed FASTA hash mismatch')
                fasta.unlink()
            continue
        if str(fasta) != native['fasta'] or sha256(fasta)!=native['fasta_sha256']:
            raise ValueError('FASTA path/hash mismatch')
        # Confirm the existing frozen source restores the exact analyzed bytes.
        opener=gzip.open if str(genome['fasta']).endswith('.gz') else open
        digest=hashlib.sha256()
        with opener(genome['fasta'],'rb') as handle:
            for chunk in iter(lambda:handle.read(4*1024*1024),b''):digest.update(chunk)
        if digest.hexdigest()!=native['fasta_sha256']:raise ValueError('source cannot restore analyzed FASTA')
        records.append({'plain':str(fasta),'plain_sha256':native['fasta_sha256'],
                        'restore_from':genome['fasta'],'restore_method':'gzip_decompress' if opener is gzip.open else 'copy',
                        'source_sha256':sha256(genome['fasta'])})
        checkpoint(root,records,'archiving')
        fasta.unlink()
    for path in sorted(root.rglob('*')):
        if not path.is_file() or path.suffix not in {'.tsv','.paf'} or path.stat().st_size<1024*1024:continue
        packed=Path(str(path)+'.gz');temporary=Path(str(packed)+'.part')
        if str(path) in restored:
            record=restored[str(path)]
            if sha256(packed)!=record['source_sha256'] or sha256(path)!=record['plain_sha256']:
                raise ValueError('resumed archive hash mismatch')
            path.unlink()
            continue
        if packed.exists():raise ValueError('archive destination exists')
        original_hash=sha256(path)
        with path.open('rb') as inp, temporary.open('wb') as raw:
            with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=1) as out:
                shutil.copyfileobj(inp,out,4*1024*1024)
        check=hashlib.sha256()
        with gzip.open(temporary,'rb') as handle:
            for chunk in iter(lambda:handle.read(4*1024*1024),b''):check.update(chunk)
        if check.hexdigest()!=original_hash:raise ValueError('gzip round-trip mismatch')
        temporary.replace(packed)
        records.append({'plain':str(path),'plain_sha256':original_hash,'plain_bytes':path.stat().st_size,
                        'restore_from':str(packed),'restore_method':'gzip_decompress','source_sha256':sha256(packed),
                        'compressed_bytes':packed.stat().st_size})
        checkpoint(root,records,'archiving')
        path.unlink()
    checkpoint(root,records,'lossless_archive')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory');archive(p.parse_args().directory)
