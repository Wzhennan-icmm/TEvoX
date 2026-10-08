#!/usr/bin/env python3
"""Verify immutable run artifacts and recorded lossless archive bindings.

Hashes current output bytes, original selected inputs and archive sources.
Archived FASTA restoration uses the previously checked uncompressed digest
bound to the unchanged compressed source; optional --roundtrip rechecks it.
"""
import argparse
import datetime
import gzip
import hashlib
import json
from pathlib import Path


def checksum(path, decompress=False):
    result=hashlib.sha256()
    with (gzip.open(path,'rb') if decompress else open(path,'rb')) as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):result.update(block)
    return result.hexdigest()


def verify(root, roundtrip=False):
    root=Path(root).resolve()
    complete=json.loads((root/'complete.json').read_text())
    run=json.loads((root/'tevox.run.json').read_text())
    signature=json.loads((root/'signature.json').read_text())
    if complete['status']!='completed' or complete['signature']!=signature:
        raise ValueError('run completion/signature mismatch')
    archive_path=root/'archive.json';archives={}
    if archive_path.exists():
        archive=json.loads(archive_path.read_text())
        if archive['status']!='lossless_archive':raise ValueError('incomplete archive')
        archives={r['plain']:r for r in archive['files']}
        if len(archives)!=len(archive['files']):raise ValueError('duplicate archive entry')
    cache={};checked=[]
    def check(path,expected):
        path=Path(path);key=str(path)
        if key not in cache:cache[key]=checksum(path)
        if cache[key]!=expected:raise ValueError('SHA256 mismatch: '+key)
        checked.append({'path':key,'sha256':expected,'bytes':path.stat().st_size})
    def original(path,expected):
        p=Path(path)
        if p.exists():check(p,expected);return
        record=archives.get(str(p))
        if record is None or record['plain_sha256']!=expected:
            raise ValueError('missing input/output without matching archive: '+str(p))
        check(record['restore_from'],record['source_sha256'])
        if roundtrip:
            result=checksum(record['restore_from'],record['restore_method']=='gzip_decompress')
            if result!=expected:raise ValueError('archive round-trip mismatch: '+str(p))
    artifact_paths=set()
    for item in complete['artifacts']:
        original(item['path'],item['sha256']);artifact_paths.add(item['path'])
    if len(run['outputs'])!=17 or len(set(run['outputs']))!=17:
        raise ValueError('unexpected output contract')
    for name in run['outputs']:
        if str(root/('tevox.'+name)) not in artifact_paths:raise ValueError('unbound output '+name)
    for item in json.loads((root/'inputs.json').read_text()):check(item['path'],item['sha256'])
    for item in run['input_files']:original(item['path'],item['sha256'])
    for record in archives.values():
        original(record['plain'],record['plain_sha256'])
    exact=json.loads((root/'exact_identity.json').read_text())
    original(root/'alignment.eqx.paf',exact['output_sha256'])
    source_alignment=signature.get('alignment',{}).get('path',str(Path(signature['manifest']['path']).parent/'alignment.paf.gz'))
    check(source_alignment,exact['source_sha256'])
    check(signature['manifest']['path'],signature['manifest']['sha256'])
    unique={item['path']:item for item in checked}
    return {'status':'PASS','directory':str(root),'checked_file_count':len(unique),
            'checked_bytes':sum(r['bytes'] for r in unique.values()),
            'roundtrip_rechecked':roundtrip,'archived_original_count':len(archives),
            'binary_sha256':signature['binary']['sha256'],'counts':run['counts'],
            'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'files':list(unique.values())}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directories',nargs='+',type=Path)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--roundtrip',action='store_true');a=p.parse_args()
    results=[]
    for directory in a.directories:
        result=verify(directory,a.roundtrip);results.append(result)
        print(directory.name,result['status'],result['checked_file_count'],flush=True)
    a.output.write_text(json.dumps({'status':'PASS','runs':results},indent=2)+'\n')
