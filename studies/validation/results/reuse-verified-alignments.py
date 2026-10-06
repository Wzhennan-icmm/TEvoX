import pathlib,json,gzip,hashlib,shutil,datetime
repo=pathlib.Path('/workspace/tevox-review-20261003/early-development/v0.5-source')
oldroot=pathlib.Path('/workspace/tevox-v05-study/results');newroot=pathlib.Path('/workspace/tevox-v05-study/final')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
for name in ['rat_corrected','maize']:
 old=oldroot/name;out=newroot/name
 signature=json.loads((old/'signature.json').read_text())
 assert signature['scripts']['paf_exact.py']==sha(repo/'studies/published-genomes-v05/paf_exact.py')
 assert signature['manifest']['sha256']==sha(signature['manifest']['path'])
 inputs=json.loads((old/'inputs.json').read_text())
 for r in inputs:assert r['sha256']==sha(r['path'])
 exact=json.loads((old/'exact_identity.json').read_text());source=pathlib.Path(signature['manifest']['path']).parent/'alignment.paf.gz'
 assert exact['source_sha256']==sha(source)
 archive=json.loads((old/'archive.json').read_text());assert archive['status']=='lossless_archive'
 record=next(r for r in archive['files'] if r['plain']==str(old/'alignment.eqx.paf'))
 assert record['plain_sha256']==exact['output_sha256'] and record['source_sha256']==sha(record['restore_from'])
 run=json.loads((old/'tevox.run.json').read_text())
 fasta_bindings=[]
 for genome in run['genomes']:
  r=next(r for r in archive['files'] if r['plain']==genome['fasta'])
  assert r['plain_sha256']==genome['fasta_sha256'] and r['source_sha256']==sha(r['restore_from'])
  fasta_bindings.append(r)
 out.mkdir(exist_ok=False)
 dest=out/'alignment.eqx.paf'
 with gzip.open(record['restore_from'],'rb') as f,dest.with_suffix('.paf.part').open('wb') as w:shutil.copyfileobj(f,w,4*1024*1024)
 assert sha(dest.with_suffix('.paf.part'))==exact['output_sha256']
 dest.with_suffix('.paf.part').replace(dest)
 (out/'exact_identity.json').write_text(json.dumps(exact,indent=2)+'\n')
 proof={'status':'VERIFIED_EXACT_ALIGNMENT_CACHE_REUSE','source_run':str(old),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'converter_sha256':signature['scripts']['paf_exact.py'],'exact_alignment':exact,'archive_record':record,'fasta_archive_bindings':fasta_bindings,'raw_inputs':inputs,'note':'Only unchanged source-derived EQX alignment reused; final TEvoX inference and all postprocessing rerun.'}
 (out/'alignment-reuse.json').write_text(json.dumps(proof,indent=2)+'\n')
 print(name,'verified cache',dest.stat().st_size,flush=True)
