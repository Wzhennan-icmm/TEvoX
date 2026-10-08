import concurrent.futures, datetime, hashlib, json, os, shutil, subprocess
from pathlib import Path
root=Path('/workspace/tevox-v05-validation-20261005')
evidence=Path('/workspace/tevox-v05-validation-20261006')
repo=Path('/workspace/tevox-review-20261003/early-development/v0.5-source')
source=evidence/'matrix-source'
shutil.copytree(repo,source,ignore=shutil.ignore_patterns('.git','build','tevox','__pycache__'))
hashes={str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.rglob('*')) if p.is_file()}
(evidence/'matrix-source-sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
script=evidence/'run-matrix.sh'
script.write_text('''#!/bin/bash
set -euo pipefail
export PATH=/opt/tevox-python/bin:$PATH
export PYTHONDONTWRITEBYTECODE=1
cd /src
cat /etc/os-release
gcc --version
python3 --version
make -j1 CC=gcc OBJDIR=/evidence/build TARGET=/evidence/tevox CFLAGS='-O2 -g -Werror' check
make OBJDIR=/evidence/build TARGET=/evidence/tevox PREFIX=/usr DESTDIR=/evidence/stage install
cmp /evidence/tevox /evidence/stage/usr/bin/tevox
/evidence/stage/usr/bin/tevox --version
for command in /evidence/stage/usr/bin/tevox-*; do "$command" --help > /dev/null; done
TEVOX_BIN=/evidence/stage/usr/bin/tevox bash tests/run_tests.sh
echo TEVOX_V05_MATRIX_PASS
''')
def run(label):
    d=evidence/'platforms'/label
    d.mkdir(parents=True,exist_ok=True)
    prepared=json.loads((root/'platforms'/label/'prepare.json').read_text())
    if prepared.get('exit_code')!=0:return {'distribution':label,'status':'PREPARATION_FAILED'}
    image=prepared['prepared_image']
    command=['docker','run','--rm','--network','none','--cpus','1','--memory','2g','--read-only',
             '--tmpfs','/tmp:rw,nosuid,nodev,exec,size=536870912',
             '--user',str(os.getuid())+':'+str(os.getgid()),
             '-v',str(source)+':/src:ro','-v',str(d)+':/evidence:rw',
             '-v',str(script)+':/run.sh:ro',image,'bash','/run.sh']
    started=datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (d/'validation.log').open('w') as handle:r=subprocess.run(command,stdout=handle,stderr=subprocess.STDOUT)
    passed=r.returncode==0 and 'TEVOX_V05_MATRIX_PASS' in (d/'validation.log').read_text()
    result=dict(distribution=label,status='PASS' if passed else 'FAIL',exit_code=r.returncode,command=command,started_utc=started,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_hashes=str(evidence/'matrix-source-sha256.json'),image_id=prepared['prepared_image_id'])
    (d/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(label,result['status'],flush=True)
    return result
labels=['centos7','centos8','rocky8','rocky9','rocky10','ubuntu20','ubuntu22','ubuntu24','ubuntu26']
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    results=list(executor.map(run,labels))
(evidence/'matrix-summary.json').write_text(json.dumps(results,indent=2)+'\n')
