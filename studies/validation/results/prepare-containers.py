import concurrent.futures, datetime, json, subprocess
from pathlib import Path

root=Path(__file__).resolve().parent
ca=Path('/workspace/tevox-build/rpm-validation/centos8/input/ca/environment-proxy-ca.crt')
repo=Path('/workspace/tevox-review-20261003/early-development/v0.5-source')
script=root/'prepare-container.sh'
script.write_text(r'''#!/bin/bash
set -euo pipefail
. /etc/os-release
if command -v update-ca-trust >/dev/null; then
    cp /proxy-ca.crt /etc/pki/ca-trust/source/anchors/environment-proxy-ca.crt
    update-ca-trust extract
else
    mkdir -p /etc/ssl/certs
    printf '\n' >> /etc/ssl/certs/ca-certificates.crt
    cat /proxy-ca.crt >> /etc/ssl/certs/ca-certificates.crt
    chmod 644 /etc/ssl/certs/ca-certificates.crt
    chmod 755 /etc/ssl /etc/ssl/certs
fi
case "$ID" in
centos)
    sed 's/gawk gzip/gawk gzip zlib-devel/g' /src/scripts/install-deps.sh > /tmp/deps.sh
    sh /tmp/deps.sh --install --centos-vault
    if [[ "$VERSION_ID" == 7* ]]; then
        mkdir -p /tmp/python-source
        cp -a /python-source/. /tmp/python-source/
        cd /tmp/python-source
        ./configure --prefix=/opt/tevox-python --without-ensurepip
        make -j2
        make install
        ln -sf python3.8 /opt/tevox-python/bin/python3
        rm -rf /tmp/python-source
    else
        sh /src/scripts/install-deps.sh --install --centos-vault
        mkdir -p /opt/tevox-python/bin
        ln -s /usr/bin/python3.8 /opt/tevox-python/bin/python3
    fi
    ;;
rocky)
    major=${VERSION_ID%%.*}
    case "$major" in 8) key=/etc/pki/rpm-gpg/RPM-GPG-KEY-rockyofficial;; *) key=/etc/pki/rpm-gpg/RPM-GPG-KEY-Rocky-$major;; esac
    mkdir -p /tmp/official-repos
    for component in BaseOS AppStream; do
        cat > /tmp/official-repos/$component.repo <<EOF
[tevox-$component]
name=Rocky $major $component
baseurl=https://dl.rockylinux.org/pub/rocky/$major/$component/\$basearch/os/
enabled=1
gpgcheck=1
sslverify=1
gpgkey=file://$key
EOF
    done
    packages='gcc make diffutils gawk python3'
    if [[ "$major" == 8 ]]; then packages="$packages python38"; fi
    dnf --setopt=reposdir=/tmp/official-repos --setopt=install_weak_deps=False -y install $packages
    dnf clean all
    if [[ "$major" == 8 ]]; then
        mkdir -p /opt/tevox-python/bin
        ln -s /usr/bin/python3.8 /opt/tevox-python/bin/python3
    fi
    ;;
ubuntu)
    printf 'Acquire::https::CaInfo "/etc/ssl/certs/ca-certificates.crt";\n' > /etc/apt/apt.conf.d/80tevox-ca
    sed -i 's@http://archive.ubuntu.com@https://archive.ubuntu.com@g; s@http://security.ubuntu.com@https://security.ubuntu.com@g' /etc/apt/sources.list /etc/apt/sources.list.d/*.sources 2>/dev/null || true
    sh /src/scripts/install-deps.sh --install
    apt-get clean
    ;;
*) exit 2;;
esac
export PATH=/opt/tevox-python/bin:$PATH
gcc --version
python3 --version
python3 -c 'import sys, csv, json, hashlib, dataclasses, math, gzip, zlib; assert sys.version_info >= (3,8)'
''')
matrix=[('centos7','tevox-final-index-centos7:local'),('centos8','tevox-final-index-centos8:local'),('rocky8','rockylinux:8'),('rocky9','rockylinux:9'),('rocky10','tevox-final-index-rocky10:local'),('ubuntu20','ubuntu:20.04'),('ubuntu22','ubuntu:22.04'),('ubuntu24','ubuntu:24.04'),('ubuntu26','ubuntu:26.04')]
def prepare(item):
    label,image=item;d=root/'platforms'/label;d.mkdir(parents=True,exist_ok=True)
    name='tevox-v05-tools-'+label;target=name+':local'
    cmd=['docker','create','--name',name,'--network','host','--cpus','2','--memory','2g',
         '-v',str(ca)+':/proxy-ca.crt:ro','-v',str(repo)+':/src:ro',
         '-v',str(root/'cpython-3.8.20')+':/python-source:ro',
         '-v',str(script)+':/prepare.sh:ro']
    for key in ['HTTPS_PROXY','HTTP_PROXY','NO_PROXY','https_proxy','http_proxy','no_proxy']:cmd+=['-e',key]
    cmd += [image,'bash','/prepare.sh']
    record={'image':image,'command':cmd,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        subprocess.run(cmd,check=True,capture_output=True)
        with (d/'prepare.log').open('w') as log:r=subprocess.run(['docker','start','--attach',name],stdout=log,stderr=subprocess.STDOUT)
        record['exit_code']=r.returncode
        if r.returncode==0:
            record['prepared_image_id']=subprocess.check_output(['docker','commit',name,target],text=True).strip()
            record['prepared_image']=target
            subprocess.run(['docker','rm',name],check=True,capture_output=True)
    except Exception as e:record['error']=str(e)
    (d/'prepare.json').write_text(json.dumps(record,indent=2)+'\n')
    print(label,record.get('exit_code'),record.get('error',''),flush=True)
    return record
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    records=list(executor.map(prepare,matrix))
(root/'platforms/preparation-summary.json').write_text(json.dumps(records,indent=2)+'\n')
