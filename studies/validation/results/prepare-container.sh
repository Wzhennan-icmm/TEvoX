#!/bin/bash
set -euo pipefail
. /etc/os-release
if command -v update-ca-trust >/dev/null; then
    cp /proxy-ca.crt /etc/pki/ca-trust/source/anchors/environment-proxy-ca.crt
    update-ca-trust extract
else
    mkdir -p /etc/ssl/certs
    printf '\n' >> /etc/ssl/certs/ca-certificates.crt
    cat /proxy-ca.crt >> /etc/ssl/certs/ca-certificates.crt
    chmod 755 /etc/ssl /etc/ssl/certs
    chmod 644 /etc/ssl/certs/ca-certificates.crt
fi
case "$ID" in
centos)
    sh /src/scripts/install-deps.sh --install --centos-vault
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
        sed 's/gcc make diffutils gawk python3/gcc make diffutils gawk python38/g' /src/scripts/install-deps.sh > /tmp/install-python.sh
        sh /tmp/install-python.sh --install --centos-vault
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
python3 -c 'import sys, csv, json, hashlib, dataclasses, math; assert sys.version_info >= (3,8)'
