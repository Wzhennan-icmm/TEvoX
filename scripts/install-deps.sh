#!/bin/sh
# Print platform-specific prerequisites; install only with explicit --install.
set -eu

usage() {
    cat <<'EOF'
Usage: sh scripts/install-deps.sh [--install] [--centos-vault]

With no arguments, print the dependency installation commands for this system.
--install runs those commands. On Linux, run it as root (or explicitly via sudo).
On macOS, --install requests Apple's Command Line Tools if they are missing.
--centos-vault explicitly selects the official archived CentOS Linux 7/8 packages.
It uses a temporary repository directory, without changing system repo files.
Package signature and TLS verification remain enabled. --help shows this text.
EOF
}

mode=print
centos_vault=no
for option in "$@"; do
    case "$option" in
        --install) mode=install ;;
        --centos-vault) centos_vault=yes ;;
        --help|-h) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done

case "$(uname -s)" in
    Darwin)
        if [ "$centos_vault" = yes ]; then
            printf '%s\n' '--centos-vault is only supported on CentOS Linux 7/8.' >&2
            exit 2
        fi
        printf '%s\n' 'macOS (Intel or Apple Silicon): Apple Command Line Tools provide Clang and GNU Make.'
        if xcode-select -p >/dev/null 2>&1; then
            printf '%s\n' 'Apple developer tools are already selected. Python 3.8+ is also required for helpers and tests.'
        elif [ "$mode" = install ]; then
            xcode-select --install
            printf '%s\n' 'Finish the Apple installer, then run make and make check.'
        else
            printf '%s\n' 'Install command: xcode-select --install'
        fi
        exit 0
        ;;
    Linux) ;;
    *) printf '%s\n' 'Unsupported OS. Install a C11 compiler and GNU Make 3.81 or newer.' >&2; exit 1 ;;
esac

if [ ! -r /etc/os-release ]; then
    printf '%s\n' 'Cannot identify this Linux distribution: /etc/os-release is missing.' >&2
    exit 1
fi
. /etc/os-release
printf 'Platform: %s (%s)\n' "${PRETTY_NAME:-${ID:-Linux}}" "$(uname -m)"

if [ "$centos_vault" = yes ]; then
    if [ "${ID:-}" != centos ]; then
        printf '%s\n' '--centos-vault is only supported on CentOS Linux 7/8.' >&2
        exit 2
    fi
    case "${NAME:-} ${PRETTY_NAME:-}" in
        *Stream*) printf '%s\n' '--centos-vault does not support CentOS Stream.' >&2; exit 2 ;;
    esac
    centos_version=${VERSION_ID:-0}
    case "${centos_version%%.*}" in
        7)
            vault_version=7.9.2009
            vault_components='os updates extras'
            vault_suffix=
            vault_key=/etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7
            ;;
        8)
            vault_version=8.5.2111
            vault_components='BaseOS AppStream extras'
            vault_suffix=os/
            vault_key=/etc/pki/rpm-gpg/RPM-GPG-KEY-centosofficial
            ;;
        *) printf '%s\n' '--centos-vault requires CentOS Linux version 7 or 8.' >&2; exit 2 ;;
    esac
fi

python_package=python3
case "${ID:-}:${VERSION_ID:-}" in
    centos:8*|rocky:8*|rhel:8*|almalinux:8*) python_package=python38 ;;
esac
case "${ID:-}" in
    ubuntu|debian)
        manager=apt-get
        printf '%s\n' 'Install commands: apt-get update && apt-get install -y --no-install-recommends build-essential diffutils gawk gzip python3'
        ;;
    centos)
        centos_version=${VERSION_ID:-0}
        case "${centos_version%%.*}" in
            7) manager=yum ;;
            *) manager=dnf ;;
        esac
        if [ "$centos_vault" = yes ]; then
            printf 'Archive: https://vault.centos.org/%s/ (%s)\n' "$vault_version" "$vault_components"
            printf 'Package signing key: %s\n' "$vault_key"
            printf 'Installer: %s, with a temporary repository directory; packages: gcc make diffutils gawk gzip %s\n' "$manager" "$python_package"
            printf '%s\n' 'Install command: sh scripts/install-deps.sh --install --centos-vault'
        else
            printf 'Install command: %s install -y gcc make diffutils gawk gzip %s\n' "$manager" "$python_package"
            printf '%s\n' 'CentOS 7 and 8 are end of life. Their default mirrors may be unavailable.' \
                'To opt into the official archive, run this script with --centos-vault.'
        fi
        ;;
    rocky|rhel|almalinux)
        manager=dnf
        printf 'Install command: dnf install -y gcc make diffutils gawk gzip %s\n' "$python_package"
        ;;
    *)
        printf '%s\n' 'No automatic installer for this distribution. Install a C11 compiler and GNU Make 3.81+.' >&2
        exit 1
        ;;
esac

case "${ID:-}:${VERSION_ID:-}" in
    centos:7*) printf '%s\n' 'Core: GCC 4.8 is supported. Helpers/tests additionally require a separately installed Python 3.8+; CentOS 7 python3 is 3.6 and is insufficient. Pass make PYTHON=/path/to/python3. See docs/PLATFORMS.md.' ;;
    *:8*) printf '%s\n' 'Use make PYTHON=python3.8 and a PATH whose python3 resolves to Python 3.8+ for installed helper commands.' ;;
esac
[ "$mode" = install ] || exit 0
if [ "$(id -u)" -ne 0 ]; then
    printf '%s\n' 'Installing system packages requires root. Run the printed commands as root,' \
        'or explicitly invoke this script with sudo and the same options.' >&2
    exit 1
fi
if ! command -v "$manager" >/dev/null 2>&1; then
    printf 'Required package manager not found: %s\n' "$manager" >&2
    exit 1
fi

if [ "$centos_vault" = yes ]; then
    if [ ! -r "$vault_key" ]; then
        printf 'Missing CentOS package signing key: %s\n' "$vault_key" >&2
        printf '%s\n' 'Restore the official centos-release package/key before using the archive.' >&2
        exit 1
    fi
    vault_reposdir=$(mktemp -d "${TMPDIR:-/tmp}/tevox-centos-repos.XXXXXX")
    trap 'rm -rf "$vault_reposdir"' 0
    trap 'exit 129' HUP
    trap 'exit 130' INT
    trap 'exit 143' TERM
    for component in $vault_components; do
        {
            printf '[tevox-centos-%s]\n' "$component"
            printf 'name=CentOS %s archive - %s\n' "$vault_version" "$component"
            printf 'baseurl=https://vault.centos.org/%s/%s/$basearch/%s\n' "$vault_version" "$component" "$vault_suffix"
            printf '%s\n' 'enabled=1' 'gpgcheck=1' 'sslverify=1' 'skip_if_unavailable=0'
            printf 'gpgkey=file://%s\n\n' "$vault_key"
        } >> "$vault_reposdir/tevox-centos-vault.repo"
    done
    "$manager" --setopt="reposdir=$vault_reposdir" install -y gcc make diffutils gawk gzip $python_package
elif [ "$manager" = apt-get ]; then
    apt-get update
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends build-essential diffutils gawk gzip $python_package
else
    "$manager" install -y gcc make diffutils gawk gzip $python_package
fi
