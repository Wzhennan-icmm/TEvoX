#!/bin/sh
# Build and install TEvoX for the current Mac without requiring sudo.
set -eu

usage() {
    cat <<'EOF'
Usage: sh scripts/install-macos.sh [--universal] [--prefix DIRECTORY]

Build and test TEvoX with Apple's tools, then install it into DIRECTORY/bin.
The default installation prefix is $HOME/.local. --universal includes both
Intel and Apple Silicon architectures. No shell profile is modified.
EOF
}

universal=no
prefix=${HOME:?HOME must be set}/.local
while [ "$#" -gt 0 ]; do
    case "$1" in
        --universal) universal=yes; shift ;;
        --prefix)
            if [ "$#" -lt 2 ] || [ -z "$2" ]; then
                printf '%s\n' 'Error: --prefix requires a directory.' >&2
                exit 2
            fi
            prefix=$2
            shift 2
            ;;
        --help|-h) usage; exit 0 ;;
        *) printf 'Error: Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ "$(uname -s)" != Darwin ]; then
    printf '%s\n' 'Error: This installer requires macOS. Use make install on Linux.' >&2
    exit 1
fi
case "$prefix" in
    /*) ;;
    *) prefix="$(pwd)/$prefix" ;;
esac
repo=$(CDPATH= cd "$(dirname "$0")/.." && pwd)
if [ "$universal" = yes ]; then
    sh "$repo/scripts/build-macos.sh" --universal
    flavor=universal
else
    sh "$repo/scripts/build-macos.sh"
    flavor=$(uname -m)
fi

binary=$repo/build/macos/$flavor/tevox
# Prepare a new executable before replacing an existing installation.
if [ -d "$prefix/bin/tevox" ]; then
    printf 'Error: Installation destination is a directory: %s/bin/tevox\n' "$prefix" >&2
    exit 1
fi
mkdir -p "$prefix/bin"
staged=$(mktemp "$prefix/bin/.tevox-install.XXXXXX")
trap 'rm -f "$staged"' 0
trap 'exit 1' 1 2 15
install -m 755 "$binary" "$staged"
"$staged" --help > /dev/null
mv -f "$staged" "$prefix/bin/tevox"
for helper in phylo score_audit benchmark export_training split_audit; do
    command_name=$(printf '%s' "$helper" | tr '_' '-')
    install -m 755 "$repo/scripts/tevox_$helper.py" "$prefix/bin/tevox-$command_name"
done
printf 'Installed: %s/bin/tevox\n'  "$prefix"
printf 'Add %s/bin to PATH if needed.\n' "$prefix"
