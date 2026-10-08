#!/bin/sh
# Build, test, and package a locally ad-hoc signed macOS command-line tool.
set -eu

usage() {
    cat <<'EOF'
Usage: sh scripts/package-macos.sh [--universal] [--revision LABEL] [--output-dir DIR]

Build and test a native or universal macOS 11.0+ executable, then create a tar.gz
archive and SHA256 file. Default directory: build/macos/packages.
Default label: source commit plus -dirty when needed, or source without Git.
Labels may contain 1-64 letters, digits, dots, underscores, or hyphens.
The package is ad-hoc signed locally, not Developer ID signed or notarized.
EOF
}

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
repo=$(CDPATH= cd "$(dirname "$0")/.." && pwd)
universal=no
revision=
output_dir=$repo/build/macos/packages
while [ "$#" -gt 0 ]; do
    case "$1" in
        --universal) universal=yes; shift ;;
        --revision)
            [ "$#" -ge 2 ] || fail '--revision requires a label.'
            revision=$2
            [ -n "$revision" ] || fail '--revision cannot be empty.'
            shift 2
            ;;
        --output-dir)
            [ "$#" -ge 2 ] || fail '--output-dir requires a directory.'
            output_dir=$2
            [ -n "$output_dir" ] || fail '--output-dir cannot be empty.'
            shift 2
            ;;
        --help|-h) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done
if [ -n "$revision" ]; then
    case "$revision" in *[!A-Za-z0-9._-]*) fail 'Invalid revision label.' ;; esac
    [ "${#revision}" -le 64 ] || fail 'Revision label must be at most 64 characters.'
fi

if [ "$universal" = yes ]; then
    sh "$repo/scripts/build-macos.sh" --universal
    build_kind=universal
else
    sh "$repo/scripts/build-macos.sh"
    build_kind=$(uname -m)
fi
build_dir=$repo/build/macos/$build_kind
if [ -z "$revision" ]; then
    revision=$(sed -n 's/^Source revision: //p' "$build_dir/build-info.txt")
fi
case "$revision" in ''|*[!A-Za-z0-9._-]*) fail 'Invalid source revision in build metadata.' ;; esac
[ "${#revision}" -le 64 ] || fail 'Revision label must be at most 64 characters.'
package_name=tevox-$revision-macos-$build_kind
archive_name=$package_name.tar.gz
mkdir -p "$output_dir"
output_dir=$(CDPATH= cd "$output_dir" && pwd)
for target in "$output_dir/$archive_name" "$output_dir/$archive_name.sha256"; do
    if command -v git >/dev/null 2>&1 && git -C "$repo" ls-files --error-unmatch -- "$target" >/dev/null 2>&1; then
        fail "Refusing to overwrite a tracked file: $target"
    fi
    [ ! -L "$target" ] || fail "Refusing to overwrite a symbolic link: $target"
    [ ! -e "$target" ] || [ -f "$target" ] || fail "Output path is not a regular file: $target"
done
work=$(mktemp -d "$output_dir/.package.XXXXXX")
trap 'rm -rf "$work"' 0
trap 'exit 1' 1 2 15
payload=$work/$package_name
mkdir -p "$payload/bin" "$payload/docs" "$payload/examples"
cp "$build_dir/tevox" "$payload/bin/tevox"
chmod 755 "$payload/bin/tevox"
cp "$build_dir/build-info.txt" "$payload/docs/build-info.txt"
cp "$repo/README.md" "$payload/docs/README.md"
cp "$repo/README.md" "$payload/README.md"
for document in "$repo"/docs/*.md; do
    [ ! -f "$document" ] || cp "$document" "$payload/docs/"
done
cp -R "$repo/tests/data/pair" "$payload/examples/pair"
cp "$repo/LICENSE" "$payload/LICENSE"
for helper in phylo score_audit benchmark export_training split_audit; do
    command_name=$(printf '%s' "$helper" | tr '_' '-')
    install -m 755 "$repo/scripts/tevox_$helper.py" "$payload/bin/tevox-$command_name"
done
cat > "$payload/INSTALL.md" <<'EOF'
# TEvoX for macOS

Requires Python 3.8+ on PATH, macOS 11.0 or newer and an architecture listed in docs/build-info.txt.
Run from this extracted directory:

    ./bin/tevox --help
    ./bin/tevox pair --genome-a A --fasta-a examples/pair/A.fa --te-a examples/pair/A.gff3 \
        --genome-b B --fasta-b examples/pair/B.fa --te-b examples/pair/B.gff3 \
        --paf examples/pair/A_B.paf --flank 20 --candidate-window 10 --output example

Install for the current user without administrator privileges:

    mkdir -p "$HOME/.local/bin"
    install -m 755 bin/tevox bin/tevox-* "$HOME/.local/bin/"
    "$HOME/.local/bin/tevox" --help

Add "$HOME/.local/bin" to PATH if desired.

This executable has an ad-hoc local signature. It is not Developer ID signed
or notarized by Apple; downloading it may trigger Gatekeeper policy. The SHA256
file checks archive integrity, not publisher identity. For distribution that
requires Gatekeeper acceptance, arrange Developer ID signing and notarization.
EOF
COPYFILE_DISABLE=1 tar -czf "$work/$archive_name" -C "$work" "$package_name"
mkdir "$work/extracted"
tar -xzf "$work/$archive_name" -C "$work/extracted"
codesign=$(xcrun --sdk macosx --find codesign)
"$codesign" --verify --strict "$work/extracted/$package_name/bin/tevox"
TEVOX_BIN=$work/extracted/$package_name/bin/tevox bash "$repo/tests/run_tests.sh"
(
    cd "$work"
    shasum -a 256 "$archive_name" > "$archive_name.sha256"
    shasum -a 256 -c "$archive_name.sha256"
)
mv -f "$work/$archive_name" "$output_dir/$archive_name"
mv -f "$work/$archive_name.sha256" "$output_dir/$archive_name.sha256"
printf 'Package: %s\nSHA256: %s\n' "$output_dir/$archive_name" "$output_dir/$archive_name.sha256"
