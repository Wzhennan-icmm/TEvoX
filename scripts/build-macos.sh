#!/bin/sh
# Native Apple developer tools only; no network or package installation.
set -eu

usage() {
    cat <<'EOF'
Usage: sh scripts/build-macos.sh [--universal]

Build and test TEvoX for macOS 11.0 or newer with Apple's Command Line Tools.
The default builds the native x86_64 or arm64 architecture. --universal builds
both architectures and runs the tests on the current Mac's native slice.
Outputs: build/macos/{x86_64|arm64|universal}/tevox and build-info.txt.
Executables are ad-hoc signed locally, not Developer ID signed or notarized.
EOF
}

python=${PYTHON:-python3}
"$python" -c 'import sys; sys.exit("Python 3.8+ required" if sys.version_info < (3, 8) else 0)'

universal=no
for option in "$@"; do
    case "$option" in
        --universal) universal=yes ;;
        --help|-h) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
repo=$(CDPATH= cd "$(dirname "$0")/.." && pwd)
[ "$(uname -s)" = Darwin ] || fail 'This builder requires macOS and Apple Command Line Tools.'
os_version=$(SYSTEM_VERSION_COMPAT=0 sw_vers -productVersion)
os_major=${os_version%%.*}
case "$os_major" in ''|*[!0-9]*) fail 'Cannot determine the macOS version.' ;; esac
[ "$os_major" -ge 11 ] || fail 'Building and running this package requires macOS 11.0 or newer.'
translated=$(sysctl -in sysctl.proc_translated 2>/dev/null || :)
if [ "$translated" = 1 ]; then
    fail 'This shell is running under Rosetta. Use a native terminal, or run arch -arm64 /bin/sh scripts/build-macos.sh with the same options.'
fi
native_arch=$(uname -m)
case "$native_arch" in x86_64|arm64) ;; *) fail "Unsupported Mac architecture: $native_arch" ;; esac
command -v xcrun >/dev/null 2>&1 || fail 'Install Apple Command Line Tools with xcode-select --install.'
clang=$(xcrun --sdk macosx --find clang)
lipo=$(xcrun --sdk macosx --find lipo)
codesign=$(xcrun --sdk macosx --find codesign)
sdk=$(xcrun --sdk macosx --show-sdk-path)
[ -d "$sdk" ] || fail "macOS SDK directory not found: $sdk"

build_kind=$native_arch
[ "$universal" = no ] || build_kind=universal
base=$repo/build/macos
output=$base/$build_kind
for target in "$output/tevox" "$output/build-info.txt"; do
    if command -v git >/dev/null 2>&1 && git -C "$repo" ls-files --error-unmatch -- "$target" >/dev/null 2>&1; then
        fail "Refusing to overwrite a tracked file: $target"
    fi
    [ ! -L "$target" ] || fail "Refusing to overwrite a symbolic link: $target"
    [ ! -e "$target" ] || [ -f "$target" ] || fail "Output path is not a regular file: $target"
done
mkdir -p "$base"
work=$(mktemp -d "$base/.build.XXXXXX")
trap 'rm -rf "$work"' 0
trap 'exit 1' 1 2 15

build_slice() {
    slice_arch=$1
    slice_dir=$work/$slice_arch
    mkdir -p "$slice_dir"
    set --
    for source in "$repo"/src/*.c; do
        object=$slice_dir/$(basename "$source" .c).o
        "$clang" -arch "$slice_arch" -isysroot "$sdk" -mmacosx-version-min=11.0 \
            -D_GNU_SOURCE -D_POSIX_C_SOURCE=200809L -I"$repo/include" \
            -std=c11 -O2 -Wall -Wextra -Werror -pedantic -c "$source" -o "$object"
        set -- "$@" "$object"
    done
    "$clang" -arch "$slice_arch" -isysroot "$sdk" -mmacosx-version-min=11.0 \
        "$@" -lm -o "$slice_dir/tevox"
}

if [ "$universal" = yes ]; then
    build_slice x86_64
    build_slice arm64
    "$lipo" -create "$work/x86_64/tevox" "$work/arm64/tevox" -output "$work/tevox"
    "$lipo" -verify_arch x86_64 arm64 "$work/tevox"
else
    build_slice "$native_arch"
    cp "$work/$native_arch/tevox" "$work/tevox"
    "$lipo" -verify_arch "$native_arch" "$work/tevox"
fi
"$codesign" --force --sign - --timestamp=none "$work/tevox"
"$codesign" --verify --strict --verbose=2 "$work/tevox"
TEVOX_BIN=$work/tevox bash "$repo/tests/run_tests.sh"

revision=source
if command -v git >/dev/null 2>&1 && commit=$(git -C "$repo" rev-parse --short=12 HEAD 2>/dev/null); then
    revision=$commit
    if [ -n "$(git -C "$repo" status --porcelain --untracked-files=normal)" ]; then
        revision=$revision-dirty
    fi
fi
{
    printf 'TEvoX macOS build\n'
    printf 'Source revision: %s\n' "$revision"
    printf 'Architectures: %s\n' "$("$lipo" -archs "$work/tevox")"
    printf 'Minimum macOS: 11.0\n'
    printf 'Build host: macOS %s (%s)\n' "$os_version" "$native_arch"
    printf 'SDK: %s\n' "$sdk"
    printf 'Compiler: %s\n' "$("$clang" --version | sed -n '1p')"
    printf 'Signing: ad-hoc; no Developer ID signature or notarization\n'
    printf 'Validation: functional suite passed on native %s slice\n' "$native_arch"
    if [ "$universal" = yes ]; then
        printf 'Other slice: structure verified; execute on its native Mac for runtime validation\n'
    fi
} > "$work/build-info.txt"
mkdir -p "$output"
mv -f "$work/tevox" "$output/tevox"
mv -f "$work/build-info.txt" "$output/build-info.txt"
printf 'Built and tested: %s\n' "$output/tevox"
