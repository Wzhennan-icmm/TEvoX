# Platform builds for the v0.5 integration branch

The C11 core requires a POSIX system, a C compiler and GNU Make 3.81+.
Tests and the five Python helper commands require Python 3.8+ (standard
library only). `--gzip-output` requires the `gzip` executable on PATH.
The published-genome study additionally requires NumPy and its frozen inputs.

| System | Core compiler | Python selection |
|---|---|---|
| CentOS Linux 7.9 | distribution GCC 4.8 | install Python 3.8+ separately |
| CentOS Linux 8.5 / Rocky 8 | distribution GCC 8 | `python38`, `make PYTHON=python3.8` |
| Rocky 9 / 10 | distribution GCC | distribution Python 3 |
| Ubuntu 20.04 / 22.04 / 24.04 / 26.04 | distribution GCC | distribution Python 3 |
| macOS 11+ Intel / Apple Silicon | Apple Clang | Python 3.8+ on PATH |

Rocky 10 x86_64 requires an x86-64-v3 CPU. Linux validation uses real
distribution userspaces in containers on the host Linux kernel; it is not a
test of each distribution's kernel or boot process. A deployment target of
macOS 11 is a build setting, not proof of execution on every macOS release.
Native Intel and ARM CI jobs exercise macOS 15.

```sh
sh scripts/install-deps.sh           # inspect package commands
sh scripts/install-deps.sh --install # explicitly install packages
make -j2
make check
make PREFIX="$HOME/.local" install
```

On CentOS Linux 7/8 use `--install --centos-vault` to select the official
archived repositories with package signature and TLS verification. These
operating systems are end of life. CentOS 7's default Python 3.6 cannot run
the helpers. Our compatibility test used CPython 3.8.20, built from upstream
tag `v3.8.20`, commit `39b2f82717a69dde7212bc39b673b0f55c99e6a3`, in a
separate prefix; it did not replace system Python. Python 3.8 is itself end
of life: this records a compatibility baseline, not a maintained security
runtime. Use a maintained Python version provided for your deployment when
available. To reproduce that legacy baseline from a verified CPython source:

Install `zlib-devel` before configuring CPython, and verify that `import
gzip, zlib` works afterward. Compression readers require this standard-library
extension. The container test scratch directory must permit execution of
the deliberately failing compressor used by the regression suite.

```sh
./configure --prefix="$HOME/.local/tevox-python" --without-ensurepip
make -j2 && make install
export PATH="$HOME/.local/tevox-python/bin:$PATH"
# In the TEvoX source directory:
make PYTHON=python3.8 check
```

For Python helper entry points, `python3` on PATH must be the selected
interpreter. A private directory containing `python3 -> /path/to/python3.8`
can supply it without modifying system commands.

## macOS

Install Apple's Command Line Tools (`xcode-select --install`) and Python
3.8+. Run from a native terminal rather than Rosetta:

```sh
make macos                # native Intel or ARM, full regression
make macos-universal      # both slices, test native slice
make macos-package        # Universal archive and checksum
make macos-install MACOS_PREFIX="$HOME/.local"
```

Packages include the core, all five Python helpers, documentation and a
pair example. Scripts verify architecture, ad-hoc signatures, extracted
packages and installation paths containing spaces. They do not provide
Developer ID signing or Apple notarization. Build provenance accompanies
each package. The native CI workflow also uploads packages as artifacts.

## Large outputs

`--gzip-output` streams all 17 TSVs through `gzip -1 -n` without creating
large temporary plain TSVs. `.run.json` is committed only when every output
and compressor succeeds, and lists the actual `.tsv.gz` names. Python
readers accept compressed tables or resolve a missing `.tsv` to `.tsv.gz`.
Use a fresh output prefix when changing compression mode; readers reject
ambiguous pairs where both files exist. Compression reduces disk usage;
the graph and evidence remain in memory, so it does not remove RAM limits.
