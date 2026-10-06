#!/bin/sh
set -eu
repo=/workspace/tevox-review-20261003/early-development/v0.5-source
build=/workspace/tevox-build/cloud-v05
study=/workspace/tevox-build/study-env
if [ ! -f "$repo/include/tevox.h" ]; then
    echo 'Restore the preserved v0.5 checkout or inspect agent/v05-integration-20261005; do not overwrite the legacy checkout.' >&2
    exit 1
fi
for tool in gcc make bash python3 gzip; do command -v "$tool" >/dev/null; done
python3 -c 'import sys,gzip,zlib; assert sys.version_info >= (3,11), "Pinned cloud study dependencies require Python 3.11+; core/helper compatibility remains Python 3.8+"'
if [ ! -x "$study/bin/python" ]; then python3 -m venv "$study"; fi
"$study/bin/python" -m pip install numpy==2.3.5 openpyxl==3.1.5
cd "$repo"
make -j2 CC=gcc PYTHON="$study/bin/python" TARGET="$build/tevox" OBJDIR="$build/obj" check
make CC=gcc TARGET="$build/tevox" OBJDIR="$build/obj" PREFIX=/usr DESTDIR="$build/stage" install
"$build/stage/usr/bin/tevox" --version
TEVOX_BIN="$build/tevox" "$study/bin/python" studies/published-genomes-v05/test_inputs.py
TEVOX_BIN="$build/tevox" "$study/bin/python" studies/published-genomes-v05/test_run_integrity.py
"$study/bin/python" studies/validation/test_truth_mapping.py
"$study/bin/python" studies/validation/test_mouse_pcr_truth.py
"$study/bin/python" studies/validation/test_synthetic_metrics.py
