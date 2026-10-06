#!/bin/bash
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
