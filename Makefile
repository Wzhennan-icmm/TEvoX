CC ?= cc
PREFIX ?= /usr/local
CPPFLAGS ?=
PROJECT_CPPFLAGS = -D_GNU_SOURCE -D_POSIX_C_SOURCE=200809L -Iinclude
CFLAGS ?= -O2 -g
WARNFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
LDFLAGS ?=
LDLIBS ?= -lm
TARGET ?= tevox
OBJDIR ?= build
PYTHON ?= python3
BINDIR ?= $(PREFIX)/bin
MACOS_PREFIX ?= $(HOME)/.local
SOURCES = $(sort $(wildcard src/*.c))
OBJECTS = $(SOURCES:src/%.c=$(OBJDIR)/%.o)
.PHONY: all clean check test asan install uninstall check-python macos macos-universal macos-package macos-install
.DELETE_ON_ERROR:
all: $(TARGET)
$(TARGET): $(OBJECTS)
	mkdir -p "$(dir $@)"
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o "$@"
$(OBJDIR)/%.o: src/%.c include/tevox.h
	mkdir -p "$(dir $@)"
	$(CC) $(PROJECT_CPPFLAGS) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) -MMD -MP -c "$<" -o "$@"
-include $(OBJECTS:.o=.d)
$(OBJDIR)/index_test: tests/index_test.c src/index.c include/tevox.h
	mkdir -p "$(dir $@)"
	$(CC) $(PROJECT_CPPFLAGS) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) $(LDFLAGS) tests/index_test.c src/index.c -o "$@"
$(OBJDIR)/delta_test: tests/delta_test.c src/aln_mummer_delta.c src/io.c src/sha256.c src/util.c include/tevox.h include/sha256.h
	mkdir -p "$(dir $@)"
	$(CC) $(PROJECT_CPPFLAGS) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) $(LDFLAGS) tests/delta_test.c src/aln_mummer_delta.c src/io.c src/sha256.c src/util.c $(LDLIBS) -o "$@"
check-python:
	$(PYTHON) -c 'import sys; sys.exit("Python 3.8+ is required for tests and helper commands" if sys.version_info < (3, 8) else 0)'
check test: check-python $(TARGET) $(OBJDIR)/index_test $(OBJDIR)/delta_test
	"$(abspath $(OBJDIR))/index_test"
	"$(abspath $(OBJDIR))/delta_test" tests/data/conflict/A.fa tests/data/conflict/A.gff3 tests/data/conflict/B.fa tests/data/reverse/B.gff3 tests/data/v04/delta/plus.delta tests/data/v04/delta/reverse.delta
	TEVOX_BIN="$(abspath $(TARGET))" PYTHON="$(PYTHON)" bash ./tests/run_tests.sh
asan:
	ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 $(MAKE) OBJDIR=build/asan TARGET=build/asan/tevox CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" LDFLAGS="-fsanitize=address,undefined" check
install: $(TARGET)
	install -d "$(DESTDIR)$(BINDIR)"
	install -m 0755 "$(TARGET)" "$(DESTDIR)$(BINDIR)/tevox"
	install -m 0755 scripts/tevox_phylo.py "$(DESTDIR)$(BINDIR)/tevox-phylo"
	install -m 0755 scripts/tevox_score_audit.py "$(DESTDIR)$(BINDIR)/tevox-score-audit"
	install -m 0755 scripts/tevox_benchmark.py "$(DESTDIR)$(BINDIR)/tevox-benchmark"
	install -m 0755 scripts/tevox_export_training.py "$(DESTDIR)$(BINDIR)/tevox-export-training"
	install -m 0755 scripts/tevox_split_audit.py "$(DESTDIR)$(BINDIR)/tevox-split-audit"
uninstall:
	rm -f "$(DESTDIR)$(BINDIR)/tevox" "$(DESTDIR)$(BINDIR)/tevox-phylo" "$(DESTDIR)$(BINDIR)/tevox-score-audit" "$(DESTDIR)$(BINDIR)/tevox-benchmark" "$(DESTDIR)$(BINDIR)/tevox-export-training" "$(DESTDIR)$(BINDIR)/tevox-split-audit"
clean:
	$(RM) $(OBJECTS) $(OBJECTS:.o=.d) "$(TARGET)" "$(OBJDIR)/index_test" "$(OBJDIR)/delta_test"
macos:
	PYTHON="$(PYTHON)" sh scripts/build-macos.sh
macos-universal:
	PYTHON="$(PYTHON)" sh scripts/build-macos.sh --universal
macos-package:
	PYTHON="$(PYTHON)" sh scripts/package-macos.sh --universal
macos-install:
	PYTHON="$(PYTHON)" sh scripts/install-macos.sh --prefix "$(MACOS_PREFIX)"
