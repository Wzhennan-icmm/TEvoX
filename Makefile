CC ?= cc
PREFIX ?= /usr/local
CPPFLAGS ?= -D_GNU_SOURCE -D_POSIX_C_SOURCE=200809L -Iinclude
CFLAGS ?= -O2 -g
WARNFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
LDFLAGS ?=
LDLIBS = -lm
TARGET = tevox
SOURCES = $(sort $(wildcard src/*.c))
OBJECTS = $(SOURCES:src/%.c=build/%.o)
.PHONY: all clean check test asan install
all: $(TARGET)
$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $@
build/%.o: src/%.c include/tevox.h | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) -MMD -MP -c $< -o $@
build:
	mkdir -p $@
-include $(OBJECTS:.o=.d)
build/index_test: tests/index_test.c src/index.c include/tevox.h | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) $(LDFLAGS) tests/index_test.c src/index.c -o $@
build/delta_test: tests/delta_test.c src/aln_mummer_delta.c src/io.c src/util.c include/tevox.h | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) $(LDFLAGS) tests/delta_test.c src/aln_mummer_delta.c src/io.c src/util.c $(LDLIBS) -o $@
check test: $(TARGET) build/index_test build/delta_test
	./build/index_test
	./build/delta_test tests/data/conflict/A.fa tests/data/conflict/A.gff3 tests/data/conflict/B.fa tests/data/reverse/B.gff3 tests/data/v04/delta/plus.delta tests/data/v04/delta/reverse.delta
	bash ./tests/run_tests.sh
asan:
	$(MAKE) clean
	$(MAKE) CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" LDFLAGS="-fsanitize=address,undefined" all build/index_test build/delta_test
	ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 ./build/index_test
	ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 ./build/delta_test tests/data/conflict/A.fa tests/data/conflict/A.gff3 tests/data/conflict/B.fa tests/data/reverse/B.gff3 tests/data/v04/delta/plus.delta tests/data/v04/delta/reverse.delta
	ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 bash ./tests/run_tests.sh
install: $(TARGET)
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 0755 $(TARGET) "$(DESTDIR)$(PREFIX)/bin/tevox"
	install -m 0755 scripts/tevox_phylo.py "$(DESTDIR)$(PREFIX)/bin/tevox-phylo"
	install -m 0755 scripts/tevox_score_audit.py "$(DESTDIR)$(PREFIX)/bin/tevox-score-audit"
clean:
	$(RM) -r build $(TARGET)
