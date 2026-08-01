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
check test: $(TARGET)
	bash ./tests/run_tests.sh
asan:
	$(MAKE) clean
	$(MAKE) CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" LDFLAGS="-fsanitize=address,undefined"
	ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 bash ./tests/run_tests.sh
install: $(TARGET)
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 0755 $(TARGET) "$(DESTDIR)$(PREFIX)/bin/tevox"
	install -m 0755 scripts/tevox_phylo.py "$(DESTDIR)$(PREFIX)/bin/tevox-phylo"
clean:
	$(RM) -r build $(TARGET)
