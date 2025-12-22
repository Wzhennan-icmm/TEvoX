CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g
TARGET = tevo
SRCDIR = src
OBJDIR = obj
SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

.PHONY: all clean install test

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ -lm

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(TARGET)

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

test: $(TARGET)
	./test/run_tests.sh

help:
	@echo "Available targets:"
	@echo "  all     - Build the TEvo analyzer"
	@echo "  clean   - Remove build files"
	@echo "  install - Install to /usr/local/bin"
	@echo "  test    - Run tests"
	@echo "  help    - Show this help message"