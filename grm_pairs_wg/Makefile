# Makefile for grm_pairs_wg
#
# Targets:
#   make            -- build grm_pairs
#   make clean      -- remove object files and binary
#   make install    -- copy binary to ~/bin/

CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c99
LDFLAGS = -lm

BIN  = grm_pairs
SRCS = main.c calc_rel_pairs.c
OBJS = $(SRCS:.c=.o)

.PHONY: all clean install

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(BIN)

install: $(BIN)
	mkdir -p ~/bin
	cp $(BIN) ~/bin/
	@echo "Installed to ~/bin/$(BIN)"
