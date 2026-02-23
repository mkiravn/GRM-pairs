# Makefile for grm_pairs
#
# Targets:
#   make          -- build grm_pairs binary
#   make clean    -- remove compiled files
#   make install  -- copy binary to ~/bin/

CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c99
LDFLAGS = -lm

SRCS = main_verbose.c calc_rel_pairs_verbose.c
OBJS = $(SRCS:.c=.o)
BIN  = grm_pairs_verbose

.PHONY: all clean install

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

main_verbose.o: main_verbose.c
calc_rel_pairs_verbose.o: calc_rel_pairs_verbose.c

clean:
	rm -f $(OBJS) $(BIN)

install: $(BIN)
	mkdir -p ~/bin
	cp $(BIN) ~/bin/
	@echo "Installed to ~/bin/$(BIN)"
