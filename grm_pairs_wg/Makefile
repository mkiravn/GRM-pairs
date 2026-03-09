# Makefile for grm_pairs
#
# Targets:
#   make          -- build grm_pairs binary
#   make clean    -- remove compiled files
#   make install  -- copy binary to ~/bin/

CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c99
LDFLAGS = -lm

SRCS = main.c calc_rel_pairs.c
OBJS = $(SRCS:.c=.o)
BIN  = grm_pairs

.PHONY: all clean install

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

main.o: main.c 
calc_rel_pairs.o: calc_rel_pairs.c 

clean:
	rm -f $(OBJS) $(BIN)

install: $(BIN)
	mkdir -p ~/bin
	cp $(BIN) ~/bin/
	@echo "Installed to ~/bin/$(BIN)"
