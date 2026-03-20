CC = gcc
CFLAGS = -O2 -Wall -std=c99
LDFLAGS = -lm
TARGETS = pheno_cov_bins sort_pheno

SRCS_PHENO = main.c
OBJS_PHENO = $(SRCS_PHENO:.c=.o)

SRCS_SORT = sort_pheno.c
OBJS_SORT = $(SRCS_SORT:.c=.o)

all: $(TARGETS)

pheno_cov_bins: $(OBJS_PHENO)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

sort_pheno: $(OBJS_SORT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS_PHENO) $(OBJS_SORT) $(TARGETS)

.PHONY: all clean
