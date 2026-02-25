CC = gcc
CFLAGS = -Wall -Wextra -O3 -std=c99 -march=native -mtune=native -funroll-loops
LDFLAGS = -lm

TARGET = grm_pheno_join
SOURCE = main.c

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: clean
