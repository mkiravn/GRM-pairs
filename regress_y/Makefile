CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -lm
TARGET = regress_y
SOURCES = main.c

$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCES)

clean:
	rm -f $(TARGET)

.PHONY: clean
