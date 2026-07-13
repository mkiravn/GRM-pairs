# Makefile for grm_shard_tool
#
# Targets:
#   make          -- build grm_shard_tool binary
#   make clean    -- remove compiled files

CXX      = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17
BIN      = grm_shard_tool

.PHONY: all clean

all: $(BIN)

$(BIN): grm_shard_tool.cpp
	$(CXX) $(CXXFLAGS) grm_shard_tool.cpp -o $(BIN)

clean:
	rm -f $(BIN)
