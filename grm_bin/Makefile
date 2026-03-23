CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra

all: grm_bin_cov

grm_bin_cov: grm_bin_cov.cpp
	$(CXX) $(CXXFLAGS) -o grm_bin_cov grm_bin_cov.cpp

clean:
	rm -f grm_bin_cov