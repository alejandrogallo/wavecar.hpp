TEST_BINARIES = tests/reader

test: $(TEST_BINARIES)

tests/%: tests/%.cpp
	g++ -std=c++11 --pedantic -Wall -I include/ -o $@ $<

