TEST_BINARIES = tests/reader

test: $(TEST_BINARIES)

tests/%: tests/%.cpp
	g++ -std=c++11 --pedantic -Wall -I include/ -o $@ $<

include/Wavecar.hpp: main.org
	emacs -q --batch --eval "(progn (require 'org) (find-file \"$<\") (org-babel-tangle))" main.org

README.html: main.html
	mv $< $@

all: README.html
