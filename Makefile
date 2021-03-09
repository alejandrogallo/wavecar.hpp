TEST_BINARIES = tests/wavecar-show

test: $(TEST_BINARIES)

$(TEST_BINARIES): include/Wavecar.hpp

CXXFLAGS = -g \
          -D_GLIBCXX_ASSERTIONS \
          -D_FORTIFY_SOURCE=2 \
          -fmax-errors=1 \
          -Werror \
          -std=c++11 \
          -pedantic \
          --all-warnings \
          -Wall \
          -I include/ \

include/Wavecar.hpp:
	emacs -q --batch --eval "(progn (require 'org) (find-file \"$<\") (org-babel-tangle))" main.org
