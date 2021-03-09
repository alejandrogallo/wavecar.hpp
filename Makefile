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

EMACSQ = emacs -q --batch

include/Wavecar.hpp: README.org
	$(EMACSQ) $< \
            --eval "(require 'org)" \
            --eval "(org-babel-tangle)" \


.deps/htmlize/htmlize.el:
	mkdir -p $(@D)
	git clone https://github.com/hniksic/emacs-htmlize $(@D)

README.html: README.org .deps/htmlize/htmlize.el
	$(EMACSQ) $< \
            --load .deps/htmlize/htmlize.el \
            --eval "(require 'org)" \
            --eval "(setq org-src-fontify-natively t)" \
            --eval "(setq org-src-fontify-natively t)" \
            --eval "(setq org-src-tab-acts-natively t)" \
            --eval "(load-theme 'tsdh-light)" \
            --eval "(org-html-export-to-html)" \

index.html: README.html
	mv $< $@
