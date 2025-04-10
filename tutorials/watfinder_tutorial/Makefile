# Makefile for tutorials
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         = letter
BUILDDIR      = _build
WORKDIR       = _workdir
BASENAME      = $${PWD\#\#*/}

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS)

.PHONY: help html clean latex pdf archive setup

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  pdf        to make LaTeX files and run them through pdflatex"

clean:
	-rm -rf $(BUILDDIR)/* $(BASENAME)_files.zip $(BASENAME)_files.tgz $(BASENAME).pdf

setup:
	mkdir -p $(WORKDIR)
	echo $(BASENAME)
	cp -f $(BASENAME)_files/* $(WORKDIR)/ 2>/dev/null || :
	ln -fs ../../_theme
	ln -fs ../../_static
	mkdir -p _static/figures
	rm -f files.txt;
	files=`ls $(BASENAME)_files`; \
	if [ -z "$$files" ]; then echo "There are no required files." > files.txt; else \
	ls -hgoL $(BASENAME)_files/ | tail -n +2 | cut -b 14-128 > files.txt; fi

html: setup
	@echo "Copy files first:" ${PWD}
	cp -p -r ipynb_file ../../_build/html/tutorials/$(BASENAME)/
	cd $(WORKDIR); $(SPHINXBUILD) -b html -d ../$(BUILDDIR)/doctrees ../ ../$(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

pdf: setup
	cd $(WORKDIR); $(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) -d ../$(BUILDDIR)/doctrees ../ ../$(BUILDDIR)/latex
	@echo
	@echo "LaTeX file is ready, see $(BUILDDIR)/latex."
	@echo "Running LaTeX files through pdflatex..."
	@echo
	make -C $(BUILDDIR)/latex all-pdf
	mv -f $(BUILDDIR)/latex/$(BASENAME).pdf .
	@echo
	@echo "PDF file is ready, see $(BASENAME).pdf."

archive:
	rm -f $(BASENAME)_files.zip $(BASENAME)_files.tgz
	files=`ls $(BASENAME)_files`; \
	if [ -n "$$files" ]; then \
	tar -chzf $(BASENAME)_files.tgz $(BASENAME)_files/*; \
	zip $(BASENAME)_files.zip $(BASENAME)_files/*; \
	echo "$(BASENAME)_files.zip and $(BASENAME)_files.tgz are created."; \
	ls -hgoL $(BASENAME)_tutorial_files/ | tail -n +2 | cut -b 14-128 > files.txt; \
	else \
	echo "There are no tutorial files." > files.txt; \
	fi

build: archive pdf

copy: build
	/bin/mv -f $(BASENAME).pdf ../../_build/html/tutorials/$(BASENAME)/
	if [ -f "$(BASENAME)_files.tgz" ]; then /bin/mv -f $(BASENAME)_files.* ../../_build/html/tutorials/$(BASENAME)/; fi
