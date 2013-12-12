 # Makefile for ProDy website

SPHINXBUILD   = sphinx-build
BUILDDIR      = _build
WORKDIR       = _workdir

.PHONY: help clean clone pull html pdf copy stats

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make HTML files of the website"
	@echo "  pdf	    to make PDF files of the manual and tutorials"
	@echo "  latest	    to make HTML/PDF files for the latest version"

clean:
	-rm -rf $(BUILDDIR)/*

clone:
	if [ ! -d "ProDy/.git" ]; then git clone git@github.com:prody/ProDy.git; fi

pull: clone
	cd ProDy; git pull

stable: pull
	cd ProDy; git checkout `git describe --tags --abbrev=0`

devel: pull
	cd ProDy; git checkout devel

link:
	ln -sf ProDy/docs manual
	ln -sf ProDy/docs/_theme _theme

workdir:
	mkdir -p $(WORKDIR)
	for tut in `find tutorials/ -type d -name "*_files"` ; do \
		if [ -d "$$tut" ]; then \
		 	ln -fs ../$$tut $(WORKDIR)/;\
		fi; \
	done

html: devel link workdir
	cd $(WORKDIR); $(SPHINXBUILD) -b html -d ../$(BUILDDIR)/doctrees ../ ../$(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

pdf: pull link
	rm -f reference/*pdf tutorials/*/*pdf tutorials/*/*files.zip tutorials/*/*files.tgz

	cd manual; make pdf
	mv manual/_build/latex/ProDy.pdf _build/html/manual/

	cd tutorials; \
	for tut in `find . -type d`; do \
		if [ -d "$$tut/$$tut"_files ]; then \
		 	cd $$tut; make copy; /bin/cp -rf $(WORKDIR)/* ../../$(WORKDIR); cd ..; \
		fi; \
	done

latest: pull
	cd ProDy; git checkout devel
	cd ProDy; make build
	cd ProDy/docs; make html pdf
	mkdir -p $(LATEST)
	rm -rf $(LATEST)/*
	mv $(BUILD)/* $(LATEST)
