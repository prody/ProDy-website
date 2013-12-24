 # Makefile for ProDy website

SPHINXBUILD   = sphinx-build
BUILDDIR      = _build
WORKDIR       = _workdir

.PHONY: help clean clone pull html pdf copy stats

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make HTML files of the website"
	@echo "  pdf	    to make PDF files of the manual and tutorials"
	@echo "  devel	    to make HTML/PDF files for the devel version"

clean:
	-rm -rf $(BUILDDIR)/*

clone:
	if [ ! -d "ProDy/.git" ]; then git clone https://github.com/prody/ProDy.git; fi

drugui:
	if [ ! -d "tutorials/drugui_tutorial/.git" ]; then git clone https://github.com/prody/DruGUI.git tutorials/drugui_tutorial; fi
	cd tutorials/drugui_tutorial; git pull

pull: clone
	cd ProDy; git checkout master; git pull origin master
	cd ProDy; git checkout devel; git pull origin devel

latest: pull
	cd ProDy; git checkout `git describe --tags --abbrev=0`; make build

devel: pull
	cd ProDy; git checkout devel; make build; cd docs; make build
	mv -f ProDy/docs/_build/html _build/html/devel

link:
	ln -sf ProDy/docs manual
	ln -sf ProDy/docs/_theme _theme

workdir:
	mkdir -p $(WORKDIR)
	cp -rf tutorials/conformational_sampling/$(WORKDIR)/* $(WORKDIR)
	cp -rf tutorials/enm_analysis/$(WORKDIR)/* $(WORKDIR)
	cp -rf tutorials/ensemble_analysis/$(WORKDIR)/* $(WORKDIR)
	cp -rf tutorials/evol_tutorial/$(WORKDIR)/* $(WORKDIR)
	cp -rf tutorials/prody_tutorial/$(WORKDIR)/* $(WORKDIR)
	cp -rf tutorials/structure_analysis/$(WORKDIR)/* $(WORKDIR)
	cp -rf tutorials/trajectory_analysis/$(WORKDIR)/* $(WORKDIR)

html: latest link drugui workdir
	cd $(WORKDIR); $(SPHINXBUILD) -b html -d ../$(BUILDDIR)/doctrees ../ ../$(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

pdf: latest link
	rm -f reference/*pdf tutorials/*/*pdf tutorials/*/*files.zip tutorials/*/*files.tgz

	cd manual; make pdf
	mv manual/_build/latex/ProDy.pdf _build/html/manual/

	make -C tutorials/conformational_sampling copy
	make -C tutorials/drugui_tutorial copy
	make -C tutorials/enm_analysis copy
	make -C tutorials/ensemble_analysis copy
	make -C tutorials/evol_tutorial copy
	make -C tutorials/nmwiz_tutorial copy
	make -C tutorials/prody_tutorial copy
	make -C tutorials/structure_analysis copy
	make -C tutorials/trajectory_analysis copy
