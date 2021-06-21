 # Makefile for ProDy website

SPHINXBUILD   = sphinx-build
BUILDDIR      = _build
WORKDIR       = _workdir
TEMPLATEDIR   = _template

.PHONY: help clean clone pull html pdf copy stats

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make HTML files of the website"
	@echo "  pdf	    to make PDF files of the manual and tutorials"
	@echo "  devel	    to make HTML/PDF files for the devel version"

clean:
	-rm -rf $(BUILDDIR)/doctrees/*

clone:
	if [ ! -d "ProDy/.git" ]; then git clone https://github.com/prody/ProDy.git; fi

drugui:
	if [ ! -d "tutorials/drugui_tutorial/.git" ]; then git clone https://github.com/prody/DruGUI.git tutorials/drugui_tutorial; fi
	cd tutorials/drugui_tutorial; git pull

pull: clone
	git checkout master; git pull origin master
	cd ProDy; git checkout master; git pull origin master
	#cd ProDy; git checkout devel; git pull origin devel

latest: pull
	cd ProDy; git checkout `git describe --tags --abbrev=0`; make build

devel: pull
	cd ProDy; git checkout devel; make build; cd docs; make build
	mv -f ProDy/docs/_build/html _build/html/devel

link:
	ln -sf ProDy/docs manual
	ln -sf ProDy/docs/_theme _theme

workdir:
	# creates workdir (where IPython directive input and output is saved)
	mkdir -p $(WORKDIR)
	mkdir -p $(WORKDIR)/conformational_sampling_files/

	# copies required files from individual tutorials
	cp -rf tutorials/conformational_sampling/$(WORKDIR)/* $(WORKDIR)
	cp -rf tutorials/conformational_sampling/conformational_sampling_files/* $(WORKDIR)/conformational_sampling_files/

	cp -rf tutorials/enm_analysis/enm_analysis_files/* $(WORKDIR)

	cp -rf tutorials/ensemble_analysis/$(WORKDIR)/* $(WORKDIR)

	cp -rf tutorials/evol_tutorial/$(WORKDIR)/* $(WORKDIR)

	cp -rf tutorials/prody_tutorial/prody_tutorial_files/* $(WORKDIR)

	cp -rf tutorials/structure_analysis/$(WORKDIR)/* $(WORKDIR)

	cp -rf tutorials/trajectory_analysis/trajectory_analysis_files/* $(WORKDIR)

	cp -rf tutorials/comd_tutorial/comd_tutorial_files/* $(WORKDIR)

	cp -rf tutorials/membrane_anm/membrane_anm_files/* $(WORKDIR)

	cp -rf tutorials/clustenmd_tutorial/clustenmd_tutorial_files/* $(WORKDIR)

	cp -rf tutorials/signdy_tutorial/signdy_tutorial_files/* $(WORKDIR)

	cp -rf tutorials/aanm/aanm_files/* $(WORKDIR)

	cp -rf tutorials/cryoem_tutorial/cryoem_tutorial_files/* $(WORKDIR)

#	cp -rf tutorials/saxs_tutorial/saxs_tutorial_files/* $(WORKDIR)

html: link drugui workdir 
	cd $(WORKDIR); $(SPHINXBUILD) -b html -d ../$(BUILDDIR)/doctrees ../ ../$(BUILDDIR)/html
	#cd -; 
	mv -f $(BUILDDIR)/html/statistics/index.html $(BUILDDIR)/html/statistics/index.php
	cp -rf $(TEMPLATEDIR)/prody_stats.* $(BUILDDIR)/html/statistics/

	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

pdf: link
	rm -f reference/*pdf tutorials/*/*pdf tutorials/*/*files.zip tutorials/*/*files.tgz

	cd manual; make pdf
	mv manual/_build/latex/ProDy.pdf _build/html/manual/

	make -C tutorials/prody_tutorial clean copy
	make -C tutorials/nmwiz_tutorial clean copy
	make -C tutorials/evol_tutorial clean copy
	make -C tutorials/drugui_tutorial clean copy
	make -C tutorials/enm_analysis clean copy
	make -C tutorials/ensemble_analysis clean copy
	make -C tutorials/structure_analysis clean copy
	make -C tutorials/trajectory_analysis clean copy
	make -C tutorials/conformational_sampling clean copy
	make -C tutorials/comd_tutorial clean copy
	make -C tutorials/membrane_anm clean copy
	make -C tutorials/mech_stiff clean copy
	make -C tutorials/perturb_response clean copy
	make -C tutorials/signdy_tutorial clean copy
	make -C tutorials/cryoem_tutorial clean copy
	make -C tutorials/clustenmd_tutorial clean copy
