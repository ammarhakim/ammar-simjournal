# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = build
PYTHON = python

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: clean html justhtml

html:
	@echo "Building RST files from Lua scripts ...."
	cd sims; ${PYTHON} ./code/mkluarst.py; cd ..
	cd sims; ${PYTHON} ./makesimindex.py; cd ..
	cd sims-2/dg-maxwell; ${PYTHON} ../../sims/code/mkluarst.py; cd ../..
	cd sims-2/dg-maxwell; ${PYTHON} ./makesimindex.py; cd ../..
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

justhtml:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

clean:
	-rm -rf $(BUILDDIR)/*

