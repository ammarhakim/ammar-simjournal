# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = cmpp
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

luarst:
	@echo "Building RST files from Lua scripts ...."
	cd source/sims; ${PYTHON} ./code/mkluarst.py; cd ../..
	cd source/sims; ${PYTHON} ./makesimindex.py; cd ../..
	cd source/sims-2/dg-maxwell; ${PYTHON} ../../sims/code/mkluarst.py; cd ../../..
	cd source/sims-2/dg-maxwell; ${PYTHON} ./makesimindex.py; cd ../../..
	cd source/sims-2/boltz-bgk; ${PYTHON} ../../sims/code/mkluarst.py; cd ../../..
	cd source/sims-2/boltz-bgk; ${PYTHON} ./makesimindex.py; cd ../../..
	cd source/sims-2/es-shock-1d; ${PYTHON} ../../sims/code/mkluarst.py; cd ../../..
	cd source/sims-2/es-shock-1d; ${PYTHON} ./makesimindex.py; cd ../../..
	cd source/sims-2/coupled-hw; ${PYTHON} ../../sims/code/mkluarst.py; cd ../../..

.PHONY: help Makefile luarst

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
