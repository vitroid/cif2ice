#Obtain BCT structure of Zeolite and make BCT.py module for GenIce.
sample:
	cif2ice BCT
%: temp_%
	cif2ice -h | python Utilities/replace.py %%usage%% "    " $< > $@
%.rst: %.md
	md2rst $<
install:
	./setup.py install
uninstall:
	pip2 uninstall cif2ice
pypi:
	make README.rst
	./setup.py check
	./setup.py sdist bdist_wheel upload
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info
	-rm README.rst
