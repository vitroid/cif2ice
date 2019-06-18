#Obtain BCT structure of Zeolite and make BCT.py module for GenIce.
sample:
	./cif2ice.x BCT
%: temp_%
	./cif2ice.x -h | python Utilities/replace.py %%usage%% "    " $< > $@
# %.rst: %.md
#	md2rst $<


test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install --index-url https://test.pypi.org/simple/ genice


install:
	./setup.py install
uninstall:
	pip uninstall cif2ice
build: README.md $(wildcard genice/*.py genice/formats/*.py genice/lattices/*.py genice/molecules/*.py)
	./setup.py sdist bdist_wheel


deploy: build
	twine upload dist/*
check:
	./setup.py check

distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info
	-rm README.rst
