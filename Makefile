PKGNAME=cif2ice

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
	pip install --index-url https://test.pypi.org/simple/ $(PKGNAME)


install:
	./setup.py install
uninstall:
	pip uninstall cif2ice
build: README.md $(wildcard cif2ice/*.py)
	./setup.py sdist bdist_wheel


deploy: build
	twine upload dist/*
check:
	./setup.py check

clean:
	-rm -rf build dist
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf *.egg-info
	-rm README.rst
