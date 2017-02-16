#MEP == CS1
#MTN == CS2, 16
#RHO == "sIII"
#BCT == "i"
CIF_FROM_ZEOLITEDB=MEP.cif MTN.cif RHO.cif BCT.cif
PY=$(patsubst %.cif, %.py, $(CIF_FROM_ZEOLITEDB))
sample:
	cif2ice BCT
prepare:
	[ -e pycifrw ] || git clone https://bitbucket.org/jamesrhester/pycifrw.git
	cd pycifrw && python2 setup.py install


