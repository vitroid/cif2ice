#cif2ice
A small utility to prepare a new ice structure from a CIF file.

The source can be
1. A `.cif` file.
2. An URL to a `.cif` file.
3. The three-letter code of Zeolites.

##Installation
Use pip.

    pip install cif2ice

##Usage

    usage: cif2ice [-h] [--rep REP REP REP] [--debug] [--quiet] [--force] name
    
    positional arguments:
      name                  CIF file, Zeolite name, or URL
    
    optional arguments:
      -h, --help            show this help message and exit
      --rep REP REP REP, -r REP REP REP
                            Repeat the unit cell in x,y, and z directions. [1,1,1]
      --debug, -D           Output debugging info.
      --quiet, -q           Do not output progress messages.
      --force, -f           Force overwrite.


##Example
In any case, the generated python module will be stored in  the private folder for [GenIce](https://github.com/vitroid/GenIce) (.genice/lattices or Library/Application Support/GenIce):

1. To obtain a Zeolite RHO structure from the Zeolite DB

        cif2ice RHO

2. To generate a python module from the `foo.cif` file:

        cif2ice foo.cif
        
3. To make the python module from a remote `.cif` file:

        cif2ice http://somewhere/bar.cif
