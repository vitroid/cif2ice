#cif2ice
A small utility to prepare a new ice structure from a CIF file.

The source can be
1. A `.cif` file.
2. An URL to a `.cif` file.
3. The three-letter code of Zeolites.

##Usage:

    usage: __main__.py [-h] [--rep REP REP REP] [--debug] [--quiet] [--force] name
    
    positional arguments:
      name                  CIF file, Zeolite name, or URL
    
    optional arguments:
      -h, --help            show this help message and exit
      --rep REP REP REP, -r REP REP REP
                            Repeat the unit cell in x,y, and z directions. [1,1,1]
      --debug, -D           Output debugging info.
      --quiet, -q           Do not output progress messages.
      --force, -f           Force overwrite.

