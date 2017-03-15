#!/usr/bin/env python2

#system modules
import os
import sys
import itertools as it
import argparse as ap
import logging
#external modules
import numpy as np
from requests import get #requests package
import validators        #validators package
#private modules
from cif2ice import read_cif

def shortest_distance(atoms):
    dmin = 1e99
    for a1,a2 in it.combinations(atoms,2):
        name1,x1,y1,z1 = a1
        name2,x2,y2,z2 = a2
        d = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        if d < dmin:
            dmin = d
    return dmin**0.5


# python format for GenIce.
def write_py(atoms, box, f, matchfunc=None):
    filtered = []
    if matchfunc is not None:
        for a in atoms:
            if matchfunc(a[0]):
                filtered.append(a)
    else:
        filtered = atoms
    dmin = shortest_distance(filtered)
    scale = 2.76 / dmin

    s = ""
    if (len(box) == 6):
        npbox = [v*scale for v in box]
        s += 'celltype = "triclinic"\n'
        s += 'cell = """\n{0} 0 0\n{1} {2} 0\n{3} {4} {5}\n"""\n'.format(*npbox)
        cell = np.array([[box[0],0,0],[box[1],box[2],0],[box[3],box[4],box[5]]])
        volume = np.linalg.det(cell)
    else:
        s += 'celltype = "rect"\n'
        s += 'cell = """\n{0} {1} {2}\n"""\n'.format(box[0]*scale,
                                                     box[1]*scale,
                                                     box[2]*scale)
        volume = np.product(box)
    s += 'waters = """\n'
    for name,x,y,z in filtered:
        s += "{0} {1} {2}\n".format(x*scale,y*scale,z*scale)
    s += '"""\n'
    s += 'coord = "absolute"\n'
    s += 'bondlen = 3\n'
    density = len(filtered)*18.0/(volume*scale**3*1e-24*6.022e23)
    s += 'density = {0}\n'.format(density)
    f.write(s)
    
def download(url, file_name):
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = get(url)
        # write to file
        file.write(response.content)
        
def getoptions():
    parser = ap.ArgumentParser(description='')
    parser.add_argument('--rep',  '-r', nargs = 3, type=int,   dest='rep',  default=[1,1,1],
                        help='Repeat the unit cell in x,y, and z directions. [1,1,1]')
    parser.add_argument('--debug', '-D', action='store_true', dest='debug',
                        help='Output debugging info.')
    parser.add_argument('--quiet', '-q', action='store_true', dest='quiet',
                        help='Do not output progress messages.')
    parser.add_argument('--force', '-f', action='store_true', dest='force',
                        help='Force overwrite.')
    parser.add_argument('name', nargs=1,
                       help='CIF file, Zeolite 3-letter code, or URL')
    return parser.parse_args()


def main():
    #prepare user's workarea
    home = os.path.expanduser("~")
    if os.path.exists(home+"/Library/Application Support"): #MacOS
        homegenice = home+"/Library/Application Support/GenIce"
    else:
        homegenice = os.path.expanduser(home + "/.genice") #Other unix
    sys.path.append(homegenice)
    try:
        os.makedirs(homegenice+"/lattices")
#        os.makedirs(homegenice+"/molecules")
    except:
        pass #just ignore.
    options = getoptions()
    if options.debug:
        logging.basicConfig(level=logging.DEBUG,
                            format="%(asctime)s %(levelname)s %(message)s")
    elif options.quiet:
        logging.basicConfig(level=logging.WARN,
                            format="%(asctime)s %(levelname)s %(message)s")
    else:
        #normal
        logging.basicConfig(level=logging.INFO,
                            format="%(asctime)s %(levelname)s %(message)s")
    logger = logging.getLogger()
    Nbox = [int(x) for x in options.rep]
    name = options.name[0]
    #input must be a file......too bad.
    if os.path.exists(name):
        fNameIn = name
        fNameOut = homegenice + "/lattices/" + os.path.basename(name)
        if fNameOut[-4:] in (".cif", ".CIF"):
            fNameOut = fNameOut[:-4]
        fNameOut += ".py"
    else:
        if validators.url(name):
            URL = name
            name = os.path.basename(name)
            if name[-4:] in (".cif", ".CIF"):
                name = name[:-4]
        else:
            URL = "http://www.iza-structure.org/IZA-SC/cif/"+name+".cif"
        fNameIn = homegenice + "/lattices/" + name + ".cif"
        fNameOut = homegenice + "/lattices/" + name + ".py"
        assert not os.path.exists(fNameIn) or options.force, "File exists: {0}. Use '--force' option to overwrite.".format(fNameIn)
        assert validators.url(URL)
        download(URL, fNameIn)

    logger.info("Input: {0}".format(fNameIn))
    logger.info("Output: {0}".format(fNameOut))
    if os.path.exists(fNameOut) and not options.force:
        logger.error("File exists: {0}. Use '--force' option to overwrite.".format(fNameOut))
        sys.exit(1)
    atoms, box = read_cif.read_and_process(fNameIn, Nbox, make_rect_box=False)
    fOut = open(fNameOut, "w")
    write_py(atoms, box, fOut, matchfunc=lambda x: x[0] != "O")


        
if __name__ == "__main__":
    main()
