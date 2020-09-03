#!/usr/bin/env python2

# system modules
import os
import sys
import itertools as it
import argparse as ap
import logging
# external modules
import numpy as np
from math import pi
from requests import get  # requests package
import validators  # validators package
# private modules
from cif2ice import read_cif
from genice2.cell import cellvectors


def shortest_distance(atoms):
    dmin = 1e99
    for a1, a2 in it.combinations(atoms, 2):
        d = np.linalg.norm(a1 - a2)
        if d < dmin:
            dmin = d
    return dmin


def is_unique(L, pos):
    for x in L:
        d = x - pos
        d -= np.floor(d + 0.5)
        if np.dot(d, d) < 0.0000001:
            return False
    return True

# python format for GenIce.


def write_py(ratoms, box, Nbox, f, matchfunc):
    logger = logging.getLogger()

    args = []
    La, Lb, Lc, alpha, beta, gamma = box
    cell0 = cellvectors(a=La, b=Lb, c=Lc, A=alpha, B=beta, C=gamma)

    filtered = []

    for a in ratoms:
        if matchfunc(a[0]):
            rpos = np.array(a[1:])
            rpos -= np.floor(rpos)
            rpos -= np.floor(rpos)
            filtered.append(rpos)

    filtered = np.array(filtered)
    dmin = shortest_distance(filtered @ cell0)
    scale = 2.76 / dmin

    La *= scale * Nbox[0]
    Lb *= scale * Nbox[1]
    Lc *= scale * Nbox[2]

    args.append("a={0}".format(La))
    args.append("b={0}".format(Lb))
    args.append("c={0}".format(Lc))

    if alpha != 90.0:
        args.append("A={0}".format(alpha))

    if beta != 90.0:
        args.append("B={0}".format(beta))

    if gamma != 90.0:
        args.append("C={0}".format(gamma))

    s = ""
    npbox = [v * scale for v in box]
    s += "from genice.cell import cellvectors\n"
    H = "cell = cellvectors("
    s += H + (",\n" + " " * len(H)).join(args) + ")\n\n"

    cell = cell0 * scale
    volume = np.linalg.det(cell)
    icell = np.linalg.inv(cell)
    repl = np.zeros([filtered.shape[0] * Nbox[0] * Nbox[1] * Nbox[2], 3])
    N = 0
    Nbox = np.array(Nbox)
    for x in range(Nbox[0]):
        for y in range(Nbox[1]):
            for z in range(Nbox[2]):
                repl[N:N + filtered.shape[0]] = filtered + np.array([x, y, z])
                N += filtered.shape[0]
    repl /= Nbox

    s += 'waters = [\n'
    for rpos in repl:
        s += "[{0}, {1}, {2}],\n".format(*rpos)
    s += ']\n'
    s += 'coord = "relative"\n'
    s += 'bondlen = 3\n'
    density = len(filtered) * 18.0 / (volume * 1e-24 * 6.022e23)
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
    parser.add_argument('--rep', '-r', nargs=3, type=int, dest='rep', default=[1, 1, 1],
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
    options = getoptions()
    if options.debug:
        logging.basicConfig(level=logging.DEBUG,
                            format="%(asctime)s %(levelname)s %(message)s")
    elif options.quiet:
        logging.basicConfig(level=logging.WARN,
                            format="%(asctime)s %(levelname)s %(message)s")
    else:
        # normal
        logging.basicConfig(level=logging.INFO,
                            format="%(asctime)s %(levelname)s %(message)s")
    logger = logging.getLogger()
    Nbox = [int(x) for x in options.rep]
    name = options.name[0]
    # input must be a file......too bad.
    if os.path.exists(name):
        fNameIn = name
        fNameOut = os.path.basename(name)
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
            URL = "http://www.iza-structure.org/IZA-SC/cif/" + name + ".cif"
        fNameIn = name + ".cif"
        fNameOut = name + ".py"
        assert not os.path.exists(fNameIn) or options.force, "File exists: {0}. Use '--force' option to overwrite.".format(fNameIn)
        assert validators.url(URL)
        download(URL, fNameIn)

    logger.info("Input: {0}".format(fNameIn))
    logger.info("Output: {0}".format(fNameOut))
    if os.path.exists(fNameOut) and not options.force:
        logger.error("File exists: {0}. Use '--force' option to overwrite.".format(fNameOut))
        sys.exit(1)
    atoms, box = read_cif.read_and_process(fNameIn, make_rect_box=False)
    fOut = open(fNameOut, "w")
    write_py(atoms, box, Nbox, fOut, matchfunc=lambda x: x[0] != "O")


if __name__ == "__main__":
    main()
