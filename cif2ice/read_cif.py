#!/usr/bin/env python2
#
# Python script that converts a CIF file (Crystallographic Information File)
# into a configuration file for Gromacs or LAMMPS.
#
# Copyright (C) 2016  Erik Lascaris
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#                              Written by Erik Lascaris (erikl-AT-bu-DOT-edu)
#                              Version 24-July-2014
#
# NOTE: this script uses the PyCifRW code from  https://bitbucket.org/jamesrhester/pycifrw/
#
#
# TODO:
# - Double-check that it produces a good Gromacs file, and indicate on
#   website how to use it in a Gromacs simulation.
# - Double-check that it produces a good LAMMPS file, and indicate on
#   website how to use it in a LAMMPS simulation.
#
# =============================================================================
# Modified for GenIce by Masakazu Matsumoto 2017.
#
#
#
# read_cif is written in Python3 style, while pycifrw is still in python2 style, so we assume read_cif is also executed by python2.
# Matsumoto add the __future__ feature extension here.
# =============================================================================
from __future__ import print_function, division

import sys

import logging

from math import *
import numpy as np

# Import python files from the "pycifrw-4.1.1-min" directory.
#sys.path.insert(0, 'pycifrw-4.1.1-min')
from CifFile import CifFile, CifBlock
from genice.cell import cellvectors

# =============================================================================
# =============================================================================


def float_with_error(x):
    """
    some value in cif accompanies error like "1.234(5)
    """
    if "?" in x:
        return 0
    pos = x.find("(")
    if pos >= 0:
        x = x[:pos]
    return float(x)


# =============================================================================

# Shows an error message and the usage.
def print_error(msg):
    logger = logging.getLogger()
    logger.error('{0}'.format(msg))
    print_usage()


# =============================================================================

# Converts an "_atom_type_label" into an element name.
def extract_element(label):
    logger = logging.getLogger()

    elem2 = ['He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti',
             'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
             'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
             'Sb', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
             'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt',
             'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
             'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

    if (label[0:2] in elem2):
        return label[0:2]

    elem1 = ['H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U']

    if (label[0] in elem1):
        return label[0]

    logger.warn('could not convert "%s" into element name!' % label)
    return label

# =============================================================================


# Read CIF file, and extract the necessary info in the form of a dictionary.
# E.g., the value of "_cell_volume" can be found with data['_cell_volume'].
def read_cif(fNameIn):
    logger = logging.getLogger()

    data = {}

    # Open the CIF file and read all the lines into a list of strings.
    try:
        f = open(fNameIn, 'r')
        lines = []
        oko_fix = False
        for line in f:
            stripped = line.strip()
            if (len(stripped) > 0):
                # fix for OKO
                # http://stackoverflow.com/questions/79968/split-a-string-by-spaces-preserving-quoted-substrings-in-python
                import shlex
                if stripped[0] == '_':
                    cols = shlex.split(stripped)  # shlex protects the quotation
                    if len(cols) > 2:
                        stripped = "{0} '{1}'".format(cols[0], " ".join(cols[1:]))
                        logger.debug("Line (Fixed): {0}".format(stripped))
                        oko_fix = True
                lines.append(stripped)
        if oko_fix:
            fixedfilename = fNameIn + ".fix"
            fixedfile = open(fixedfilename, "w")
            fixedfile.write("\n".join(lines) + "\n")
            fixedfile.close()
            fNameIn = fixedfilename
    except BaseException:
        logger.error("Failed to open CIF file '{0}'".format(fNameIn))
        sys.exit()

    # Use the CifFile parser to extract the data.  Although there might be
    # multiple data blocks, we'll only use the first one.
    cif_file = CifFile(fNameIn)

    for db in cif_file:
        data_block = db
        break

    try:

        # Extract some parameters, and convert them to floats.
        data['_cell_length_a'] = float_with_error(data_block['_cell_length_a'])
        data['_cell_length_b'] = float_with_error(data_block['_cell_length_b'])
        data['_cell_length_c'] = float_with_error(data_block['_cell_length_c'])
        data['_cell_angle_alpha'] = float_with_error(data_block['_cell_angle_alpha'])
        data['_cell_angle_beta'] = float_with_error(data_block['_cell_angle_beta'])
        data['_cell_angle_gamma'] = float_with_error(data_block['_cell_angle_gamma'])
        if '_cell_volume' in data_block:
            data['_cell_volume'] = float_with_error(data_block['_cell_volume'])

        # Get the symbolic operations that define the space group.  In a CIF file
        # that's the part that looks like:
        #
        # loop_
        # _symmetry_equiv_pos_as_xyz
        #   'x,y,z'
        #   'y,x,2/3-z'
        #   '-y,x-y,2/3+z'
        #   '-x,-x+y,1/3-z'
        #   '-x+y,-x,1/3+z'
        #   'x-y,-y,-z'
        #
        # In some cases it's called "_space_group_symop_operation_xyz" apparently?!?!
        data['_symmetry_equiv_pos_as_xyz'] = []

        try:
            xyz = data_block["_symmetry_equiv_pos_as_xyz"]

        except KeyError:
            try:
                xyz = data_block["_space_group_symop_operation_xyz"]
            except KeyError:
                logger.error("Missing item in CIF file: need either '_symmetry_equiv_pos_as_xyz' or '_space_group_symop_operation_xyz'.")
                sys.exit()

        # Copy the x,y,z symmetry group operations.  Remove the quotes if there
        # are any.
        for op_xyz in xyz:

            if (op_xyz[0] == '\''):
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz[1:-1])
            else:
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz)

        # Add x,y,z of the atoms to "data", but make sure to convert
        # e.g. "0.1549(8)" to "0.1549".
        data['_atom_site_label'] = data_block['_atom_site_label']

        data['_atom_site_fract_x'] = []
        for str_x in data_block['_atom_site_fract_x']:
            data['_atom_site_fract_x'].append(float_with_error(str_x))

        data['_atom_site_fract_y'] = []
        for str_y in data_block['_atom_site_fract_y']:
            data['_atom_site_fract_y'].append(float_with_error(str_y))

        data['_atom_site_fract_z'] = []
        for str_z in data_block['_atom_site_fract_z']:
            data['_atom_site_fract_z'].append(float_with_error(str_z))

    except KeyError as e:
        logger.error("Error!  Missing item in file.")
        logger.error(e)
        sys.exit()

    # print "ALL DATA:"
    # print data
    # print

    # Return the extracted data.
    return data


def read_and_process(fNameIn, make_rect_box=False):
    logger = logging.getLogger()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Read input file.

    # Make sure an input file was given.
    if (fNameIn == ''):
        print_error('no input file given.  Use:  -i filename')

    # Open the CIF file and read the data.
    data = read_cif(fNameIn)

    # Extract lengths and angles from the CIF file.
    La = float(data['_cell_length_a'])
    Lb = float(data['_cell_length_b'])
    Lc = float(data['_cell_length_c'])
    alpha = float(data['_cell_angle_alpha'])
    beta = float(data['_cell_angle_beta'])
    gamma = float(data['_cell_angle_gamma'])
    volume = 0
    if '_cell_volume' in data:
        volume = float(data['_cell_volume'])

    # Extract the symmetry operations.  This will be a list of strings such as:
    #    ['x,y,z', 'y,x,2/3-z', '-y,x-y,2/3+z', '-x,-x+y,1/3-z', ... ]
    ops = data['_symmetry_equiv_pos_as_xyz']

    # For proper evaluation, we need to convert "2/3" into "2./3", etc. to prevent
    # integer division which would turn e.g. 2/3 into 0.
    for i in range(len(ops)):
        ops[i] = ops[i].replace("0/", "0./")  # also for e.g. 10/9
        ops[i] = ops[i].replace("1/", "1./")
        ops[i] = ops[i].replace("2/", "2./")
        ops[i] = ops[i].replace("3/", "3./")
        ops[i] = ops[i].replace("4/", "4./")
        ops[i] = ops[i].replace("5/", "5./")
        ops[i] = ops[i].replace("6/", "6./")
        ops[i] = ops[i].replace("7/", "7./")
        ops[i] = ops[i].replace("8/", "8./")
        ops[i] = ops[i].replace("9/", "9./")
    #    ops[i] = ops[i].replace("/", "./")

    # Get the atom labels and coordinates.
    labels = data['_atom_site_label']
    fX = [float(s) for s in data['_atom_site_fract_x']]
    fY = [float(s) for s in data['_atom_site_fract_y']]
    fZ = [float(s) for s in data['_atom_site_fract_z']]

    # Create a list of 4-tuples, where each tuple is an atom:
    #   [ ('Si', 0.4697, 0.0, 0.0),  ('O', 0.4135, 0.2669, 0.1191),  ... ]
    atoms = [(labels[i], fX[i], fY[i], fZ[i]) for i in range(len(labels))]

    # Make sure that all atoms lie within the unit cell.  Also convert names such
    # as 'Oa1' into 'O'.
    for i in range(len(atoms)):
        (name, xn, yn, zn) = atoms[i]
        xn = (xn + 10.0) % 1.0
        yn = (yn + 10.0) % 1.0
        zn = (zn + 10.0) % 1.0
        name = extract_element(name)
        atoms[i] = (name, xn, yn, zn)

    # Update the user.
    logger.info('Loaded a CIF file with {0} atom coordinates and {1} symmetry operations.'.format(len(atoms), len(ops)))

    # Just for reference, here is a typical example of a CIF file:
    """
    _cell_length_a 4.916
    _cell_length_b 4.916
    _cell_length_c 5.4054
    _cell_angle_alpha 90
    _cell_angle_beta 90
    _cell_angle_gamma 120
    _cell_volume 113.131
    _exptl_crystal_density_diffrn      2.646
    _symmetry_space_group_name_H-M 'P 32 2 1'
    loop_
    _space_group_symop_operation_xyz
      'x,y,z'
      'y,x,2/3-z'
      '-y,x-y,2/3+z'
      '-x,-x+y,1/3-z'
      '-x+y,-x,1/3+z'
      'x-y,-y,-z'
    loop_
    _atom_site_label
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    Si   0.46970   0.00000   0.00000
    O   0.41350   0.26690   0.11910
    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Use symmetry operations to create the unit cell.

    # The CIF file consists of a few atom positions plus several "symmetry
    # operations" that indicate the other atom positions within the unit cell.  So
    # using these operations, create copies of the atoms until no new copies can be
    # made.

    # Two atoms are on top of each other if they are less than "eps" away.
    eps = 0.01  # in Angstrom

    # For each atom, apply each symmetry operation to create a new atom.
    imax = len(atoms)
    i = 0
    while (i < imax):

        label, x, y, z = atoms[i]

        for op in ops:

            # Python is awesome: calling e.g. eval('x,y,1./2+z') will convert the
            # string into a 3-tuple using the current values for x,y,z!
            logger.debug("OP: {0}".format(op))
            op = op.lower()
            xn, yn, zn = eval(op)

            # Make sure that the new atom lies within the unit cell.
            xn = (xn + 10.0) % 1.0
            yn = (yn + 10.0) % 1.0
            zn = (zn + 10.0) % 1.0

            # Check if the new position is actually new, or the same as a previous
            # atom.
            new_atom = True
            for at in atoms:
                if (abs(at[1] - xn) < eps and abs(at[2] - yn) < eps and abs(at[3] - zn) < eps):
                    new_atom = False

                    # Check that this is the same atom type.
                    if (at[0] != label):
                        print_error('invalid CIF file: atom of type %s overlaps with atom of type %s' % (at[0], label))

            # If the atom is new, add it to the list!
            if (new_atom):
                atoms.append((label, xn, yn, zn))  # add a 4-tuple

        # Update the loop iterator.
        i = i + 1
        imax = len(atoms)

    # Sort the atoms according to type alphabetically.
    atoms = sorted(atoms, key=lambda at: at[0])
    atoms.reverse()

    # Done with creating the unit cell.  Update the user.
    logger.info('Created a unit cell consisting of %d atoms.' % len(atoms))

    logger.info('Fractional coordinates:')
    for atom in atoms:
        logger.info('%10s  %.3f  %.3f  %.3f' % atom)

    # If the user wants us to create a copy of the current CIF file, with
    # additional atoms, then do that.  Note that the atoms here have *fractional*
    # coordinates!
    # if (fNameOut.endswith('.cif')):

    ##    write_cif(fNameIn, atoms, fNameOut)

    ##    logger.info('Done writing extended CIF file (%d atoms in total).' % len(atoms))
    # exit(0)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Convert the fractional coordinates into real coordinates.

    # The primitive vectors a,b,c are such that
    #
    #   cos(alpha) = b.c / |b||c|
    #   cos(beta)  = a.c / |a||c|
    #   cos(gamma) = a.b / |a||b|
    #
    # with the convention
    #
    #   a = La*xhat
    #   b = bx*xhat + by*yhat
    #   c = cx*xhat + cy*yhat + cz*zhat
    #
    # Determine the box size.
    box = (La, Lb, Lc, alpha, beta, gamma)
    cell = cellvectors(La, Lb, Lc, alpha, beta, gamma)
    V = np.linalg.det(cell)
    logging.info("Cell volume: {0} (calc.), {1} (data)".format(V, volume))
    return atoms, box
