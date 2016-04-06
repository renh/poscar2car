#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
import argparse
import datetime
import time

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
                    help='input file name, in VASP POSCAR format',
                    default='POSCAR'
)

parser.add_argument('-o', '--output',
                    help='output file name, in DMol car format',
                    default='geom.car'
)

parser.add_argument('-d', '--duplicate',
                    help='duplicate images in the a, b, and c directions, a 3-digit int (eg. 333 <=> 3x3x3)',
                    type=int, default=111
)
args = parser.parse_args()


VASP_POSCAR = args.input
CAR_file = args.output
try:
    fh = open(VASP_POSCAR)
    fh.close()
except:
    raise IOError('Failed to open file {}'.format(VASP_POSCAR))

duplicate = args.duplicate
if duplicate < 111 or duplicate > 999:
    raise ValueError("Invalid argument for duplicate")
nim = []
dump = str(duplicate)
for i in range(len(dump)):
    nim.append(int(dump[i]))
nim = np.array(nim)


# read POSCAR contents
with open(VASP_POSCAR, 'r') as fh:

    # I assume the first line contains the ion names
    l = fh.readline()
    IonTypes = l.split()

    # lattice constants
    lc = float(fh.readline())

    # lattice vectors (3x3)
    a = []
    for i in range(3):
        l = fh.readline()
        a.append([float(x) for x in l.split()])
    a = np.array(a)

    a *= lc

    # read number of ions per type
    l = fh.readline()
    IonsPerType = [int(x) for x in l.split()]
    assert(len(IonTypes) == len(IonsPerType))
    NumberIons = int(np.sum(IonsPerType))
    IonList = []
    for ion, n in zip(IonTypes, IonsPerType):
        IonList.extend([ion]*n)

    # read the next line, wheter seletive dynamics
    selective_dynamics = False
    l = fh.readline()
    l = l.strip()
    if l.startswith('s') or l.startswith('S'):
        selective_dynamics = True
        l = fh.readline()
        l = l.strip()

    # whether coordinates in Cartesian or Direct
    if l.startswith('d') or l.startswith('D'):
        coord_fmt = 'direct'
    elif l.startswith('c') or l.startswith('C'):
        coord_fmt = 'cart'
    else:
        raise ValueError("Invalid value before coordiantes")


    # read coordiantes, only coordiantes
    coord = []
    for iion in range(NumberIons):
        l = fh.readline()
        coord.append([float(x) for x in l.split()[:3]])
    coord = np.array(coord)

    # transform into cartesian if direct in POSCAR
    if coord_fmt == 'direct':
        coord = np.dot(coord, a)

    # read done. leaving other info, if any

try:
    fh = open(CAR_file, 'w')
except:
    raise IOError('Can not open file {}'.format(CAR_file))

# write car header
fh.write('!BIOSYM archive 3\n')
fh.write('PBC=ON\n')
fh.write('postcar2car generated CAR file, duplication = ({}, {}, {})\n'.format(*nim))
str_time = datetime.datetime.fromtimestamp(time.time())
fh.write(str_time.strftime('!DATE %c\n'))

# write PBC info
def length(v):
    return np.sqrt(np.sum(v*v))
la = []
for i in range(3):
    la.append(length(a[i]))
la = np.array(la)

angles = []
for i in range(3):
    n = (i+1) % 3
    m = (i+2) % 3
    angles.append(np.arccos(
        np.dot(a[n], a[m]) / (la[n] * la[m])
    ))

angles = np.array(angles)
angles = angles / (4.0*np.arctan(1.0)) * 180.0
#print(angles)
#alpha = np.arccos(
#    np.dot(a[1], a[2]) / (la[1]*la[2])
#)
#beta = np.arccos(
#    np.dot(a[2], a[0]) / (la[2]*la[0])
#)
#gamma = np.arccos(
#    np.dot(a[0], a[1]) / (la[0]*la[1])
#)
#print(alpha, beta, gamma)
pbc = "PBC{:10.4f}{:10.4f}{:10.4f}".format(*(la*nim))
pbc += "{:10.4f}{:10.4f}{:10.4f} (P1)\n".format(*angles)
fh.write(pbc)

# write coordinates
# loop over duplicated images
def write_coord(coord, ion_start_num):
    for iion in range(NumberIons):
        label = "{}{}".format(IonList[iion], iion+1+ion_start_num)
        line = "{:<5s}".format(label)
        line += "{:15.9f}{:15.9f}{:15.9f}".format(*(coord[iion]))
        line += " XXXX 1      xx      {:<2s}  0.000\n".format(IonList[iion])
        fh.write(line)

ion_start_num = 0
for i0 in range(nim[0]):
    for i1 in range(nim[1]):
        for i2 in range(nim[2]):
            R = np.zeros([3,3])
            R[0] = a[0] * i0
            R[1] = a[1] * i1
            R[2] = a[2] * i2
            trans = np.sum(R, axis=0)
            this_coord = np.zeros([NumberIons,3])
            for idir in range(3):
                this_coord[:,idir] = coord[:,idir] + trans[idir]
            write_coord(this_coord, ion_start_num)
            ion_start_num += NumberIons
            
            

fh.write('end\nend\n')
fh.close()
