import os
import re
import sys
import warnings
import math as m
import numpy as np
from ase.io import read, write
from ase import Atoms

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from status import status

def rebuild(file, substitute, x = [99999, 99999], y = [99999, 99999],\
            z = [99999, 99999], loc = './', subloc = './'):
    """
        The method to replace the geometry structure of LAMMPS Input File from an XYZ File.
        rebuild(file, substitute, x = [xmin, xmax], y = [ymin, ymax], z = [zmin, zmax], loc = './', subloc = './')
        file: Your molecule system name on your .data file.
        substitute: XYZ File used to replace.
        x, y, z: Caption of the box.
        loc: Data File Location. The default is your current location.
        subloc: XYZ File Location. The default is your current location.
        Example:
            Input:
                from MCPoly.lmpset import DATAtoXYZ
                from MCPoly.lmpset import rebuild
                from MCPoly.lmpset import mould
                
                atoms = mould('Atoms1')
                atoms.cube(1, 0, 1, 0, 2, 0) # Create the file Atoms1_112.data, see in lmpset.cube.
                DATAtoXYZ('Atoms1_112') # Create Atoms1_112.xyz, see in lmpset.DATAtoXYZ
                # After changing the geometry structure by GaussView or Avogrado and save it as Atoms1_replace.xyz
                rebuild('Atoms1_112', 'Atoms1_replace')
            
            Atoms1_replace.xyz:
                106
                
                O          9.64263        4.67845        3.83021
                C          9.97863        4.59047        2.66070
                N         11.06345        5.23157        2.07395
                O          9.31893        3.82840        1.73881
                C         12.08548        6.09614        2.75528
                H         11.25704        4.90681        1.13587
                C          8.15829        2.96807        2.00836
                ......
            
            Output in Atoms1_112.data:
                ......
                
                Atoms
                
                1 1 1 -0.4804 9.64263 4.67845 3.83021
                2 1 2 0.7694 9.97863 4.59047 2.66070
                3 1 3 -1.1201 11.06345 5.23157 2.07395
                4 1 4 -0.3798 9.31893 3.82840 1.73881
                5 1 5 0.1753 12.08548 6.09614 2.75528
                6 1 6 0.5150 11.25704 4.90681 1.13587
                7 1 5 0.0718 8.15829 2.96807 2.00836
                ......
            
    """
    path = os.getcwd()
    atoms = status('{0}'.format(substitute), loc = subloc)
    pos = atoms.atom()
    xmax = -99999
    xmin = 99999
    ymax = -99999
    ymin = 99999
    zmax = -99999
    zmin = 99999
    for num in pos:
        if num[1]  >  xmax:
            xmax = num[1]
        if num[2] > ymax:
            ymax = num[2]
        if num[3] > zmax:
            zmax = num[3]
        if num[1] < xmin:
            xmin = num[1]
        if num[2] < ymin:
            ymin = num[2]
        if num[3] < zmin:
            zmin = num[3]
    atoms_number = len(pos)
    os.chdir(path)
    readfile = open(loc+'{0}.data'.format(file), 'r')
    try:
        f = open(loc+'pre.data', 'x')
    except:
        f = open(loc+'pre.data', 'w')
    c = 0
    i = 0
    for line in readfile:
        x0 = re.search(' xlo xhi', line)
        y0 = re.search(' ylo yhi', line)
        z0 = re.search(' zlo zhi', line)
        a = re.search('Atoms', line)
        if x0:
            if x == [99999, 99999]:
                f.write('{0:.5f} {1:.5f} xlo xhi\n'.format(xmin, xmax))
            else:
                f.write('{0:.5f} {1:.5f} xlo xhi\n'.format(x[0], x[1]))
        elif y0:
            if y == [99999, 99999]:
                f.write('{0:.5f} {1:.5f} ylo yhi\n'.format(ymin, ymax))
            else:
                f.write('{0:.5f} {1:.5f} ylo yhi\n'.format(y[0], y[1]))
        elif z0:
            if z == [99999, 99999]:
                f.write('{0:.5f} {1:.5f} zlo zhi\n'.format(zmin, zmax))
            else:
                f.write('{0:.5f} {1:.5f} zlo zhi\n'.format(z[0], z[1]))
        elif c == 1:
            need = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(need) >= 6:
                c = 2
            if len(need) == 7:
                f.write('{0} {1} {2} {3} {4:.5f} {5:.5f} {6:.5f}\n'.format(need[0],\
                                                                           need[1], need[2], need[3], pos[i][1], pos[i][2], pos[i][3]))
                i = i + 1
            elif len(need) == 6:
                f.write('{0} {1} {2} {3:.5f} {4:.5f} {5:.5f}\n'.format(need[0],\
                                                                       need[1], need[2], pos[i][1], pos[i][2], pos[i][3]))
                i = i + 1
        elif c == 2:
            need = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(need) == 7:
                f.write('{0} {1} {2} {3} {4:.5f} {5:.5f} {6:.5f}\n'.format(need[0],\
                                                                           need[1], need[2], need[3], pos[i][1], pos[i][2], pos[i][3]))
            elif len(need) == 6:
                f.write('{0} {1} {2} {3:.5f} {4:.5f} {5:.5f}\n'.format(need[0],\
                                                                       need[1], need[2], pos[i][1], pos[i][2], pos[i][3]))
            i = i + 1
            if atoms_number == i:
                c = 0
        elif a:
            c = 1
            f.write(line)
            f.write('\n')
        else:
            f.write(line)
    readfile.close()
    f.close()
    opath = os.getcwd()
    os.chdir(loc)
    os.system('rm {0}.data'.format(file))
    os.system('cp pre.data {0}.data'.format(file))
    os.system('rm pre.data')
    os.chdir(path)