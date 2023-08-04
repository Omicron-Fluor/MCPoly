import os
import re
from ase import Atoms
from ase.io import write

def DataToMolTxt(file, types = {0:0}, loc = './', savename = ''):
    """
    The method to change LAMMPS Data File into LAMMPS Molecule Text File.
    DataToMolTxt(file, types = {0:0, ...}, loc = './', savename = '')
    file: Your molecule system name on your .data file.
    types: Change the Atom Type number to suit the main LAMMPS input files.
    loc: File Location. The default is your current location.
    savename: Name of the saved XYZ File. The default is name of origin LAMMPS Data File.
    
    Example 1:
        Input:
            from MCPoly.lmpset import DATAtoXYZ
            DATAtomolTXT('Poly1')
            
        Output in Poly1.txt:
            #Polymer: BioPolymer1

            84 atoms
            85 bonds
            158 angles
            217 dihedrals
            48 impropers
            
            Coords
            
            1 1.3822 0.027 0.8364
                
                ...
       
            Types
    
            1 1
            2 2
            3 1
            4 2
            5 2
            6 3
            7 2
            8 4
            9 5
            
                ...
    
    Example 2:
        Input:
            from MCPoly.lmpset import DATAtomolTXT
            DataToMolTxt('Poly1', types = {1:11}, savename = 'Poly1_rami')
            
        Output in Poly1_rami:
            #Polymer: BioPolymer1

            84 atoms
            85 bonds
            158 angles
            217 dihedrals
            48 impropers
            
            Coords
            
            1 1.3822 0.027 0.8364
                
                ...
       
            Types
    
            1 11
            2 2
            3 11
            4 2
            5 2
            6 3
            7 2
            8 4
            9 5
            
                ...
    """
    opath = os.getcwd()
    os.chdir(loc)
    k = 0
    w1 = 0
    n1 = 0
    l = []
    pastline = ''
    real = [0, 0, 0]
    f = open(file+'.data', 'r')
    try:
        if savename == '':
            txt = open(file+'.txt', 'x')
        else:
            txt = open(savename+'.txt', 'x')
    except:
        if savename == '':
            txt = open(file+'.txt', 'w')
        else:
            txt = open(savename+'.txt', 'w')
    for line in f:
        del1 = re.search(r' types', line)
        del2 = re.search(r' [x-z]lo [x-z]hi', line)
        a1 = re.search(r' atoms', line)
        m1 = re.search(r'Masses', line)
        c1 = re.search(r'Atoms', line)
        #c2 = re.search(r'Bonds', line)
        #c3 = re.search(r'Angles', line)
        #c4 = re.search(r'Dihedrals', line)
        #c5 = re.search(r'Impropers', line)
        if del1:
            continue
        elif del2:
            continue
        elif pastline == line and line == '\n':
            continue
        elif a1:
            b1 = re.search(r'[0-9]+', line)
            num = b1.group(0)
            txt.write(line)
        elif c1:
            w1 = 1
        elif w1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(l1) == 7:
                l.append(l1)
                if l1[0] == num:
                    w1 = 2
        elif w1 == 2:
            txt.write('Coords\n\n')
            for i4 in range(len(l)):
                txt.write('{0} {1} {2} {3}\n'.format(i4+1, eval(l[i4][-3]),\
                                                     eval(l[i4][-2]), eval(l[i4][-1])))
            txt.write('\n')
            txt.write('Types\n\n')
            for i4 in range(len(l)):
                try:
                    ju = types[eval(l[i4][2])]
                    txt.write('{0} {1}\n'.format(i4+1, ju))
                except:
                    txt.write('{0} {1}\n'.format(i4+1, eval(l[i4][2])))
            txt.write('\n')
            txt.write('Charges\n\n')
            for i4 in range(len(l)):
                txt.write('{0} {1}\n'.format(i4+1, eval(l[i4][3])))
            txt.write('\n')
            w1 = 0
        elif n1 == 1:
            l0 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l0 == []:
                continue
            else:
                n1 = 2
        elif n1 == 2:
            l0 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l0 == []:
                n1 = 0
                continue
            else:
                continue
        elif m1:
            n1 = 1
        else:
            txt.write(line)
        pastline = line
    f.close()
    txt.close()