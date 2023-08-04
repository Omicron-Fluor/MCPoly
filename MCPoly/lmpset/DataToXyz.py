import os
import re
from ase import Atoms
from ase.io import write

def DataToXyz(file, loc = './', savename = ''):
    """
    The method to change LAMMPS Data File into XYZ File. Powered by ASE.
    DataToXyz(file, loc = './', savename = '')
    file: Your molecule system name on your .data file.
    loc: File Location. The default is your current location.
    savename: Name of the saved XYZ File. The default is name of origin LAMMPS Data File.
    Example:
        Input:
            from MCPoly.lmpset import DATAtoXYZ
            DATAtoXYZ('Poly1_Chain')
        Output in Poly1_Chain.xyz:
            82
            Properties = species:S:1:pos:R:3 pbc = "F F F"
            H       34.26389000       7.81856000       4.72112000
            C       33.60940000       7.29160000       5.40582000
            H       33.81841000       6.25379000       5.63907000
            C       32.30303000       7.89744000       5.79799000
            H       31.98877000       7.45674000       6.75085000
            C       31.23940000       7.64473000       4.74277000
            
            ...
    """
    opath = os.getcwd()
    os.chdir(loc)
    elementsort = []
    elements = []
    xyzs = []
    w1 = 0
    n1 = 0
    l = []
    f = open(file + '.data', 'r')
    for line in f:
        a1 = re.search(r' atoms', line)
        m1 = re.search(r'Masses', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            num = b1.group(0)
        if n1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l1 == []:
                continue
            else:
                n1 = 2
        if n1 == 2:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l1 == []:
                n1 = 0
                continue
            mass = eval(l1[-1])
            if abs(mass - 1) < 0.5:
                elementsort.append('H')
            elif abs(mass - 2) < 0.5:
                elementsort.append('D')
            elif abs(mass - 3) < 0.5:
                elementsort.append('T')
            elif abs(mass - 4) < 0.5:
                elementsort.append('He')
            elif abs(mass - 7) < 0.5:
                elementsort.append('Li')
            elif abs(mass - 9) < 0.5:
                elementsort.append('Be')
            elif abs(mass - 10.8) < 0.5:
                elementsort.append('B')
            elif abs(mass - 12) < 0.5:
                elementsort.append('C')
            elif abs(mass - 14) < 0.5:
                elementsort.append('N')
            elif abs(mass - 16) < 0.5:
                elementsort.append('O')
            elif abs(mass - 19) < 0.5:
                elementsort.append('F')
            elif abs(mass - 22.9) < 0.5:
                elementsort.append('Na')
            elif abs(mass - 24.3) < 0.5:
                elementsort.append('Mg')
            elif abs(mass - 27) < 0.5:
                elementsort.append('Al')
            elif abs(mass - 28.1) < 0.5:
                elementsort.append('Si')
            elif abs(mass - 31) < 0.5:
                elementsort.append('P')
            elif abs(mass - 32) < 0.5:
                elementsort.append('S')
            elif abs(mass - 35.45) < 0.5:
                elementsort.append('Cl')
            elif abs(mass - 80) < 0.5:
                elementsort.append('Br')
            elif abs(mass - 127) < 0.5:
                elementsort.append('I')
            elif abs(mass - 39.1) < 0.5:
                elementsort.append('K')
            elif abs(mass - 40.1) < 0.5:
                elementsort.append('Ca')
            elif abs(mass - 47.86) < 0.2:
                elementsort.append('Ti')
            elif abs(mass - 54.93) < 0.2:
                elementsort.append('Mn')
            elif abs(mass - 55.84) < 0.2:
                elementsort.append('Fe')
            elif abs(mass - 58.93) < 0.12:
                elementsort.append('Co')
            elif abs(mass - 58.69) < 0.12:
                elementsort.append('Ni')
            elif abs(mass - 63.55) < 0.2:
                elementsort.append('Cu')
            elif abs(mass - 65.4) < 0.2:
                elementsort.append('Zn')
            elif abs(mass - 101.07) < 0.2:
                elementsort.append('Ru')
            elif abs(mass - 102.91) < 0.2:
                elementsort.append('Rh')
            elif abs(mass - 106.42) < 0.2:
                elementsort.append('Pd')
            elif abs(mass - 107.87) < 0.2:
                elementsort.append('Ag')
            elif abs(mass - 112.41) < 0.2:
                elementsort.append('Cd')
            elif abs(mass - 192.22) < 0.2:
                elementsort.append('Ir')
            elif abs(mass - 195.08) < 0.2:
                elementsort.append('Pt')
            elif abs(mass - 196.97) < 0.2:
                elementsort.append('Au')
            elif abs(mass - 200.59) < 0.2:
                elementsort.append('Hg')
        if m1:
            n1 = 1
        c1 = re.search('Atoms', line)
        if c1:
            w1 = 1
        elif w1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(l1) >= 6:
                l.append(l1)
                if len(l) == eval(num):
                    w1 = 2
        elif w1 == 2:
            for i4 in range(len(l)):
                elements.append(elementsort[eval(l[i4][2]) - 1])
                xyzs.append([eval(l[i4][4]), eval(l[i4][5]), eval(l[i4][6])])
            w1 = 0
    f.close()
    molecule = Atoms(elements, positions = xyzs)
    if savename == '':
        write(file + '.xyz', molecule)
    else:
        write(savename + '.xyz', molecule)
    os.chdir(opath)