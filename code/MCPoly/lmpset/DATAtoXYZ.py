import math as m
import numpy as np
from ase import Atoms
from ase.io import write
import os
import re

def DATAtoXYZ(file,loc='./',savename=''):
    """
    The method to change LAMMPS Data File into XYZ File. Powered by ASE.
    chain(file,loc='./',savename='')
    file: Your molecule system name on your .data file.
    loc: File Location. The default is your current location.
    savename: Name of the saved XYZ File. The default is name of origin LAMMPS Data File.
    Example:
        Input:
            from MCPoly.lmpset import DATAtoXYZ
            DATAtoXYZ('Poly1_Chain')
        Output in Poly1_Chain.xyz:
            82
            Properties=species:S:1:pos:R:3 pbc="F F F"
            H       34.26389000       7.81856000       4.72112000
            C       33.60940000       7.29160000       5.40582000
            H       33.81841000       6.25379000       5.63907000
            C       32.30303000       7.89744000       5.79799000
            H       31.98877000       7.45674000       6.75085000
            C       31.23940000       7.64473000       4.74277000
            
            ...
            84 11 79 1
            
            ...
    """
    opath=os.getcwd()
    os.chdir(loc)
    elementsort=[]
    elements=[]
    XYZs=[]
    k=0
    w1=0
    n1=0
    H=[]
    l=[]
    real=[0,0,0]
    f=open(file+'.data','r')
    for line in f:
        a1=re.search(r' atoms',line)
        m1=re.search(r'Masses',line)
        if a1:
            b1=re.search(r'[0-9]+',line)
            num=b1.group(0)
        if n1==1:
            l1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            if l1==[]:
                continue
            else:
                n1=2
        if n1==2:
            l1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            if l1==[]:
                n1=0
                continue
            if abs(eval(l1[-1])-1)<0.5:
                elementsort.append('H')
            elif abs(eval(l1[-1])-2)<0.5:
                elementsort.append('D')
            elif abs(eval(l1[-1])-3)<0.5:
                elementsort.append('T')
            elif abs(eval(l1[-1])-4)<0.5:
                elementsort.append('He')
            elif abs(eval(l1[-1])-7)<0.5:
                elementsort.append('Li')
            elif abs(eval(l1[-1])-9)<0.5:
                elementsort.append('Be')
            elif abs(eval(l1[-1])-10.8)<0.5:
                elementsort.append('B')
            elif abs(eval(l1[-1])-12)<0.5:
                elementsort.append('C')
            elif abs(eval(l1[-1])-14)<0.5:
                elementsort.append('N')
            elif abs(eval(l1[-1])-16)<0.5:
                elementsort.append('O')
            elif abs(eval(l1[-1])-19)<0.5:
                elementsort.append('F')
            elif abs(eval(l1[-1])-22.9)<0.5:
                elementsort.append('Na')
            elif abs(eval(l1[-1])-24.3)<0.5:
                elementsort.append('Mg')
            elif abs(eval(l1[-1])-27)<0.5:
                elementsort.append('Al')
            elif abs(eval(l1[-1])-28.1)<0.5:
                elementsort.append('Si')
            elif abs(eval(l1[-1])-31)<0.5:
                elementsort.append('P')
            elif abs(eval(l1[-1])-32)<0.5:
                elementsort.append('S')
            elif abs(eval(l1[-1])-35.45)<0.5:
                elementsort.append('Cl')
            elif abs(eval(l1[-1])-80)<0.5:
                elementsort.append('Br')
            elif abs(eval(l1[-1])-127)<0.5:
                elementsort.append('I')
            elif abs(eval(l1[-1])-39.1)<0.5:
                elementsort.append('K')
            elif abs(eval(l1[-1])-40.1)<0.5:
                elementsort.append('Ca')
            elif abs(eval(l1[-1])-47.86)<0.2:
                elementsort.append('Ti')
            elif abs(eval(l1[-1])-54.93)<0.2:
                elementsort.append('Mn')
            elif abs(eval(l1[-1])-55.84)<0.2:
                elementsort.append('Fe')
            elif abs(eval(l1[-1])-58.93)<0.12:
                elementsort.append('Co')
            elif abs(eval(l1[-1])-58.69)<0.12:
                elementsort.append('Ni')
            elif abs(eval(l1[-1])-63.55)<0.2:
                elementsort.append('Cu')
            elif abs(eval(l1[-1])-65.4)<0.2:
                elementsort.append('Zn')
            elif abs(eval(l1[-1])-101.07)<0.2:
                elementsort.append('Ru')
            elif abs(eval(l1[-1])-102.91)<0.2:
                elementsort.append('Rh')
            elif abs(eval(l1[-1])-106.42)<0.2:
                elementsort.append('Pd')
            elif abs(eval(l1[-1])-107.87)<0.2:
                elementsort.append('Ag')
            elif abs(eval(l1[-1])-112.41)<0.2:
                elementsort.append('Cd')
            elif abs(eval(l1[-1])-192.22)<0.2:
                elementsort.append('Ir')
            elif abs(eval(l1[-1])-195.08)<0.2:
                elementsort.append('Pt')
            elif abs(eval(l1[-1])-196.97)<0.2:
                elementsort.append('Au')
            elif abs(eval(l1[-1])-200.59)<0.2:
                elementsort.append('Hg')
        if m1:
            n1=1
        c1=re.search('Atoms',line)
        if c1:
            w1=1
        elif w1==1:
            l1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            if len(l1)>=6:
                l.append(l1)
                if l1[0]==num:
                    w1=2
        elif w1==2:
            for i4 in range(len(l)):
                elements.append(elementsort[eval(l[i4][2])-1])
                XYZs.append([eval(l[i4][-3]),eval(l[i4][-2]),eval(l[i4][-1])])
            w1=0
    f.close()
    molecule = Atoms(elements, positions=XYZs)
    if savename=='':
        write(file+'.xyz',molecule)
    else:
        write(savename+'.xyz',molecule)
    os.chdir(opath)