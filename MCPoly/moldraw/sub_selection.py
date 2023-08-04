import math as m
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.visualize import view

def sub_selection(substitution, atoms, i, state=1, place=0, strict={1:'X'}, shuffle=0, angle=-1):
    '''This function will change the hydrogen into other substitution.''' 
    
    if substitution in ('H', 'F', 'Cl', 'Br', 'I'):
        return Halogen(substitution, atoms, i, state)
    elif substitution == 'OH':
        return OH(atoms, i, state, angle)
    elif substitution in ('OMe', 'OCH3'):
        return OMe(atoms, i, state, angle)
    elif substitution == 'SH':
        return SH(atoms, i, state, angle)
    elif substitution in ('SMe', 'SCH3'):
        return SMe(atoms, i, state, angle)
    elif substitution == 'NH2':
        return NH2(atoms, i, state)
    elif substitution in ('NMe2', 'N(CH3)2'):
        return NMe2(atoms, i, state)
    elif substitution in ('CH3', 'Me'):
        return CH3(atoms, i, state)
    elif substitution == 'C2H3':
        return C2H3(atoms, i, state, shuffle)
    elif substitution in ('C2H5', 'Et', 'i-C3H7', 'i-Pr',\
                          'iPr', 't-C4H9', 't-Bu', 'tBu'):
        if substitution in ('C2H5', 'Et'):
            return CChainR(atoms, i, 1, state)
        if substitution in ('i-C3H7', 'i-Pr', 'iPr'):
            return CChainR(atoms, i, 2, state)
        if substitution in ('t-C4H9', 't-Bu', 'tBu'):
            return CChainR(atoms, i, 3, state)
    elif substitution in ('C3H7', 'n-C3H7', 'Pr', 'n-Pr', 'nPr',\
                          'C4H9', 'n-C4H9', 'nBu', 'Bu', 'n-Bu'): 
        if substitution in ('C3H7', 'n-C3H7', 'n-Pr', 'nPr', 'Pr'):
            return CChainL(atoms, i, 2, state)
        if substitution in ('C4H9', 'n-C4H9', 'n-Bu', 'nBu', 'Bu'):
            return CChainL(atoms, i, 3, state)
    elif substitution in ('CF3', 'CCl3', 'CBr3', 'CI3'):
        return trihalo(substitution, atoms, i, state)
    elif substitution == 'NO2':
        return NO2(atoms, i, state)
    elif substitution == 'CHO':
        return CHO(atoms, i, state)
    elif substitution == 'C=NH':
        return C2N(atoms, i, state, shuffle)
    elif substitution in ('COCH3', 'COMe'):
        return COMe(atoms, i, state)
    elif substitution == 'COOH':
        return COOH(atoms, i, state)
    elif substitution in ('COOMe', 'COOCH3'):
        return COOMe(atoms, i, state)
    elif substitution == 'CN':
        return CN(atoms, i, state)
    elif substitution == 'Ph':
        return Ph(atoms, i, state, shuffle)
    elif substitution in ('Py', 'pyridine', 'Pyridine'):
        return Py(atoms, i, state, place, shuffle)
    elif substitution in ('PL', 'pyrrole', 'Pyrrole'):
        return Pyrrole(atoms, i, state, place, shuffle)
    elif substitution in ('cycloC3H5', 'cycloPr', 'cyclo3'):
        return ring3(atoms, i, state, strict)
    elif substitution in ('cycloC5H9', 'cycloPe', 'cyclo5'):
        return ring5(atoms, i, state, strict, shuffle)
    elif substitution in ('cycloC6H11', 'cycloHe', 'cyclo6'):
        return ring6(atoms, i, state, strict, shuffle)
    elif substitution in ('SOMe', 'SOCH3'):
        return SOMe(atoms, i, state)
    elif substitution in ('SO2Me', 'SO2CH3'):
        return SO2_('CH3', atoms, i, state)
    elif substitution == 'SO3H':
        return SO2_('OH', atoms, i, state)
    elif substitution in ('double5', 'd5', 'D5', 'bicyclo[2.2.1]'):
        return D(5, atoms, i, state, shuffle)
    elif substitution in ('double6', 'd6', 'D6', 'bicyclo[2.2.2]'):
        return D(6, atoms, i, state, shuffle)
    elif substitution == 'Bpin':
        return Bpin(atoms, i, state, shuffle)
    else:
        assert 'Error! No such substitution!'
    
def Halogen(substitution, atoms, i, state=1):
    elements = atoms.get_chemical_symbols()
    elements[i] = substitution
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    if state == 1:
        if substitution == 'F':
            atoms.set_distance(i, Cnear, 1.3605, fix=1)
        elif substitution == 'Cl':
            atoms.set_distance(i, Cnear, 1.781, fix=1)
        elif substitution == 'Br':
            atoms.set_distance(i, Cnear, 1.945, fix=1)
        elif substitution == 'I':
            atoms.set_distance(i, Cnear, 2.210, fix=1)
        return atoms
    elif state == 2 or state == 6:
        if substitution == 'F':
            if state == 6:
                atoms.set_distance(i, Cnear, 1.354, fix=1)
            if state == 2:
                atoms.set_distance(i, Cnear, 1.340, fix=1)
        elif substitution == 'Cl':
            atoms.set_distance(i, Cnear, 1.725, fix=1)
        elif substitution == 'Br':
            atoms.set_distance(i, Cnear, 1.870, fix=1)
        elif substitution == 'I':
            atoms.set_distance(i, Cnear, 2.080, fix=1)
        return atoms

def OH(atoms, i, state=1, angle=-1):
    atoms.append('H')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-1) + [1]
    x = atoms.get_positions()
    x[-1][0] = x[-1][0] + 50
    x[-1][1] = x[-1][1] + 50
    x[-1][2] = x[-1][2] + 50
    atoms.set_positions(x)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'O'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    if state == 1:
        atoms.set_distance(i, Cnear, 1.410, fix=1)
    elif state == 2 or state == 6:
        atoms.set_distance(i, Cnear, 1.364, fix=1)
    elif state == 4:
        atoms.set_distance(i, Cnear, 1.327, fix=1)
    elif state == 164:
        atoms.set_distance(i, Cnear, 1.670, fix=1)
    else:
        atoms.set_distance(i, Cnear, 1.410, fix=1)
    atoms.set_distance(total_atoms-1, i, 0.945, fix=1)
    atoms.set_angle(Cnear, i, total_atoms-1, 103.5)
    if angle > 0:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-1, angle, mask=mask)
    if state == 2 or state == 6:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-1, 5)
    return atoms

def OMe(atoms, i, state=1, angle=-1):
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    total_atoms = len(atoms)
    atoms = OH(atoms, i, state, angle)
    total_atoms = len(atoms)
    atoms = CH3(atoms, total_atoms-1, 8)
    #print(Cnear, CnearC)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-4) + [1]*4
    angle=0
    atoms.set_distance(i, total_atoms-1, 2, mask=mask)
    atoms.set_dihedral(CnearC, Cnear, i, total_atoms-4, 105, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-4, total_atoms-1, 88, mask=mask)
    mask[-3:] = [0, 0, 0]
    atoms.set_distance(i, total_atoms-4, 1.43)
    return atoms

def SH(atoms, i, state=1, angle=-1):
    atoms.append('H')
    total_atoms = len(atoms)
    x = atoms.get_positions()
    x[-1][0] = x[-1][0] + 50
    x[-1][1] = x[-1][1] + 50
    x[-1][2] = x[-1][2] + 50
    atoms.set_positions(x)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'S'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    if state == 16:
        atoms.set_distance(i, Cnear, 2.057, fix=1)
    else:
        atoms.set_distance(i, Cnear, 1.810, fix=1)
    atoms.set_angle(Cnear, i, total_atoms-1, 101.5)
    if angle > 0:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-1, angle)
    atoms.set_distance(total_atoms-1, i, 1.336, fix=1)
    return atoms

def SMe(atoms, i, state=1, angle=-1):
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    atoms = SH(atoms, i, state, angle)
    total_atoms = len(atoms)
    atoms = CH3(atoms, total_atoms-1, 16)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-4) + [1]*4
    angle = 0
    atoms.set_dihedral(CnearC, Cnear, i, total_atoms-4, 105, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-4, total_atoms-1, 88, mask=mask)
    mask[-3:] = [0, 0, 0]
    atoms.set_distance(i, total_atoms-4, 1.82)
    return atoms

def NH2(atoms, i, state=1):
    atoms.append('H')
    atoms.append('H')
    total_atoms = len(atoms)
    x = atoms.get_positions()
    x[-2][0] = x[-2][0] + 50
    x[-2][1] = x[-2][1] + 50
    x[-2][2] = x[-2][2] + 50
    x[-1][0] = x[-1][0] - 50
    x[-1][1] = x[-1][1] - 50
    x[-1][2] = x[-1][2] - 50
    atoms.set_positions(x)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'N'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    mask = [0]*(total_atoms-2) + [1]*2
    if state == 1:
        atoms.set_distance(i, Cnear, 1.448, fix=1)
    elif state == 2 or state == 6:
        atoms.set_distance(i, Cnear, 1.340, fix=1)
    else:
        atoms.set_distance(i, Cnear, 1.448, fix=1)
    atoms.set_distance(total_atoms-2, i, 1.0100, fix=1)
    atoms.set_angle(Cnear, i, total_atoms-2, 110.5)
    if state == 1:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-2, 70, mask=mask)
    elif state == 2 or state == 6:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-2, -175, mask=mask)
    mask[-2] = 0
    atoms.set_angle(total_atoms-2, i, total_atoms-1, 110.5)
    atoms.set_distance(total_atoms-1, i, 1.0100, fix=1)
    atoms.set_angle(Cnear, i, total_atoms-1, 110.5)
    if state == 1:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-1, 179, mask=mask)
    elif state == 2 or state == 6:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-1, 70, mask=mask)
    #print(atoms.get_positions())
    return atoms

def NMe2(atoms, i):
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    atoms = NH2(atoms, i)
    total_atoms = len(atoms)
    atomA = total_atoms-2
    atomB = total_atoms-1
    mask = [0] * total_atoms
    mask[atomA] = 1
    atoms.set_distance(i, atomA, 1.448, fix=1, mask=mask)
    mask[atomA] = 0
    atoms = CH3(atoms, atomA, 7)
    total_atoms = len(atoms)
    mask = [0] * total_atoms
    mask[atomB] = 1
    atoms.set_distance(i, atomB, 1.448, fix=1, mask=mask)
    mask[atomB] = 0
    atoms = CH3(atoms, atomB, 7)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-8) + [1]*8
    atoms.set_dihedral(CnearC, Cnear, i, total_atoms-8, 95, mask=mask)
    return atoms

def C2H3(atoms, i, state=1, shuffle=0):
    atoms.append('H')
    atoms.append('C')
    atoms.append('H')
    atoms.append('H')
    Cnear = near(i, atoms)
    total_atoms = len(atoms)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'C'
    if state == 1:
        atoms.set_distance(i, Cnear, 1.5095, fix=1)
    elif state == 2:
        atoms.set_distance(i, Cnear, 1.4600, fix=1)
    elif state == 6:
        atoms.set_distance(i, Cnear, 1.4330, fix=1)
    elif state == 82:
        atoms.set_distance(i, Cnear, 1.4440, fix=1)
    elif state == 162:
        atoms.set_distance(i, Cnear, 1.7700, fix=1)
    elif state == 164:
        atoms.set_distance(i, Cnear, 1.7600, fix=1)
    else:
        atoms.set_distance(i, Cnear, 1.5095, fix=1)
    CnearH = near(Cnear, atoms, 'H', ignore=[Cnear, i])
    #print(CnearH)
    atoms.set_chemical_symbols(elements)
    mask = [0]*(total_atoms-4) + [1]*4
    x = atoms.get_positions()
    x[-4][0] = x[-4][0] + 50
    x[-4][1] = x[-4][1] + 50
    x[-4][2] = x[-4][2] + 50
    for j in range(-3, 0):
        x[j][0] = x[j][0] - 50
        x[j][1] = x[j][1] - 50
        x[j][2] = x[j][2] - 50
    atoms.set_positions(x)
    atoms.set_angle(Cnear, i, total_atoms-4, 120, mask=mask)
    atoms.set_distance(total_atoms-4, i, 1.0800, fix=1, mask=mask)
    if shuffle == 1:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-4, 5.0, mask=mask)
    else:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-4, -175.0, mask=mask)
    mask[-4] = 0
    atoms.set_angle(Cnear, i, total_atoms-3, 120, mask=mask)
    atoms.set_distance(total_atoms-3, i, 1.3400, fix=1, mask=mask)
    if shuffle == 1:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-3, -175.0, mask=mask)
    else:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-3, 5.0, mask=mask)
    mask[-3] = 0
    atoms.set_angle(i, total_atoms-3, total_atoms-2, 120, mask=mask)
    atoms.set_distance(total_atoms-2, total_atoms-3, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-3, total_atoms-2, 0.5, mask=mask)
    mask[-2] = 0
    atoms.set_angle(i, total_atoms-3, total_atoms-1, 120, mask=mask)
    atoms.set_distance(total_atoms-1, total_atoms-3, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-3, total_atoms-1, -179.5, mask=mask)
    return atoms

def CH3(atoms, i, state=1, ring=False):
    atoms.append('H')
    atoms.append('H')
    atoms.append('H')
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    total_atoms = len(atoms)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'C'
    if state == 1:
        atoms.set_distance(i, Cnear, 1.5290, fix=1)
    elif state == 2:
        atoms.set_distance(i, Cnear, 1.5095, fix=1)
    elif state == 3:
        atoms.set_distance(i, Cnear, 1.4700, fix=1)
    elif state == 6:
        atoms.set_distance(i, Cnear, 1.5105, fix=1)
    elif state == 7:
        atoms.set_distance(i, Cnear, 1.4480, fix=1)
    elif state == 8:
        atoms.set_distance(i, Cnear, 1.4100, fix=1)
    elif state == 82:
        atoms.set_distance(i, Cnear, 1.5220, fix=1)
    elif state == 16:
        atoms.set_distance(i, Cnear, 1.8100, fix=1)
    elif state == 162:
        atoms.set_distance(i, Cnear, 1.7900, fix=1)
    elif state == 164:
        atoms.set_distance(i, Cnear, 1.7700, fix=1)
    else:
        atoms.set_distance(i, Cnear, 1.5290, fix=1)
    atoms.set_chemical_symbols(elements)
    mask = [0]*(total_atoms-3) + [1]*3
    x = atoms.get_positions()
    x[-3][0] = x[-3][0] - 50
    x[-3][1] = x[-3][1] - 50
    x[-3][2] = x[-3][2] - 50
    x[-2][0] = x[-2][0] - 50
    x[-2][1] = x[-2][1] - 50
    x[-2][2] = x[-2][2] - 50
    atoms.set_positions(x)
    atoms.set_angle(Cnear, i, total_atoms-3, 109.5, mask=mask)
    atoms.set_distance(total_atoms-3, i, 1.0900, fix=1, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-3, 109.5, mask=mask)
    if ring == True:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-3, 120, mask=mask)
    else:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-3, 179, mask=mask)
    mask[-3] = 0
    atoms.set_angle(total_atoms-3, i, total_atoms-2, 104.5, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-2, 109.5, mask=mask)
    atoms.set_distance(total_atoms-2, i, 1.0900, fix=1, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-2, 109.5, mask=mask)
    if ring == True:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-2, 0.5, mask=mask)
    else:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-2, 59, mask=mask)
    mask[-2] = 0
    x = atoms.get_positions()
    x[-1][0] = x[-1][0] + 50
    x[-1][1] = x[-1][1] + 50
    x[-1][2] = x[-1][2] + 50
    atoms.set_positions(x)
    atoms.set_angle(Cnear, i, total_atoms-1, 109.5, mask=mask)
    atoms.set_angle(total_atoms-3, i, total_atoms-1, 109.5, mask=mask)
    atoms.set_distance(total_atoms-1, i, 1.0900, fix=1, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-1, 109.5, mask=mask)
    if ring == True:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-1, -120, mask=mask)
    else:
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-1, -60, mask=mask)
    return atoms

def CChainL(atoms, i, chain, state=1):
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    atoms = CH3(atoms, i, state)
    for i in range(chain):
        total_atoms = len(atoms)
        mask = [0] * total_atoms
        total_atoms = len(atoms)
        atom = total_atoms - 2
        mask[atom] = 1
        atoms.set_distance(i, atom, 1.5290, fix=1, mask=mask)
        mask[atom] = 0
        atoms = CH3(atoms, atom)
    return atoms

def CChainR(atoms, i, chain, state=1):
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    atoms = CH3(atoms, i, state)
    #view(atoms)
    total_atoms = len(atoms)
    atomA = total_atoms - 3
    atomB = total_atoms - 2
    atomC = total_atoms - 1
    mask = [0] * total_atoms
    if chain >= 1:
        mask[atomB] = 1
        atoms.set_distance(i, atomB, 1.5290, fix=1, mask=mask)
        mask[atomB] = 0
        atoms = CH3(atoms, atomB)
        if chain >= 2:
            total_atoms = len(atoms)
            mask = [0] * total_atoms
            mask[atomC] = 1
            atoms.set_distance(i, atomC, 1.5290, fix=1, mask=mask)
            mask[atomC] = 0
            atoms = CH3(atoms, atomC)
            if chain == 3:
                total_atoms = len(atoms)
                mask = [0] * total_atoms
                mask[atomA] = 1
                atoms.set_distance(i, atomA, 1.5290, fix=1, mask=mask)
                mask[atomA] = 0
                atoms = CH3(atoms, atomA)
    total_atoms = len(atoms)
    mask = [0] * total_atoms
    if chain == 1:
        mask[-6:]=[1]*6
        atoms.set_dihedral(CnearC, Cnear, i, atomB, 110, mask=mask)
    if chain == 2:
        mask[-9:]=[1]*9
        atoms.set_dihedral(CnearC, Cnear, i, atomB, 95, mask=mask)
    return atoms

def trihalo(substitution, atoms, i, state=1):
    import re
    atoms = CH3(atoms, i, state)
    halo=re.search(r'[A-Z]+[a-z]*', substitution)
    for i in range(-3, 0):
        Halogen(halo.group(0)[1:], atoms, i)
    return atoms

def NO2(atoms, i, state=1):
    atoms.append('O')
    atoms.append('O')
    total_atoms = len(atoms)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'N'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    CnearH = near(Cnear, atoms, 'H')
    if state == 1:
        atoms.set_distance(i, Cnear, 1.4900, fix=1)
    elif state == 6:
        atoms.set_distance(i, Cnear, 1.4600, fix=1)
    atoms.set_dihedral(CnearH, Cnear, i, total_atoms-2, -175.0, mask=mask)
    atoms.set_distance(total_atoms-2, i, 1.2250, fix=1)
    atoms.set_angle(Cnear, i, total_atoms-2, 119.5)
    atoms.set_angle(total_atoms-2, i, total_atoms-1, 119.5)
    atoms.set_distance(total_atoms-1, i, 1.2250, fix=1)
    atoms.set_angle(Cnear, i, total_atoms-1, 120.5)
    #print(atoms.get_positions())
    return atoms

def CHO(atoms, i, state=1):
    atoms.append('O')
    atoms.append('H')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-2) + [1]*2
    elements = atoms.get_chemical_symbols()
    elements[i] = 'C'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    CnearH = near(Cnear, atoms, 'H')
    if state == 1:
        atoms.set_distance(i, Cnear, 1.5220, fix=1)
    elif state == 2:
        atoms.set_distance(i, Cnear, 1.4440, fix=1)
    elif state == 6:
        atoms.set_distance(i, Cnear, 1.4900, fix=1)
    elif state == 7:
        atoms.set_distance(i, Cnear, 1.3350, fix=1)
    elif state == 8:
        atoms.set_distance(i, Cnear, 1.3640, fix=1)
    elif state == 82:
        atoms.set_distance(i, Cnear, 1.5100, fix=1)
    x = atoms.get_positions()
    #Separate into a side
    x[-2][0] = x[-2][0] - 50
    x[-2][1] = x[-2][1] - 50
    x[-2][2] = x[-2][2] - 50
    atoms.set_positions(x)
    atoms.set_angle(Cnear, i, total_atoms-2, 120.5, mask=mask)
    atoms.set_distance(total_atoms-2, i, 1.2290, fix=1, mask=mask)
    atoms.set_dihedral(CnearH, Cnear, i, total_atoms-2, 55, mask=mask)
    mask[-2] = 0
    x = atoms.get_positions()
    #Separate into another side
    x[-1][0] = x[-1][0] + 50
    x[-1][1] = x[-1][1] + 50
    x[-1][2] = x[-1][2] + 50
    atoms.set_positions(x)
    for k in range(10):
        atoms.set_angle(total_atoms-2, i, total_atoms-1, 120, mask=mask)
        atoms.set_angle(Cnear, i, total_atoms-1, 120, mask=mask)
    atoms.set_distance(total_atoms-1, i, 1.0900, fix=1, mask=mask)
    return atoms

def C2N(atoms, i, state=1, shuffle=0):
    atoms.append('N')
    atoms.append('H')
    atoms.append('H')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-3) + [1]*3
    elements = atoms.get_chemical_symbols()
    elements[i] = 'C'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    CnearH = near(Cnear, atoms, 'H', ignore=[Cnear, i])
    if state == 1:
        atoms.set_distance(i, Cnear, 1.5220, fix=1)
    elif state == 2:
        atoms.set_distance(i, Cnear, 1.4440, fix=1)
    elif state == 6:
        atoms.set_distance(i, Cnear, 1.5220, fix=1)
    elif state == 7:
        atoms.set_distance(i, Cnear, 1.3350, fix=1)
    elif state == 8:
        atoms.set_distance(i, Cnear, 1.3640, fix=1)
    elif state == 82:
        atoms.set_distance(i, Cnear, 1.5100, fix=1)
    x = atoms.get_positions()
    #Separate into a side
    x[-3][0] = x[-3][0] - 50
    x[-3][1] = x[-3][1] - 50
    x[-3][2] = x[-3][2] - 50
    atoms.set_positions(x)
    atoms.set_angle(Cnear, i, total_atoms-3, 120.5, mask=mask)
    atoms.set_distance(total_atoms-3, i, 1.3580, fix=1, mask=mask)
    if shuffle == 0:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-3, 55, mask=mask)
    elif shuffle == 1:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-3, 179.5, mask=mask)
    mask[-3] = 0
    x = atoms.get_positions()
    x[-2][0] = x[-2][0] + 50
    x[-2][1] = x[-2][1] + 50
    x[-2][2] = x[-2][2] + 50
    atoms.set_positions(x)
    atoms.set_angle(i, total_atoms-3, total_atoms-2, 120.5, mask=mask)
    atoms.set_distance(total_atoms-2, total_atoms-3, 1.0100, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-3, total_atoms-2, 179.5, mask=mask)
    mask[-2] = 0
    x = atoms.get_positions()
    #Separate into another side
    x[-1][0] = x[-1][0] + 50
    x[-1][1] = x[-1][1] + 50
    x[-1][2] = x[-1][2] + 50
    atoms.set_positions(x)
    if shuffle == 1:
        atoms.set_dihedral(Cnear, i, total_atoms-3, total_atoms-1, 55, mask=mask)
    for k in range(10):
        atoms.set_angle(total_atoms-3, i, total_atoms-1, 120, mask=mask)
        atoms.set_angle(Cnear, i, total_atoms-1, 120, mask=mask)
    atoms.set_distance(total_atoms-1, i, 1.0100, fix=1, mask=mask)
    return atoms

def COMe(atoms, i, state=1):
    atoms = CHO(atoms, i, state)
    total_atoms = len(atoms)
    atoms = CH3(atoms, total_atoms-1, 82)
    total_atoms = len(atoms)
    return atoms
    
def COOH(atoms, i, state=1):
    atoms.append('O')
    atoms.append('O')
    atoms.append('H')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-3) + [1]*3
    elements = atoms.get_chemical_symbols()
    elements[i] = 'C'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    CnearH = near(Cnear, atoms, 'H')
    if state == 1:
        atoms.set_distance(i, Cnear, 1.5220, fix=1)
    elif state == 2:
        atoms.set_distance(i, Cnear, 1.4440, fix=1)
    elif state == 6:
        atoms.set_distance(i, Cnear, 1.4900, fix=1)
    elif state == 7:
        atoms.set_distance(i, Cnear, 1.3350, fix=1)
    elif state == 8:
        atoms.set_distance(i, Cnear, 1.3640, fix=1)
    elif state == 82:
        atoms.set_distance(i, Cnear, 1.5095, fix=1)
    x = atoms.get_positions()
    #Separate into a side
    x[-3][0] = x[-3][0] - 50
    x[-3][1] = x[-3][1] - 50
    x[-3][2] = x[-3][2] - 50
    atoms.set_positions(x)
    atoms.set_angle(Cnear, i, total_atoms-3, 120.5, mask=mask)
    atoms.set_distance(total_atoms-3, i, 1.3640, fix=1, mask=mask)
    atoms.set_dihedral(CnearH, Cnear, i, total_atoms-3, 10.4, mask=mask)
    mask[-3] = 0
    x = atoms.get_positions()
    #Separate into another side
    x[-2][0] = x[-2][0] + 50
    x[-2][1] = x[-2][1] + 50
    x[-2][2] = x[-2][2] + 50
    atoms.set_positions(x)
    for k in range(10):
        atoms.set_angle(total_atoms-3, i, total_atoms-2, 120, mask=mask)
        atoms.set_angle(Cnear, i, total_atoms-2, 120, mask=mask)
    atoms.set_distance(total_atoms-2, i, 1.2290, fix=1, mask=mask)
    mask[-2] = 0
    atoms.set_distance(total_atoms-1, total_atoms-3, 0.945, fix=1)
    atoms.set_angle(i, total_atoms-3, total_atoms-1, 120.5)
    atoms.set_dihedral(total_atoms-2, i, total_atoms-3, total_atoms-1, 5)
    #print(atoms.get_positions())
    #print(atoms.get_chemical_symbols())
    return atoms

def COOMe(atoms, i, state=1):
    atoms = COOH(atoms, i, state)
    total_atoms = len(atoms)
    atoms = CH3(atoms, total_atoms-1, 8)
    return atoms

def CN(atoms, i, state=1):
    atoms.append('N')
    total_atoms = len(atoms)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'C'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    if state == 1:
        atoms.set_distance(i, Cnear, 1.4700, fix=1)
    elif state == 6:
        atoms.set_distance(i, Cnear, 1.4510, fix=1)
    atoms.set_distance(total_atoms-1, i, 1.157, fix=1)
    atoms.set_angle(Cnear, i, total_atoms-1, 179.5)
    #print(atoms.get_positions())
    return atoms

def Ph(atoms, i, state=1, shuffle=0):
    for j in range(5):
        atoms.append('C')
    for j in range(5):
        atoms.append('H')
    Cnear = near(i, atoms)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-10) + [1]*10
    elements = atoms.get_chemical_symbols()
    elements[i] = 'C'
    atoms.set_chemical_symbols(elements)
    x = atoms.get_positions()
    #Separate into a side
    for j in range(-10, -8):
        x[j][0] = x[j][0] - 50
        x[j][1] = x[j][1] - 50
        x[j][2] = x[j][2] - 50
    for j in range(-5, -3):
        x[j][0] = x[j][0] - 50
        x[j][1] = x[j][1] - 50
        x[j][2] = x[j][2] - 50
    for j in range(-8, -5):
        x[j][0] = x[j][0] + 50
        x[j][1] = x[j][1] + 50
        x[j][2] = x[j][2] + 50
    for j in range(-3, 0):
        x[j][0] = x[j][0] + 50
        x[j][1] = x[j][1] + 50
        x[j][2] = x[j][2] + 50
    atoms.set_positions(x)
    #if shuffle == 0:
    CnearH = near(Cnear, atoms, 'H', ignore=[Cnear, i])
    #elif shuffle == 1:
    #    CnearH = near(Cnear, atoms, 'C')
    #print(Cnear, CnearH)
    if state == 1:
        atoms.set_distance(i, Cnear, 1.5105, fix=1)
    elif state == 2:
        atoms.set_distance(i, Cnear, 1.4330, fix=1)
    elif state == 6:
        atoms.set_distance(i, Cnear, 1.5100, fix=1)
    elif state == 7:
        atoms.set_distance(i, Cnear, 1.3400, fix=1)
    elif state == 8:
        atoms.set_distance(i, Cnear, 1.3640, fix=1)
    elif state == 16:
        atoms.set_distance(i, Cnear, 1.7700, fix=1)
    elif state == 82:
        atoms.set_distance(i, Cnear, 1.4900, fix=1)
    atoms.set_angle(Cnear, i, total_atoms-10, 120, mask=mask)
    atoms.set_distance(total_atoms-10, i, 1.4, fix=1, mask=mask)
    if shuffle == 0:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-10, 5.4, mask=mask)
    elif shuffle == 1:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-10, 95.4, mask=mask)
    mask[-10] = 0
    atoms.set_angle(Cnear, i, total_atoms-6, 120, mask=mask)
    atoms.set_distance(total_atoms-6, i, 1.4, fix=1, mask=mask)
    if shuffle == 0:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-6, -175.0, mask=mask)
    elif shuffle == 1:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-6, -85.0, mask=mask)
    mask[-6] = 0
    atoms.set_angle(i, total_atoms-10, total_atoms-9, 120, mask=mask)
    atoms.set_distance(total_atoms-9, total_atoms-10, 1.4, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-10, total_atoms-9, 179.5, mask=mask)
    mask[-9] = 0
    atoms.set_angle(i, total_atoms-6, total_atoms-7, 120, mask=mask)
    atoms.set_distance(total_atoms-7, total_atoms-6, 1.4, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-6, total_atoms-7, 179.5, mask=mask)
    mask[-7] = 0
    atoms.set_angle(total_atoms-6, total_atoms-7, total_atoms-8, 120, mask=mask)
    atoms.set_distance(total_atoms-8, total_atoms-7, 1.4, fix=1, mask=mask)
    atoms.set_angle(total_atoms-7, total_atoms-8, total_atoms-9, 120, mask=mask)
    atoms.set_dihedral(i, total_atoms-6, total_atoms-7, total_atoms-8,\
                       0.5, mask=mask)
    mask[-8] = 0
    atoms.set_angle(i, total_atoms-10, total_atoms-5, 120, mask=mask)
    atoms.set_distance(total_atoms-5, total_atoms-10, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-10, total_atoms-5, 0.5, mask=mask)
    mask[-5] = 0
    atoms.set_angle(i, total_atoms-6, total_atoms-1, 120, mask=mask)
    atoms.set_distance(total_atoms-1, total_atoms-6, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-6, total_atoms-1, 0.5, mask=mask)
    mask[-1] = 0
    atoms.set_angle(total_atoms-10, total_atoms-9, total_atoms-4, 120, mask=mask)
    atoms.set_distance(total_atoms-4, total_atoms-9, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-5, total_atoms-10, total_atoms-9,\
                       total_atoms-4, 0.5, mask=mask)
    mask[-4] = 0
    atoms.set_angle(total_atoms-6, total_atoms-7, total_atoms-2, 120, mask=mask)
    atoms.set_distance(total_atoms-2, total_atoms-7, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-1, total_atoms-6, total_atoms-7,\
                       total_atoms-2, 0.5, mask=mask)
    mask[-2] = 0
    atoms.set_angle(total_atoms-9, total_atoms-8, total_atoms-3, 120, mask=mask)
    atoms.set_distance(total_atoms-3, total_atoms-8, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-4, total_atoms-9, total_atoms-8,\
                       total_atoms-3, 0.5, mask=mask)
    return atoms

def Py(atoms, i, state=1, place=4, shuffle=0):
    atoms = Ph(atoms, i, state, shuffle)
    total_atoms = len(atoms)
    mask = [0] * total_atoms

def Pyrrole(atoms, i, state=1, place=0, shuffle=0):
    if place == 0:
        for j in range(4):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear = near(i, atoms)
        total_atoms = len(atoms)
        mask = [0]*(total_atoms-8) + [1]*8
        elements = atoms.get_chemical_symbols()
        elements[i] = 'N'
        if shuffle == 0:
            CnearH = near(Cnear, atoms, 'H')
        elif shuffle == 1:
            CnearH = near(Cnear, atoms, 'C')
        atoms.set_chemical_symbols(elements)
        if state == 1:
            atoms.set_distance(i, Cnear, 1.475, fix=1, ignore=i)
        elif state == 6:
            atoms.set_distance(i, Cnear, 1.440, fix=1, ignore=i)
        x = atoms.get_positions()
        #Separate into a side
        for j in range(-8, -6):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-4, -2):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-6, -4):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        for j in range(-2, 0):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        atoms.set_positions(x)
        atoms.set_angle(Cnear, i, total_atoms-8, 125.35, mask=mask)
        atoms.set_distance(total_atoms-8, i, 1.381, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-8, 5.4, mask=mask)
        mask[-8] = 0
        atoms.set_angle(Cnear, i, total_atoms-5, 125.35, mask=mask)
        atoms.set_distance(total_atoms-5, i, 1.381, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-5, -175.0, mask=mask)
        mask[-5] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-7, 107.7, mask=mask)
        atoms.set_distance(total_atoms-7, total_atoms-8, 1.367, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-7, 179.5, mask=mask)
        mask[-7] = 0
        atoms.set_angle(i, total_atoms-5, total_atoms-6, 107.7, mask=mask)
        atoms.set_distance(total_atoms-6, total_atoms-5, 1.367, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-6, -179.5, mask=mask)
        mask[-6] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-4, 126, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-8, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-4, 0.5, mask=mask)
        mask[-4] = 0
        atoms.set_angle(i, total_atoms-5, total_atoms-1, 126, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-5, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-1, -0.5, mask=mask)
        mask[-1] = 0
        atoms.set_angle(total_atoms-8, total_atoms-7, total_atoms-3, 126, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-7, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-8, total_atoms-7,\
                           total_atoms-3, 179.5, mask=mask)
        mask[-3] = 0
        atoms.set_angle(total_atoms-5, total_atoms-6, total_atoms-2, 126, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-6, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-5, total_atoms-6,\
                           total_atoms-2, -179.5, mask=mask)
        mask[-2] = 0
    elif place == 1:
        atoms.append('N')
        for j in range(3):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear = near(i, atoms)
        total_atoms = len(atoms)
        mask = [0]*(total_atoms-8) + [1]*8
        elements = atoms.get_chemical_symbols()
        elements[i] = 'C'
        if shuffle == 0:
            CnearH = near(Cnear, atoms, 'H', ignore=i)
        elif shuffle == 1:
            CnearH = near(Cnear, atoms, 'C', ignore=i)
        atoms.set_chemical_symbols(elements)
        if state == 1:
            atoms.set_distance(i, Cnear, 1.5100, fix=1)
        elif state == 6:
            atoms.set_distance(i, Cnear, 1.4600, fix=1)
        x = atoms.get_positions()
        #Separate into a side
        for j in range(-8, -6):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-4, -2):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-6, -4):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        for j in range(-2, 0):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        atoms.set_positions(x)
        atoms.set_angle(Cnear, i, total_atoms-8, 120.8, mask=mask)
        atoms.set_distance(total_atoms-8, i, 1.381, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-8, 6.0, mask=mask)
        mask[-8] = 0
        atoms.set_angle(Cnear, i, total_atoms-5, 131.58, mask=mask)
        atoms.set_distance(total_atoms-5, i, 1.367, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-5, -175.0, mask=mask)
        mask[-5] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-7, 109.8, mask=mask)
        atoms.set_distance(total_atoms-7, total_atoms-8, 1.381, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-7, 179.5, mask=mask)
        mask[-7] = 0
        atoms.set_angle(i, total_atoms-5, total_atoms-6, 107.3, mask=mask)
        atoms.set_distance(total_atoms-6, total_atoms-5, 1.424, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-6, -179.5, mask=mask)
        mask[-6] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-4, 126, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-8, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-4, 0.5, mask=mask)
        mask[-4] = 0
        atoms.set_angle(i, total_atoms-5, total_atoms-1, 126, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-5, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-1, -0.5, mask=mask)
        mask[-1] = 0
        atoms.set_angle(total_atoms-8, total_atoms-7, total_atoms-3, 126, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-7, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-8, total_atoms-7, total_atoms-3,\
                           179.5, mask=mask)
        mask[-3] = 0
        atoms.set_angle(total_atoms-5, total_atoms-6, total_atoms-2, 126, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-6, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-5, total_atoms-6, total_atoms-2,\
                           -179.5, mask=mask)
        mask[-2] = 0
    elif place == -1:
        atoms.append('N')
        for j in range(3):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear = near(i, atoms)
        total_atoms = len(atoms)
        mask = [0]*(total_atoms-8) + [1]*8
        elements = atoms.get_chemical_symbols()
        elements[i] = 'C'
        if shuffle == 0:
            CnearH = near(Cnear, atoms, 'H', ignore=i)
        elif shuffle == 1:
            CnearH = near(Cnear, atoms, 'C', ignore=i)
        atoms.set_chemical_symbols(elements)
        if state == 1:
            atoms.set_distance(i, Cnear, 1.5100, fix=1)
        elif state == 6:
            atoms.set_distance(i, Cnear, 1.4600, fix=1)
        x = atoms.get_positions()
        #Separate into a side
        for j in range(-8, -6):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-4, -2):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-6, -4):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        for j in range(-2, 0):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        atoms.set_positions(x)
        atoms.set_angle(Cnear, i, total_atoms-5, 131.58, mask=mask)
        atoms.set_distance(total_atoms-5, i, 1.367, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-5, 6.0, mask=mask)
        mask[-5] = 0
        atoms.set_angle(Cnear, i, total_atoms-8, 120.8, mask=mask)
        atoms.set_distance(total_atoms-8, i, 1.381, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-8, -175.0, mask=mask)
        mask[-8] = 0
        atoms.set_angle(i, total_atoms-5, total_atoms-6, 107.3, mask=mask)
        atoms.set_distance(total_atoms-6, total_atoms-5, 1.424, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-6, 179.5, mask=mask)
        mask[-6] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-7, 109.8, mask=mask)
        atoms.set_distance(total_atoms-7, total_atoms-8, 1.381, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-7, -179.5, mask=mask)
        mask[-7] = 0
        atoms.set_angle(i, total_atoms-7, total_atoms-4, 126, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-5, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-4, 0.5, mask=mask)
        mask[-4] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-1, 126, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-8, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-1, -0.5, mask=mask)
        mask[-1] = 0
        atoms.set_angle(total_atoms-5, total_atoms-6, total_atoms-3, 126, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-6, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-5, total_atoms-6,\
                           total_atoms-3, 179.5, mask=mask)
        mask[-3] = 0
        atoms.set_angle(total_atoms-8, total_atoms-7, total_atoms-2, 126, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-7, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-8, total_atoms-7,\
                           total_atoms-2, -179.5, mask=mask)
        mask[-2] = 0
    elif place == 2:
        atoms.append('C')
        atoms.append('N')
        for j in range(2):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear = near(i, atoms)
        total_atoms = len(atoms)
        mask = [0]*(total_atoms-8) + [1]*8
        elements = atoms.get_chemical_symbols()
        elements[i] = 'C'
        if shuffle == 0:
            CnearH = near(Cnear, atoms, 'H', ignore=i)
        elif shuffle == 1:
            CnearH = near(Cnear, atoms, 'C', ignore=i)
        atoms.set_chemical_symbols(elements)
        if state == 1:
            atoms.set_distance(i, Cnear, 1.5100, fix=1)
        elif state == 6:
            atoms.set_distance(i, Cnear, 1.4600, fix=1)
        x = atoms.get_positions()
        #Separate into a side
        for j in range(-8, -6):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-4, -2):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-6, -4):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        for j in range(-2, 0):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        atoms.set_positions(x)
        atoms.set_angle(Cnear, i, total_atoms-8, 125.4, mask=mask)
        atoms.set_distance(total_atoms-8, i, 1.367, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-8, 5.4, mask=mask)
        mask[-8] = 0
        atoms.set_angle(Cnear, i, total_atoms-5, 126.75, mask=mask)
        atoms.set_distance(total_atoms-5, i, 1.424, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-5, -175.0, mask=mask)
        mask[-5] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-7, 107.7, mask=mask)
        atoms.set_distance(total_atoms-7, total_atoms-8, 1.381, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-7, 179.5, mask=mask)
        mask[-7] = 0
        atoms.set_angle(i, total_atoms-5, total_atoms-6, 107.3, mask=mask)
        atoms.set_distance(total_atoms-6, total_atoms-5, 1.367, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-6, -179.5, mask=mask)
        mask[-6] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-4, 126, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-8, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-4, 0.5, mask=mask)
        mask[-4] = 0
        atoms.set_angle(i, total_atoms-5, total_atoms-1, 126, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-5, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-1, -0.5, mask=mask)
        mask[-1] = 0
        atoms.set_angle(total_atoms-8, total_atoms-7, total_atoms-3, 126, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-7, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-8, total_atoms-7,\
                           total_atoms-3, 179.5, mask=mask)
        mask[-3] = 0
        atoms.set_angle(total_atoms-5, total_atoms-6, total_atoms-2, 126, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-6, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-5, total_atoms-6,\
                           total_atoms-2, -179.5, mask=mask)
        mask[-2] = 0
    elif place == -2:
        atoms.append('C')
        atoms.append('N')
        for j in range(2):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear = near(i, atoms)
        total_atoms = len(atoms)
        mask = [0]*(total_atoms-8) + [1]*8
        elements = atoms.get_chemical_symbols()
        elements[i] = 'C'
        if shuffle == 0:
            CnearH = near(Cnear, atoms, 'H', ignore=i)
        elif shuffle == 1:
            CnearH = near(Cnear, atoms, 'C', ignore=i)
        atoms.set_chemical_symbols(elements)
        if state == 1:
            atoms.set_distance(i, Cnear, 1.5100, fix=1)
        elif state == 6:
            atoms.set_distance(i, Cnear, 1.4600, fix=1)
        x = atoms.get_positions()
        #Separate into a side
        for j in range(-8, -6):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-4, -2):
            x[j][0] = x[j][0] - 50
            x[j][1] = x[j][1] - 50
            x[j][2] = x[j][2] - 50
        for j in range(-6, -4):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        for j in range(-2, 0):
            x[j][0] = x[j][0] + 50
            x[j][1] = x[j][1] + 50
            x[j][2] = x[j][2] + 50
        atoms.set_positions(x)
        atoms.set_angle(Cnear, i, total_atoms-8, 125.4, mask=mask)
        atoms.set_distance(total_atoms-8, i, 1.367, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-8, 175.0, mask=mask)
        mask[-8] = 0
        atoms.set_angle(Cnear, i, total_atoms-5, 126.75, mask=mask)
        atoms.set_distance(total_atoms-5, i, 1.424, fix=1, mask=mask)
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-5, -5.4, mask=mask)
        mask[-5] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-7, 107.7, mask=mask)
        atoms.set_distance(total_atoms-7, total_atoms-8, 1.381, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-7, -179.5, mask=mask)
        mask[-7] = 0
        atoms.set_angle(i, total_atoms-5, total_atoms-6, 107.3, mask=mask)
        atoms.set_distance(total_atoms-6, total_atoms-5, 1.367, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-6, 179.5, mask=mask)
        mask[-6] = 0
        atoms.set_angle(i, total_atoms-5, total_atoms-4, 126, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-5, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-5, total_atoms-4, 0.5, mask=mask)
        mask[-4] = 0
        atoms.set_angle(i, total_atoms-8, total_atoms-1, 126, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-8, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms-8, total_atoms-1, -0.5, mask=mask)
        mask[-1] = 0
        atoms.set_angle(total_atoms-8, total_atoms-7, total_atoms-3, 126, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-7, 1.0800, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-8, total_atoms-7,\
                           total_atoms-3, 179.5, mask=mask)
        mask[-3] = 0
        atoms.set_angle(total_atoms-5, total_atoms-6, total_atoms-2, 126, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-6, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms-5, total_atoms-6,\
                           total_atoms-2, -179.5, mask=mask)
        mask[-2] = 0
    return atoms

def ring3(atoms, i, state=1, strict={1:'X'}):
    for j in range(2):
        atoms.append('C')
    for j in range(5):
        atoms.append('H')
    total_atoms = len(atoms)
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    mask = [0]*(total_atoms-7) + [1]*7
    elements = atoms.get_chemical_symbols()
    elements[i] = 'C'
    try:
        if strict[1] == 'N':
            elements[i] = 'N'
        else:
            elements[i] = 'C'
    except:
        elements[i] = 'C'
    atoms.set_chemical_symbols(elements)
    if elements[i] == 'C':
        if state == 1:
            atoms.set_distance(i, Cnear, 1.5290, fix=1)
        elif state == 2 or state == 82:
            atoms.set_distance(i, Cnear, 1.5095, fix=1)
        elif state == 3:
            atoms.set_distance(i, Cnear, 1.4700, fix=1)
        elif state == 6:
            atoms.set_distance(i, Cnear, 1.5105, fix=1)
        elif state == 7:
            atoms.set_distance(i, Cnear, 1.4480, fix=1)
        elif state == 8:
            atoms.set_distance(i, Cnear, 1.4100, fix=1)
        elif state == 16:
            atoms.set_distance(i, Cnear, 1.8100, fix=1)
        elif state == 162:
            atoms.set_distance(i, Cnear, 1.7900, fix=1)
        elif state == 164:
            atoms.set_distance(i, Cnear, 1.7700, fix=1)
    elif elements[i] == 'N':
        if state == 1:
            atoms.set_distance(i, Cnear, 1.448, fix=1)
        elif state == 2 or state == 6:
            atoms.set_distance(i, Cnear, 1.340, fix=1)
        else:
            atoms.set_distance(i, Cnear, 1.448, fix=1)
    x = atoms.get_positions()
    x[-5][0] = x[-5][0] - 50
    x[-5][1] = x[-5][1] - 50
    x[-5][2] = x[-5][2] - 50
    x[-7][0] = x[-7][0] - 50
    x[-7][1] = x[-7][1] - 50
    x[-7][2] = x[-7][2] - 50
    x[-6][0] = x[-6][0] + 50
    x[-6][1] = x[-6][1] + 50
    x[-6][2] = x[-6][2] + 50
    for j in range(-4, -2):
        x[j][0] = x[j][0] - 50
        x[j][1] = x[j][1] - 50
        x[j][2] = x[j][2] - 50
    for j in range(-2, 0):
        x[j][0] = x[j][0] + 50
        x[j][1] = x[j][1] + 50
        x[j][2] = x[j][2] + 50
    atoms.set_positions(x)
    if elements[i] == 'C':
        atoms.set_angle(Cnear, i, total_atoms-5, 109.5, mask=mask)
        atoms.set_distance(total_atoms-5, i, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-5, -150, mask=mask)
        mask[-5] = 0
    else:
        mask[-5] = 0
    if elements[i] == 'N':
        atoms.set_angle(Cnear, i, total_atoms-7, 122.5, mask=mask)
        atoms.set_distance(total_atoms-7, i, 1.4480, fix=1, mask=mask)
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-7, 0.5, mask=mask)
        mask[-7] = 0
        atoms.set_angle(Cnear, i, total_atoms-6, 122.5, mask=mask)
        atoms.set_distance(total_atoms-6, i, 1.4480, fix=1, mask=mask)
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-6, 73.5, mask=mask)
        mask[-6] = 0
    elif elements[i] == 'C':
        atoms.set_angle(Cnear, i, total_atoms-7, 122.5, mask=mask)
        atoms.set_distance(total_atoms-7, i, 1.5290, fix=1, mask=mask)
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-7, 0.5, mask=mask)
        mask[-7] = 0
        atoms.set_angle(Cnear, i, total_atoms-6, 122.5, mask=mask)
        atoms.set_distance(total_atoms-6, i, 1.5290, fix=1, mask=mask)
        atoms.set_dihedral(CnearC, Cnear, i, total_atoms-6, 73.5, mask=mask)
        mask[-6] = 0
    atoms.set_angle(i, total_atoms-7, total_atoms-4, 122.5, mask=mask)
    atoms.set_distance(total_atoms-4, total_atoms-7, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-7, total_atoms-4, 5, mask=mask)
    mask[-4] = 0
    atoms.set_angle(i, total_atoms-7, total_atoms-3, 122.5, mask=mask)
    atoms.set_distance(total_atoms-3, total_atoms-7, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-7, total_atoms-3, -140, mask=mask)
    mask[-3] = 0
    atoms.set_angle(i, total_atoms-6, total_atoms-2, 122.5, mask=mask)
    atoms.set_distance(total_atoms-2, total_atoms-6, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-6, total_atoms-2, -5, mask=mask)
    mask[-2] = 0
    atoms.set_angle(i, total_atoms-6, total_atoms-1, 122.5, mask=mask)
    atoms.set_distance(total_atoms-1, total_atoms-6, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-6, total_atoms-1, 140, mask=mask)
    mask[-1] = 0
    for j in range(2):
        try:
            ring_element = strict[j+2]
        except: 
            continue
        elements = atoms.get_chemical_symbols()
        if ring_element == 'O':
            elements[-7+j] = 'O'
            atoms.set_chemical_symbols(elements)
            mask = [0] * total_atoms
            mask[-7+j] = 1
            if j == 0:
                for l in range(10):
                    atoms.set_distance(total_atoms-6, total_atoms-7, 1.4100, mask=mask)
                    atoms.set_distance(i, total_atoms-7, 1.4100, mask=mask)
            elif j == 1:
                for l in range(10):
                    atoms.set_distance(i, total_atoms-6, 1.4100, mask=mask)
                    atoms.set_distance(total_atoms-7, total_atoms-6, 1.4100, mask=mask)
        elif ring_element == 'S':
            elements[-7+j] = 'S'
            atoms.set_chemical_symbols(elements)
            mask = [0] * total_atoms
            mask[-7+j] = 1
            if j == 0:
                for l in range(10):
                    atoms.set_distance(total_atoms-6, total_atoms-7, 1.8100, mask=mask)
                    atoms.set_distance(i, total_atoms-7, 1.8100, mask=mask)
            elif j == 1:
                for l in range(10):
                    atoms.set_distance(i, total_atoms-6, 1.8100, mask=mask)
                    atoms.set_distance(total_atoms-7, total_atoms-6, 1.8100, mask=mask)
        elif ring_element == 'N':
            elements[-7+j] = 'N'
            atoms.set_chemical_symbols(elements)
            mask = [0] * total_atoms
            mask[-7+j] = 1
            if j == 0:
                for l in range(10):
                    atoms.set_distance(total_atoms-6, total_atoms-7, 1.4480, mask=mask)
                    atoms.set_distance(i, total_atoms-7, 1.4480, mask=mask)
                mask[-7] = 0
                mask[-3] = 1
                for l in range(10):
                    atoms.set_angle(i, total_atoms-7, total_atoms-3, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-6, total_atoms-7, total_atoms-3, 109.5, mask=mask)
                atoms.set_distance(total_atoms-3, total_atoms-7, 1.0100, fix=1, mask=mask)
                mask[-3] = 0
            elif j == 1:
                for l in range(10):
                    atoms.set_distance(i, total_atoms-6, 1.4480, mask=mask)
                    atoms.set_distance(total_atoms-7, total_atoms-6, 1.4480, mask=mask)
                mask[-6] = 0
                mask[-1] = 1
                for l in range(10):
                    atoms.set_angle(i, total_atoms-6, total_atoms-1, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-7, total_atoms-6, total_atoms-1, 109.5, mask=mask)
                atoms.set_distance(total_atoms-1, total_atoms-6, 1.0100, fix=1, mask=mask)
                mask[-1] = 0
        mask[-7+j] = 0
    try:
        if strict[1] == 'N':
            del atoms[-5]
        else:
            pass
    except:
        pass
    deleted=[]
    for j in range(2):
        try:
            ring_element = strict[j+2]
        except: 
            continue
        if ring_element == 'N':
            if j + 2 == 2:
                deleted.append(total_atoms-4)
            if j + 2 == 3:
                deleted.append(total_atoms-2)
        else:
            if j + 2 == 2:
                deleted.append(total_atoms-3)
                deleted.append(total_atoms-4)
            if j + 2 == 3:
                deleted.append(total_atoms-1)
                deleted.append(total_atoms-2)
    deleted.sort(reverse=True)
    for d in deleted:
        del atoms[d]
    return atoms

def ring5(atoms, i, state=1, strict={1:'X'}, shuffle=0):
    Cnear = near(i, atoms)
    try:
        if strict[1] == 'N':
            pass
    except:
        strict[1] = 'X'
    if strict[1] == 'N':
        atoms = NH2(atoms, i, state)
    else:
        atoms = CH3(atoms, i, state)
    for j in range(2):
        atoms.append('C')
    for j in range(8):
        atoms.append('H')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-10) + [1]*10
    elements = atoms.get_chemical_symbols()
    if shuffle == 0:
        box=[-12, -11]
    elif shuffle == 1:
        box=[-13, -12]
    elif shuffle == 2:
        box=[-11, -13]
    for j in box:
        elements[total_atoms+j] = 'C'
        atoms.set_chemical_symbols(elements)
        if strict[1] == 'N':
            atoms.set_distance(total_atoms+j, i, 1.4480, fix=1)
        else:
            atoms.set_distance(total_atoms+j, i, 1.5290, fix=1)
    x = atoms.get_positions()
    x[-10][0] = x[-10][0] - 50
    x[-10][1] = x[-10][1] - 50
    x[-10][2] = x[-10][2] - 50
    x[-9][0] = x[-9][0] + 50
    x[-9][1] = x[-9][1] + 50
    x[-9][2] = x[-9][2] + 50
    for j in range(-8, -4):
        x[j][0] = x[j][0] - 50
        x[j][1] = x[j][1] - 50
        x[j][2] = x[j][2] - 50
    for j in range(-4, 0):
        x[j][0] = x[j][0] + 50
        x[j][1] = x[j][1] + 50
        x[j][2] = x[j][2] + 50
    if strict[1] == 'N':
        atoms.set_positions(x)
        atoms.set_angle(i, total_atoms+box[0], total_atoms-10, 104.98, mask=mask)
        atoms.set_distance(total_atoms-10, total_atoms+box[0], 1.5290, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[0], total_atoms-10,\
                           164.01, mask=mask)
        mask[-10] = 0
        atoms.set_angle(i, total_atoms+box[1], total_atoms-9, 104.98, mask=mask)
        atoms.set_distance(total_atoms-9, total_atoms+box[1], 1.5290, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[1], total_atoms-9,\
                           -164.01, mask=mask)
        mask[-9] = 0
    else:
        atoms.set_positions(x)
        atoms.set_angle(i, total_atoms+box[0], total_atoms-10, 99.37, mask=mask)
        atoms.set_distance(total_atoms-10, total_atoms+box[0], 1.5290, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[0], total_atoms-10,\
                           -161, mask=mask)
        mask[-10] = 0
        atoms.set_angle(i, total_atoms+box[1], total_atoms-9, 99.37, mask=mask)
        atoms.set_distance(total_atoms-9, total_atoms+box[1], 1.5290,\
                           fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[1], total_atoms-9,\
                           161, mask=mask)
        mask[-9] = 0
    atoms.set_angle(i, total_atoms+box[0], total_atoms-8, 109.5, mask=mask)
    atoms.set_distance(total_atoms-8, total_atoms+box[0], 1.0900,\
                       fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[0], total_atoms-8,\
                       75, mask=mask)
    mask[-8] = 0
    atoms.set_angle(i, total_atoms+box[0], total_atoms-7, 109.5, mask=mask)
    atoms.set_distance(total_atoms-7, total_atoms+box[0], 1.0900,\
                       fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[0], total_atoms-7,\
                       -45, mask=mask)
    mask[-7] = 0
    atoms.set_angle(i, total_atoms+box[1], total_atoms-6, 109.5, mask=mask)
    atoms.set_distance(total_atoms-6, total_atoms+box[1], 1.0900,\
                       fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[1], total_atoms-6,\
                       -75, mask=mask)
    mask[-6] = 0
    atoms.set_angle(i, total_atoms+box[1], total_atoms-5, 109.5, mask=mask)
    atoms.set_distance(total_atoms-5, total_atoms+box[1], 1.0900,\
                       fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[1], total_atoms-5,\
                       45, mask=mask)
    mask[-5] = 0
    if strict[1] == 'N':
        atoms.set_angle(total_atoms+box[0], total_atoms-10,\
                        total_atoms-4, 109.5, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-10, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0], total_atoms-10,\
                           total_atoms-4, -140.5, mask=mask)
        mask[-4] = 0
        atoms.set_angle(total_atoms+box[0], total_atoms-10,\
                        total_atoms-3, 109.5, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-10, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0], total_atoms-10,\
                           total_atoms-3, 95, mask=mask)
        mask[-3] = 0
        atoms.set_angle(total_atoms+box[1], total_atoms-9,\
                        total_atoms-2, 109.5, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-9, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1], total_atoms-9,\
                           total_atoms-2, 140.5, mask=mask)
        mask[-2] = 0
        atoms.set_angle(total_atoms+box[1], total_atoms-9,\
                        total_atoms-1, 109.5, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-9, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1], total_atoms-9,\
                           total_atoms-1, -95, mask=mask)
        mask[-1] = 0
    else:
        atoms.set_angle(total_atoms+box[0], total_atoms-10,\
                        total_atoms-4, 109.5, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-10, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0], total_atoms-10,\
                           total_atoms-4, -90.5, mask=mask)
        mask[-4] = 0
        atoms.set_angle(total_atoms+box[0], total_atoms-10,\
                        total_atoms-3, 109.5, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-10, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0], total_atoms-10,\
                           total_atoms-3, 145, mask=mask)
        mask[-3] = 0
        atoms.set_angle(total_atoms+box[1], total_atoms-9,\
                        total_atoms-2, 109.5, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-9, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1], total_atoms-9,\
                           total_atoms-2, 90.5, mask=mask)
        mask[-2] = 0
        atoms.set_angle(total_atoms+box[1], total_atoms-9,\
                        total_atoms-1, 109.5, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-9, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1], total_atoms-9,\
                           total_atoms-1, -145, mask=mask)
        mask[-1] = 0
    l=0
    for j in range(4):
        try:
            ring_element = strict[j+2]
        except: 
            continue
        if ring_element == 'O':
            mask.append(1)
            atoms.append(ring_element)
            if j + 2 == 2:
                atoms.set_distance(total_atoms+l, total_atoms+box[0], 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-10, 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, i, 1.4100, fix=1, mask=mask)
            elif j + 2 == 3:
                atoms.set_distance(total_atoms+l, total_atoms-10, 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[0], 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-9, 1.4100, fix=1, mask=mask)
            elif j + 2 == 4:
                atoms.set_distance(total_atoms+l, total_atoms-9, 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[1], 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-10, 1.4100, fix=1, mask=mask)
            elif j + 2 == 5:
                atoms.set_distance(total_atoms+l, total_atoms+box[1], 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, i, 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-9, 1.4100, fix=1, mask=mask)
            mask[-1] = 0
            l = l + 1
        elif ring_element == 'S':
            mask.append(1)
            atoms.append(ring_element)
            if j + 2 == 4:
                atoms.set_distance(total_atoms+l, total_atoms-9, 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[1], 1.8100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-10, 1.8100, fix=1, mask=mask)
            elif j + 2 == 3:
                atoms.set_distance(total_atoms+l, total_atoms-10, 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[0], 1.8100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-9, 1.8100, fix=1, mask=mask)
            elif j + 2 == 5:
                atoms.set_distance(total_atoms+l, total_atoms+box[1], 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, i, 1.8100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-9, 1.8100, fix=1, mask=mask)
            elif j + 2 == 2:
                atoms.set_distance(total_atoms+l, total_atoms+box[0], 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-10, 1.8100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, i, 1.8100, fix=1, mask=mask)
            mask[-1] = 0
            l = l + 1
        elif ring_element == 'N':
            mask.append(1)
            mask.append(1)
            atoms.append(ring_element)
            atoms.append('H')
            x = atoms.get_positions()
            x[-1][0] = x[-1][0] - 50
            x[-1][1] = x[-1][1] - 50
            x[-1][2] = x[-1][2] - 50
            atoms.set_positions(x)
            if j + 2 == 4:
                atoms.set_distance(total_atoms+l, total_atoms-9, 0, fix=1, mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[1], 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-10, 1.4480, fix=1, mask=mask)
                mask[-2] = 0
                for k in range(10):
                    atoms.set_angle(total_atoms-10, total_atoms+l, \
                                    total_atoms+l+1, 109.5, mask=mask)
                    atoms.set_angle(total_atoms+box[1], total_atoms+l,\
                                    total_atoms+l+1, 109.5, mask=mask)
                atoms.set_distance(total_atoms+l+1, total_atoms+l, 1.0100, fix=1, mask=mask)
            elif j + 2 == 3:
                atoms.set_distance(total_atoms+l, total_atoms-10, 0, fix=1, mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[0], 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-9, 1.4480, fix=1, mask=mask)
                mask[-2] = 0
                for k in range(10):
                    atoms.set_angle(total_atoms+box[0], total_atoms+l,\
                                    total_atoms+l+1, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-9, total_atoms+l,\
                                    total_atoms+l+1, 109.5, mask=mask)
                atoms.set_distance(total_atoms+l+1, total_atoms+l, 1.0100, fix=1, mask=mask)
            elif j + 2 == 5:
                atoms.set_distance(total_atoms+l, total_atoms+box[1], 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, i, 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-9, 1.4480, fix=1, mask=mask)
                mask[-2] = 0
                for k in range(10):
                    atoms.set_angle(i, total_atoms+l, total_atoms+l+1, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-9, total_atoms+l,\
                                    total_atoms+l+1, 109.5, mask=mask)
                atoms.set_distance(total_atoms+l+1, total_atoms+l, 1.0100, fix=1, mask=mask)
            elif j + 2 == 2:
                atoms.set_distance(total_atoms+l, total_atoms+box[0], 0, fix=1, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-10, 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, i, 1.4480, fix=1, mask=mask)
                mask[-2] = 0
                for k in range(10):
                    atoms.set_angle(i, total_atoms+l, total_atoms+l+1, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-10, total_atoms+l,\
                                    total_atoms+l+1, 109.5, mask=mask)
                atoms.set_distance(total_atoms+l+1, total_atoms+l, 1.0100, fix=1, mask=mask)
            mask[-1] = 0
            l = l + 2
    deleted = []
    for j in range(5):
        try:
            ring_element = strict[j+2]
        except: 
            continue
        if j + 2 == 4:
            deleted.append(total_atoms-1)
            deleted.append(total_atoms-2)
            deleted.append(total_atoms-9)
        if j + 2 == 3:
            deleted.append(total_atoms-3)
            deleted.append(total_atoms-4)
            deleted.append(total_atoms-10)
        if j + 2 == 5:
            deleted.append(total_atoms-5)
            deleted.append(total_atoms-6)
            deleted.append(total_atoms+box[1])
        if j + 2 == 2:
            deleted.append(total_atoms-7)
            deleted.append(total_atoms-8)
            deleted.append(total_atoms+box[0])
    deleted.sort(reverse=True)
    for d in deleted:
        del atoms[d]
    return atoms

def ring6(atoms, i, state=1, strict={1:'X'}, shuffle=0):
    Cnear = near(i, atoms)
    try:
        if strict[1] == 'N':
            pass
    except:
        strict[1] = 'X'
    if strict[1] == 'N':
        atoms = NH2(atoms, i, state)
    else:
        atoms = CH3(atoms, i, state)
    for j in range(3):
        atoms.append('C')
    for j in range(10):
        atoms.append('H')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-13) + [1]*13
    elements = atoms.get_chemical_symbols()
    if shuffle == 0:
        box=[-15, -14]
    elif shuffle == 1:
        box=[-16, -15]
    elif shuffle == 2:
        box=[-14, -16]
    for j in box:
        elements[total_atoms+j] = 'C'
        atoms.set_chemical_symbols(elements)
        if strict[1] == 'N':
            atoms.set_distance(total_atoms+j, i, 1.4480, fix=1)
        else:
            atoms.set_distance(total_atoms+j, i, 1.5290, fix=1)
    x = atoms.get_positions()
    x[-13][0] = x[-13][0] - 50
    x[-13][1] = x[-13][1] - 50
    x[-13][2] = x[-13][2] - 50
    x[-12][0] = x[-12][0] - 50
    x[-12][1] = x[-12][1] - 50
    x[-12][2] = x[-12][2] - 50
    x[-11][0] = x[-11][0] + 50
    x[-11][1] = x[-11][1] + 50
    x[-11][2] = x[-11][2] + 50
    for j in range(-10, -5):
        x[j][0] = x[j][0] - 50
        x[j][1] = x[j][1] - 50
        x[j][2] = x[j][2] - 50
    for j in range(-5, 0):
        x[j][0] = x[j][0] + 50
        x[j][1] = x[j][1] + 50
        x[j][2] = x[j][2] + 50
    atoms.set_positions(x)
    if strict[1] == 'N':
        atoms.set_angle(i, total_atoms+box[0], total_atoms-13, 115.7, mask=mask)
        atoms.set_distance(total_atoms-13, total_atoms+box[0],\
                           1.5290, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[0],\
                           total_atoms-13, 170, mask=mask)
        mask[-13] = 0
        atoms.set_angle(i, total_atoms+box[1], total_atoms-12, 115.7, mask=mask)
        atoms.set_distance(total_atoms-12, total_atoms+box[1],\
                           1.5290, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[1],\
                           total_atoms-12, -170, mask=mask)
        mask[-12] = 0
        atoms.set_angle(total_atoms+box[0], total_atoms-13,\
                        total_atoms-11, 113.7, mask=mask)
        for k in range(10):
            atoms.set_distance(total_atoms-11, total_atoms-12, 1.5290, fix=1, mask=mask)
            atoms.set_distance(total_atoms-11, total_atoms-13, 1.5290, fix=1, mask=mask)
        mask[-11] = 0
    else:
        atoms.set_angle(i, total_atoms+box[0], total_atoms-13, 107.4, mask=mask)
        atoms.set_distance(total_atoms-13, total_atoms+box[0],\
                           1.5290, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[0],\
                           total_atoms-13, 170, mask=mask)
        mask[-13] = 0
        atoms.set_angle(i, total_atoms+box[1], total_atoms-12, 107.4, mask=mask)
        atoms.set_distance(total_atoms-12, total_atoms+box[1],\
                           1.5290, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[1],\
                           total_atoms-12, -170, mask=mask)
        mask[-12] = 0
        atoms.set_angle(total_atoms+box[0], total_atoms-13,\
                        total_atoms-11, 107.4, mask=mask)
        atoms.set_distance(total_atoms-11, total_atoms-13, 1.5290, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0], total_atoms-13,\
                           total_atoms-11, 55, mask=mask)
        mask[-11] = 0
    atoms.set_angle(i, total_atoms+box[0], total_atoms-10, 109.5, mask=mask)
    atoms.set_distance(total_atoms-10, total_atoms+box[0],\
                       1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[0],\
                       total_atoms-10, -60, mask=mask)
    mask[-10] = 0
    atoms.set_angle(i, total_atoms+box[0], total_atoms-9, 109.5, mask=mask)
    atoms.set_distance(total_atoms-9, total_atoms+box[0], 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[0], total_atoms-9, 60, mask=mask)
    mask[-9] = 0
    atoms.set_angle(i, total_atoms+box[1], total_atoms-8, 109.5, mask=mask)
    atoms.set_distance(total_atoms-8, total_atoms+box[1], 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[1], total_atoms-8, -60, mask=mask)
    mask[-8] = 0
    atoms.set_angle(i, total_atoms+box[1], total_atoms-7, 109.5, mask=mask)
    atoms.set_distance(total_atoms-7, total_atoms+box[1], 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[1], total_atoms-7, 60, mask=mask)
    mask[-7] = 0
    if strict[1] == 'N':
        atoms.set_angle(total_atoms+box[0], total_atoms-13,\
                        total_atoms-6, 109.5, mask=mask)
        atoms.set_distance(total_atoms-6, total_atoms-13, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0], total_atoms-13,\
                           total_atoms-6, 165, mask=mask)
        mask[-6] = 0
        atoms.set_angle(total_atoms+box[0], total_atoms-13,\
                        total_atoms-5, 109.5, mask=mask)
        atoms.set_distance(total_atoms-5, total_atoms-13, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0], total_atoms-13,\
                           total_atoms-5, 59.5, mask=mask)
        mask[-5] = 0
        atoms.set_angle(total_atoms+box[1], total_atoms-12,\
                        total_atoms-4, 109.5, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-12, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1], total_atoms-12,\
                           total_atoms-4, -165, mask=mask)
        mask[-4] = 0
        atoms.set_angle(total_atoms+box[1], total_atoms-12,\
                        total_atoms-3, 109.5, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-12, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1], total_atoms-12,\
                           total_atoms-3, -59.5, mask=mask)
        mask[-3] = 0
    else:
        atoms.set_angle(total_atoms+box[0], total_atoms-13,\
                        total_atoms-6, 109.5, mask=mask)
        atoms.set_distance(total_atoms-6, total_atoms-13, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0], total_atoms-13,\
                           total_atoms-6, -75, mask=mask)
        mask[-6] = 0
        atoms.set_angle(total_atoms+box[0], total_atoms-13,\
                        total_atoms-5, 109.5, mask=mask)
        atoms.set_distance(total_atoms-5, total_atoms-13, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0], total_atoms-13,\
                           total_atoms-5, 179.5, mask=mask)
        mask[-5] = 0
        atoms.set_angle(total_atoms+box[1], total_atoms-12,\
                        total_atoms-4, 109.5, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-12, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1], total_atoms-12,\
                           total_atoms-4, 75, mask=mask)
        mask[-4] = 0
        atoms.set_angle(total_atoms+box[1], total_atoms-12,\
                        total_atoms-3, 109.5, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-12, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1], total_atoms-12,\
                           total_atoms-3, -179.5, mask=mask)
        mask[-3] = 0
    atoms.set_angle(total_atoms-13, total_atoms-11,\
                    total_atoms-2, 109.5, mask=mask)
    atoms.set_distance(total_atoms-2, total_atoms-11, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms+box[0], total_atoms-13,\
                       total_atoms-11, total_atoms-2, -179.5, mask=mask)
    mask[-2] = 0
    atoms.set_angle(total_atoms-13, total_atoms-11,\
                    total_atoms-1, 109.5, mask=mask)
    atoms.set_distance(total_atoms-1, total_atoms-11, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms+box[0], total_atoms-12,\
                       total_atoms-11, total_atoms-1, 300, mask=mask)
    mask[-1] = 0
    l = 0
    for j in range(5):
        try:
            ring_element = strict[j+2]
        except: 
            continue
        if ring_element == 'O':
            mask.append(1)
            atoms.append(ring_element)
            if j + 2 == 4:
                atoms.set_distance(total_atoms+l, total_atoms-11, 0.1, fix=1, mask=mask)
                atoms.set_angle(total_atoms+box[0], total_atoms-13,\
                                total_atoms+l, 100, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms-13, 1.4100, fix=1, mask=mask)
                atoms.set_dihedral(i, total_atoms+box[0], total_atoms-13,\
                                   total_atoms+l, 47.9, mask=mask)
            elif j + 2 == 3:
                atoms.set_distance(total_atoms+l, total_atoms-13, 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms+box[1], 5, mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[0], 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-11, 1.4100, fix=1, mask=mask)
            elif j + 2 == 5:
                atoms.set_distance(total_atoms+l, total_atoms-12, 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms+box[0], 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[1], 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-11, 1.4100, fix=1, mask=mask)
            elif j + 2 == 2:
                atoms.set_distance(total_atoms+l, total_atoms+box[0], 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms-12, 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-13, 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, i, 1.4100, fix=1, mask=mask)
            elif j + 2 == 6:
                atoms.set_distance(total_atoms+l, total_atoms+box[1], 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms-13, 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-12, 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, i, 1.4100, fix=1, mask=mask)
            l = l + 1
        elif ring_element == 'S':
            mask.append(1)
            atoms.append(ring_element)
            if j + 2 == 4:
                atoms.set_distance(total_atoms+l, total_atoms-11, 0.1, fix=1, mask=mask)
                atoms.set_angle(total_atoms+box[0], total_atoms-13,\
                                total_atoms+l, 127.35, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms-13, 1.8100, fix=1, mask=mask)
                atoms.set_dihedral(i, total_atoms+box[0], total_atoms-13,\
                                   total_atoms+l, 47.9, mask=mask)
            elif j + 2 == 3:
                atoms.set_distance(total_atoms+l, total_atoms-13, 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms+box[1], 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[0], 1.8100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-11, 1.8100, fix=1, mask=mask)
            elif j + 2 == 5:
                atoms.set_distance(total_atoms+l, total_atoms-12, 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms+box[0], 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[1], 1.8100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-11, 1.8100, fix=1, mask=mask)
            elif j + 2 == 2:
                atoms.set_distance(total_atoms+l, total_atoms+box[0], 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms-12, 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-13, 1.8100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, i, 1.8100, fix=1, mask=mask)
            elif j + 2 == 6:
                atoms.set_distance(total_atoms+l, total_atoms+box[1], 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms-13, 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-12, 1.8100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, i, 1.8100, fix=1, mask=mask)
            l = l + 1
        elif ring_element == 'N':
            mask.append(1)
            mask.append(1)
            atoms.append(ring_element)
            atoms.append('H')
            x = atoms.get_positions()
            x[-1][0] = x[-1][0] - 50
            x[-1][1] = x[-1][1] - 50
            x[-1][2] = x[-1][2] - 50
            atoms.set_positions(x)
            if j + 2 == 4:
                atoms.set_distance(total_atoms+l, total_atoms-11, 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, i, 5, mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-13, 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-12, 1.4480, fix=1, mask=mask)
                mask[-2] = 0
                atoms.set_angle(total_atoms-13, total_atoms+l, total_atoms+l+1, 109.5, mask=mask)
                atoms.set_distance(total_atoms+l+1, total_atoms+l, 1.0100, fix=1, mask=mask)
                atoms.set_dihedral(total_atoms+l+1, total_atoms+l,\
                                   total_atoms-13, total_atoms+box[0], -150, mask=mask)
            elif j + 2 == 3:
                atoms.set_distance(total_atoms+l, total_atoms-13, 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms-14, 5, mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[0], 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-11, 1.4480, fix=1, mask=mask)
                mask[-2] = 0
                for k in range(10):
                    atoms.set_angle(total_atoms+box[0], total_atoms+l,\
                                    total_atoms+l+1, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-11, total_atoms+l, total_atoms+l+1, 109.5, mask=mask)
                atoms.set_distance(total_atoms+l+1, total_atoms+l, 1.0100, fix=1, mask=mask)
            elif j + 2 == 5:
                atoms.set_distance(total_atoms+l, total_atoms-12, 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms+box[0], 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms+box[1], 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, total_atoms-11, 1.4480, fix=1, mask=mask)
                mask[-2] = 0
                for k in range(10):
                    atoms.set_angle(total_atoms+box[0], total_atoms+l,\
                                    total_atoms+l+1, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-11, total_atoms+l, total_atoms+l+1, 109.5, mask=mask)
                atoms.set_distance(total_atoms+l+1, total_atoms+l, 1.0100, fix=1, mask=mask)
            elif j + 2 == 2:
                atoms.set_distance(total_atoms+l, total_atoms+box[0], 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms-12, 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-13, 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, i, 1.4480, fix=1, mask=mask)
                mask[-2] = 0
                for k in range(10):
                    atoms.set_angle(i, total_atoms+l, total_atoms+l+1, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-13, total_atoms+l,\
                                    total_atoms+l+1, 109.5, mask=mask)
                atoms.set_distance(total_atoms+l+1, total_atoms+l, 1.0100, fix=1, mask=mask)
            elif j + 2 == 6:
                atoms.set_distance(total_atoms+l, total_atoms+box[1], 0, fix=1, mask=mask)
                atoms.set_distance(total_atoms+l, total_atoms-13, 5, mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l, total_atoms-12, 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms+l, i, 1.4480, fix=1, mask=mask)
                mask[-2] = 0
                for k in range(10):
                    atoms.set_angle(i, total_atoms+l, total_atoms+l+1, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-12, total_atoms+l,\
                                    total_atoms+l+1, 109.5, mask=mask)
                atoms.set_distance(total_atoms+l+1, total_atoms+l, 1.0100, fix=1, mask=mask)
            l = l + 2
        mask[-1] = 0
    deleted=[]
    for j in range(5):
        try:
            ring_element = strict[j+2]
        except: 
            continue
        if j + 2 == 4:
            deleted.append(total_atoms-1)
            deleted.append(total_atoms-2)
            deleted.append(total_atoms-11)
        if j + 2 == 3:
            deleted.append(total_atoms-5)
            deleted.append(total_atoms-6)
            deleted.append(total_atoms-13)
        if j + 2 == 5:
            deleted.append(total_atoms-3)
            deleted.append(total_atoms-4)
            deleted.append(total_atoms-12)
        if j + 2 == 2:
            deleted.append(total_atoms-9)
            deleted.append(total_atoms-10)
            deleted.append(total_atoms+box[0])
        if j + 2 == 6:
            deleted.append(total_atoms-7)
            deleted.append(total_atoms-8)
            deleted.append(total_atoms+box[1])
    deleted.sort(reverse=True)
    for d in deleted:
        del atoms[d]
    return atoms
    
def SOMe(atoms, i, state=1):
    atoms.append('O')
    atoms.append('H')
    total_atoms = len(atoms)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'S'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    mask = [0]*(total_atoms-2) + [1]*2
    x = atoms.get_positions()
    x[-1][0] = x[-1][0] - 50
    x[-1][1] = x[-1][1] - 50
    x[-1][2] = x[-1][2] - 50
    if state == 1:
        atoms.set_distance(i, Cnear, 1.7900, fix=1)
    elif state == 2 or state == 6:
        atoms.set_distance(i, Cnear, 1.7700, fix=1)
    atoms.set_distance(total_atoms-2, i, 1.5300, fix=1, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-2, 107.0, mask=mask)
    atoms.set_dihedral(CnearC, Cnear, i, total_atoms-2, 65, mask=mask)
    mask[-2] = 0
    atoms.set_angle(total_atoms-2, i, total_atoms-1, 107.0, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-1, 96.0, mask=mask)
    atoms.set_distance(total_atoms-1, i, 1.7900, fix=1, mask=mask)
    atoms.set_dihedral(CnearC, Cnear, i, total_atoms-1, 175, mask=mask)
    atoms = CH3(atoms, total_atoms-1, state=162)
    return atoms

def SO2_(addition, atoms, i, state=1):
    atoms.append('O')
    atoms.append('O')
    atoms.append('H')
    total_atoms = len(atoms)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'S'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    CnearC = near(Cnear, atoms, 'C')
    mask = [0]*(total_atoms-3) + [1]*3
    x = atoms.get_positions()
    x[-3][0] = x[-3][0] - 50
    x[-3][1] = x[-3][1] - 50
    x[-3][2] = x[-3][2] - 50
    x[-2][0] = x[-2][0] + 50
    x[-2][1] = x[-2][1] + 50
    x[-2][2] = x[-2][2] + 50
    x[-1][0] = x[-1][0] - 50
    x[-1][1] = x[-1][1] - 50
    x[-1][2] = x[-1][2] - 50
    if state == 1:
        atoms.set_distance(i, Cnear, 1.7700, fix=1)
    elif state == 2 or state == 6:
        atoms.set_distance(i, Cnear, 1.7600, fix=1)
    atoms.set_distance(total_atoms-3, i, 1.5300, fix=1, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-3, 108.9, mask=mask)
    atoms.set_dihedral(CnearC, Cnear, i, total_atoms-2, 65, mask=mask)
    mask[-3] = 0
    atoms.set_distance(total_atoms-2, i, 1.5300, fix=1, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-2, 108.9, mask=mask)
    atoms.set_dihedral(CnearC, Cnear, i, total_atoms-2, -70, mask=mask)
    mask[-2] = 0
    atoms.set_angle(total_atoms-2, i, total_atoms-1, 108.9, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-1, 97.67, mask=mask)
    atoms.set_distance(total_atoms-1, i, 1.7700, fix=1, mask=mask)
    atoms.set_dihedral(CnearC, Cnear, i, total_atoms-1, 175, mask=mask)
    if addition == 'CH3' or addition == 'Me':
        atoms = CH3(atoms, total_atoms-1, state=164)
    if addition == 'OH':
        atoms = OH(atoms, total_atoms-1, state=164)
    #print(atoms.get_positions())
    return atoms

def D(member, atoms, i, state=1, shuffle=0):
    Cnear = near(i, atoms)
    atoms = CH3(atoms, i, state)
    for j in range(3):
        atoms.append('C')
    for j in range(8):
        atoms.append('H')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-11) + [1]*11
    elements = atoms.get_chemical_symbols()
    if shuffle == 0:
        box = [-13, -12]
    elif shuffle == 1:
        box = [-14, -13]
    elif shuffle == 2:
        box = [-12, -14]
    for j in box:
        mask[j] = 1
        elements[total_atoms+j] = 'C'
        atoms.set_chemical_symbols(elements)
        atoms.set_distance(total_atoms+j, i, 1.5290, fix=1)
        atoms.set_angle(Cnear, i, total_atoms+j, 113, mask=mask)
        mask[j] = 0
    x = atoms.get_positions()
    x[-11][0] = x[-11][0] - 50
    x[-11][1] = x[-11][1] - 50
    x[-11][2] = x[-11][2] - 50
    x[-10][0] = x[-10][0] + 50
    x[-10][1] = x[-10][1] + 50
    x[-10][2] = x[-10][2] + 50
    x[-9][0] = x[-9][0] - 50
    x[-9][1] = x[-9][1] - 50
    x[-9][2] = x[-9][2] - 50
    for j in range(-8, -4):
        x[j][0] = x[j][0] - 50
        x[j][1] = x[j][1] - 50
        x[j][2] = x[j][2] - 50
    for j in range(-4, 0):
        x[j][0] = x[j][0] + 50
        x[j][1] = x[j][1] + 50
        x[j][2] = x[j][2] + 50
    atoms.set_positions(x)
    atoms.set_angle(i, total_atoms+box[0], total_atoms-11, 102.3, mask=mask)
    atoms.set_distance(total_atoms-11, total_atoms+box[0],\
                       1.5290, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[0],\
                       total_atoms-11, 235, mask=mask)
    mask[-11] = 0
    atoms.set_angle(i, total_atoms+box[1], total_atoms-10, 115.5, mask=mask)
    atoms.set_distance(total_atoms-10, total_atoms+box[1],\
                       1.5290, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[1],\
                       total_atoms-10, 165, mask=mask)
    mask[-10] = 0
    atoms.set_angle(total_atoms+box[0], total_atoms-11,\
                    total_atoms-9, 110.5, mask=mask)
    atoms.set_dihedral(i, total_atoms+box[0], total_atoms-11,\
                       total_atoms-9, 290, mask=mask)
    for k in range(10):
        atoms.set_distance(total_atoms-9, total_atoms-11, 1.5290, fix=1, mask=mask)
        atoms.set_distance(total_atoms-9, total_atoms-10, 1.5290, fix=1, mask=mask)
    mask[-9] = 0
    atoms.set_angle(i, total_atoms+box[0], total_atoms-8, 109.5, mask=mask)
    atoms.set_distance(total_atoms-8, total_atoms+box[0],\
                       1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[0],\
                       total_atoms-8, -60, mask=mask)
    mask[-8] = 0
    atoms.set_angle(i, total_atoms+box[0], total_atoms-7, 109.5, mask=mask)
    atoms.set_distance(total_atoms-7, total_atoms+box[0],\
                       1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[0], total_atoms-7, 60, mask=mask)
    mask[-7] = 0
    atoms.set_angle(i, total_atoms+box[1], total_atoms-6, 109.5, mask=mask)
    atoms.set_distance(total_atoms-6, total_atoms+box[1],\
                       1.0900, fix=1, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms+box[1],\
                       total_atoms-6, -60, mask=mask)
    mask[-6] = 0
    atoms.set_angle(total_atoms+box[0], total_atoms-11,\
                    total_atoms-5, 109.5, mask=mask)
    atoms.set_distance(total_atoms-5, total_atoms-11, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(i, total_atoms+box[0], total_atoms-11,\
                       total_atoms-5, 179.5, mask=mask)
    mask[-5] = 0
    atoms.set_angle(total_atoms+box[1], total_atoms-10,\
                    total_atoms-4, 109.5, mask=mask)
    atoms.set_distance(total_atoms-4, total_atoms-10, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(i, total_atoms+box[1], total_atoms-10,\
                       total_atoms-4, 75, mask=mask)
    mask[-4] = 0
    atoms.set_angle(total_atoms+box[1], total_atoms-10,\
                    total_atoms-3, 109.5, mask=mask)
    atoms.set_distance(total_atoms-3, total_atoms-10, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(i, total_atoms+box[1], total_atoms-10,\
                       total_atoms-3, -179.5, mask=mask)
    mask[-3] = 0
    atoms.set_angle(total_atoms-11, total_atoms-9,\
                    total_atoms-2, 109.5, mask=mask)
    atoms.set_distance(total_atoms-2, total_atoms-9, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms+box[0], total_atoms-11,\
                       total_atoms-9, total_atoms-2, -179.5, mask=mask)
    mask[-2] = 0
    atoms.set_angle(total_atoms-11, total_atoms-9,\
                    total_atoms-1, 109.5, mask=mask)
    atoms.set_distance(total_atoms-1, total_atoms-9, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms+box[0], total_atoms-11,\
                       total_atoms-9, total_atoms-1, 300, mask=mask)
    mask[-1] = 0
    if member == 5:
        for j in range(1):
            atoms.append('C')
        for j in range(2):
            atoms.append('H')
        total_atoms = len(atoms)
        mask = [0]*(total_atoms-3) + [1]*3
        atoms.set_angle(i, total_atoms+box[1]-3, total_atoms-3, 109.5, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms+box[1]-3,\
                           1.5290, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[1]-3,\
                           total_atoms-3, 60, mask=mask)
        for k in range(10):
            atoms.set_distance(total_atoms-3, total_atoms+box[1]-3,\
                               1.5290, fix=1, mask=mask)
            atoms.set_distance(total_atoms-3, total_atoms-14, 1.5290, fix=1, mask=mask)
        mask[-3] = 0
        atoms.set_angle(total_atoms+box[1]-3, total_atoms-3, total_atoms-2,\
                        113.2, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-3, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1]-3, total_atoms-3,\
                           total_atoms-2, 185, mask=mask)
        mask[-2] = 0
        atoms.set_angle(total_atoms+box[1]-3, total_atoms-3,\
                        total_atoms-1, 113.2, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-3, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1]-3, total_atoms-3,\
                           total_atoms-1, 300, mask=mask)
    elif member == 6:
        for j in range(2):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        total_atoms = len(atoms)
        mask = [0]*(total_atoms-6) + [1]*6
        atoms.set_angle(i, total_atoms+box[1]-6, total_atoms-6, 110.5, mask=mask)
        atoms.set_distance(total_atoms-6, total_atoms+box[1]-6,\
                           1.5290, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i, total_atoms+box[1]-6,\
                           total_atoms-6, 67.5, mask=mask)
        mask[-6] = 0
        atoms.set_angle(total_atoms+box[0]-6, total_atoms-17,\
                        total_atoms-5, 110.5, mask=mask)
        atoms.set_distance(total_atoms-5, total_atoms-17, 1.5290, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[0]-6, total_atoms-17,\
                           total_atoms-5, 67.5, mask=mask)
        for k in range(10):
            atoms.set_distance(total_atoms-5, total_atoms-6, 1.5290, fix=1, mask=mask)
            atoms.set_distance(total_atoms-5, total_atoms-17, 1.5290, fix=1, mask=mask)
        mask[-5] = 0
        atoms.set_angle(total_atoms+box[1]-6, total_atoms-6,\
                        total_atoms-4, 110.2, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-6, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1]-6,\
                           total_atoms-6, total_atoms-4, 185, mask=mask)
        mask[-4] = 0
        atoms.set_angle(total_atoms+box[1]-6, total_atoms-6,\
                        total_atoms-3, 110.2, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-6, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(i, total_atoms+box[1]-6,\
                           total_atoms-6, total_atoms-3, 300, mask=mask)
        mask[-3] = 0
        atoms.set_angle(total_atoms-17, total_atoms-5, total_atoms-2, 110.2, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-5, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms+box[0]-6, total_atoms-17,\
                           total_atoms-5, total_atoms-2, 185, mask=mask)
        mask[-2] = 0
        atoms.set_angle(total_atoms-17, total_atoms-5, total_atoms-1, 110.2, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-5, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms+box[0]-6, total_atoms-17,\
                           total_atoms-5, total_atoms-1, 60, mask=mask)
    return atoms
    
def Bpin(atoms, i, state=1, shuffle=0):
    atoms.append('O')
    atoms.append('O')
    total_atoms = len(atoms)
    elements = atoms.get_chemical_symbols()
    elements[i] = 'B'
    atoms.set_chemical_symbols(elements)
    Cnear = near(i, atoms)
    mask = [0]*(total_atoms-2) + [1]*2
    CnearH = near(Cnear, atoms, 'H')
    x = atoms.get_positions()
    x[-2][0] = x[-2][0] + 50
    x[-2][1] = x[-2][1] + 50
    x[-2][2] = x[-2][2] + 50
    x[-1][0] = x[-1][0] - 50
    x[-1][1] = x[-1][1] - 50
    x[-1][2] = x[-1][2] - 50
    if state == 1:
        atoms.set_distance(i, Cnear, 1.5950, fix=1)
    atoms.set_distance(total_atoms-2, i, 1.4700, fix=1, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-2, 120.0, mask=mask)
    if shuffle == 0:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-2, 20, mask=mask)
    elif shuffle == 1:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-2, -20, mask=mask)
    elif shuffle == 5:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-2, 0.5, mask=mask)
    mask[-2] = 0
    atoms.set_distance(total_atoms-1, i, 1.4700, fix=1, mask=mask)
    atoms.set_angle(Cnear, i, total_atoms-1, 120.0, mask=mask)
    if shuffle == 0:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-1, 200, mask=mask)
    elif shuffle == 1:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-1, 160, mask=mask)
    elif shuffle == 5:
        atoms.set_dihedral(CnearH, Cnear, i, total_atoms-1, 180.5, mask=mask)
    atoms.append('C')
    atoms.append('C')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-2) + [1]*2
    x = atoms.get_positions()
    x[-2][0] = x[-2][0] + 50
    x[-2][1] = x[-2][1] + 50
    x[-2][2] = x[-2][2] + 50
    x[-1][0] = x[-1][0] - 50
    x[-1][1] = x[-1][1] - 50
    x[-1][2] = x[-1][2] - 50
    atoms.set_distance(total_atoms-2, total_atoms-4, 1.4100, fix=1, mask=mask)
    atoms.set_angle(i, total_atoms-4, total_atoms-2, 95.0, mask=mask)
    atoms.set_dihedral(Cnear, i, total_atoms-4, total_atoms-2, 125, mask=mask)
    mask[-2] = 0
    atoms.set_dihedral(Cnear, i, total_atoms-3, total_atoms-1, 235, mask=mask)
    for k in range(10):
        atoms.set_distance(total_atoms-1, total_atoms-3, 1.4100, fix=1, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-2, 1.5290, fix=1, mask=mask)
    for i in range(4):
        atoms.append('H')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-4) + [1]*4
    x = atoms.get_positions()
    x[-4][0] = x[-4][0] - 50
    x[-4][1] = x[-4][1] - 50
    x[-4][2] = x[-4][2] - 50
    x[-3][0] = x[-3][0] + 50
    x[-3][1] = x[-3][1] + 50
    x[-3][2] = x[-3][2] + 50
    x[-2][0] = x[-2][0] - 50
    x[-2][1] = x[-2][1] - 50
    x[-2][2] = x[-2][2] - 50
    x[-1][0] = x[-1][0] + 50
    x[-1][1] = x[-1][1] + 50
    x[-1][2] = x[-1][2] + 50
    atoms.set_dihedral(i, total_atoms-8, total_atoms-6,\
                       total_atoms-4, 120, mask=mask)
    atoms.set_angle(total_atoms-8, total_atoms-6, total_atoms-4, 111.5, mask=mask)
    atoms.set_distance(total_atoms-4, total_atoms-6, 1.0900, fix=1)
    mask[-4] = 0
    atoms.set_dihedral(i, total_atoms-8, total_atoms-6,\
                       total_atoms-3, 240, mask=mask)
    atoms.set_angle(total_atoms-8, total_atoms-6, total_atoms-3, 111.5, mask=mask)
    atoms.set_distance(total_atoms-3, total_atoms-6, 1.0900, fix=1)
    mask[-3] = 0
    atoms.set_dihedral(i, total_atoms-7, total_atoms-5,\
                       total_atoms-2, 120, mask=mask)
    atoms.set_angle(total_atoms-7, total_atoms-5, total_atoms-2, 111.5, mask=mask)
    atoms.set_distance(total_atoms-2, total_atoms-5, 1.0900, fix=1)
    mask[-2] = 0
    atoms.set_dihedral(i, total_atoms-7, total_atoms-5,\
                       total_atoms-1, 240, mask=mask)
    atoms.set_angle(total_atoms-7, total_atoms-5, total_atoms-1, 111.5, mask=mask)
    atoms.set_distance(total_atoms-1, total_atoms-5, 1.0900, fix=1)
    if (atoms.get_distance(total_atoms-4, total_atoms-5) <= 1.3 or
        atoms.get_distance(total_atoms-3, total_atoms-5) <= 1.3 or 
        atoms.get_distance(total_atoms-2, total_atoms-6) <= 1.3 or
        atoms.get_distance(total_atoms-1, total_atoms-6)<=1.3):
        mask = [0]*(total_atoms-4) + [1]*4
        x = atoms.get_positions()
        x[-4][0] = x[-4][0] - 50
        x[-4][1] = x[-4][1] - 50
        x[-4][2] = x[-4][2] - 50
        x[-3][0] = x[-3][0] + 50
        x[-3][1] = x[-3][1] + 50
        x[-3][2] = x[-3][2] + 50
        x[-2][0] = x[-2][0] - 50
        x[-2][1] = x[-2][1] - 50
        x[-2][2] = x[-2][2] - 50
        x[-1][0] = x[-1][0] + 50
        x[-1][1] = x[-1][1] + 50
        x[-1][2] = x[-1][2] + 50
        atoms.set_dihedral(i, total_atoms-8, total_atoms-6,\
                           total_atoms-4, 300, mask=mask)
        atoms.set_angle(total_atoms-8, total_atoms-6, total_atoms-4, 111.5, mask=mask)
        atoms.set_distance(total_atoms-4, total_atoms-6, 1.0900, fix=1)
        mask[-4] = 0
        atoms.set_dihedral(i, total_atoms-8, total_atoms-6,\
                           total_atoms-3, 179.5, mask=mask)
        atoms.set_angle(total_atoms-8, total_atoms-6, total_atoms-3, 111.5, mask=mask)
        atoms.set_distance(total_atoms-3, total_atoms-6, 1.0900, fix=1)
        mask[-3] = 0
        atoms.set_dihedral(i, total_atoms-7, total_atoms-5,\
                           total_atoms-2, 60, mask=mask)
        atoms.set_angle(total_atoms-7, total_atoms-5, total_atoms-2, 111.5, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-5, 1.0900, fix=1)
        mask[-2] = 0
        atoms.set_dihedral(i, total_atoms-7, total_atoms-5,\
                           total_atoms-1, 179.5, mask=mask)
        atoms.set_angle(total_atoms-7, total_atoms-5, total_atoms-1, 111.5, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-5, 1.0900, fix=1)
    return atoms
    
def near(i, atoms, target='ALL', ignore=[], limit=True, reverse=False):
    elements = atoms.get_chemical_symbols()
    total_atoms = len(atoms)
    j = 0
    targeted_distances=[]
    for element in elements:
        if j in ignore:
            j = j + 1
            continue
        if element == target:
            d = atoms.get_distances(i, j)
            if limit == True:
                if d[0] >= 0.5:
                    targeted_distances.append(d[0])
            else:
                targeted_distances.append(d[0])
        if target == 'Not H' or target == 'notH' or target == 'not H':
            if element != 'H':
                d = atoms.get_distances(i, j)
                if limit == True:
                    if d[0] >= 0.5:
                        targeted_distances.append(d[0])
                else:
                    targeted_distances.append(d[0])
        j = j + 1
    if reverse == True:
        _near = -1
    else:
        _near = 0
    if target == 'all' or target == 'All' or target == 'ALL':
        targeted_distances = atoms.get_distances(i, range(len(atoms)))
        targeted_distances = np.delete(targeted_distances,\
                                       np.where(targeted_distances == 0.0))
    distances = atoms.get_distances(i, range(len(atoms)))
    #print(distances, targeted_distances)
    for distance in distances:
        if total_atoms + _near in ignore:
            continue
        if reverse == True:
            if total_atoms + _near in ignore:
                _near = _near - 1
                continue
            if abs(min(targeted_distances) - atoms.get_distances(i, total_atoms+_near)) <= 0.001:
                #print(_near, atoms.get_distances(i, _near))
                return total_atoms + _near
            _near = _near - 1
        else:
            if _near in ignore:
                _near = _near + 1
                continue
            if abs(min(targeted_distances) - atoms.get_distances(i, _near)) <= 0.001:
                #print(_near, atoms.get_distances(i, _near))
                return _near
            #print(min(targeted_distances), atoms.get_distances(i, _near))
            _near = _near + 1