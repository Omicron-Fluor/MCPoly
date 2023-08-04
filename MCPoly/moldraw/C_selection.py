import math as m
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.visualize import view

def C_selection(substitution, atoms, Cnear, state=[1, 1]):
    if substitution in ('CO', 'C=O'):
        return CO(atoms, Cnear, state)
    if substitution in ('CNH', 'C=NH'):
        return CNH(atoms, Cnear, state)
    if substitution == 'O':
        return O(atoms, Cnear, state)
    if substitution == 'S':
        return S(atoms, Cnear, state)
    if substitution == 'N':
        return N(atoms, Cnear, state)
                        
def CO(atoms, Cnear, state=[1, 1]):
    atoms.append('O')
    total_atoms = len(atoms)
    x = atoms.get_positions()
    x[-1][0] = x[-1][0] - 50
    x[-1][1] = x[-1][1] - 50
    x[-1][2] = x[-1][2] - 50
    atoms.set_positions(x)
    CnearC1 = near(Cnear, atoms, 'C')
    CnearC2 = near(Cnear, atoms, 'C', ignore=[CnearC1], reverse=True)
    #print(CnearC1, CnearC2)
    mask = [0]*(total_atoms-1)+[1]
    mask[Cnear] = 1
    for k in range(10):
        if state[0] == 1:
            atoms.set_distance(Cnear, CnearC1, 1.5220, fix=1, mask=mask)
        elif state[0] == 2:
            atoms.set_distance(Cnear, CnearC1, 1.4440, fix=1, mask=mask)
        elif state[0] == 6:
            atoms.set_distance(Cnear, CnearC1, 1.4900, fix=1, mask=mask)
        elif state[0] == 7:
            atoms.set_distance(Cnear, CnearC1, 1.3350, fix=1, mask=mask)
        elif state[0] == 8:
            atoms.set_distance(Cnear, CnearC1, 1.3640, fix=1, mask=mask)
        elif state[0] == 82:
            atoms.set_distance(Cnear, CnearC1, 1.5100, fix=1, mask=mask)
        if state[1] == 1:
            atoms.set_distance(Cnear, CnearC2, 1.5220, fix=1, mask=mask)
        elif state[1] == 2:
            atoms.set_distance(Cnear, CnearC2, 1.4440, fix=1, mask=mask)
        elif state[1] == 6:
            atoms.set_distance(Cnear, CnearC2, 1.4900, fix=1, mask=mask)
        elif state[1] == 7:
            atoms.set_distance(Cnear, CnearC2, 1.3350, fix=1, mask=mask)
        elif state[1] == 8:
            atoms.set_distance(Cnear, CnearC2, 1.3640, fix=1, mask=mask)
        elif state[1] == 82:
            atoms.set_distance(Cnear, CnearC2, 1.5100, fix=1, mask=mask)
    mask[Cnear] = 0
    for k in range(10):
        atoms.set_angle(CnearC1, Cnear, total_atoms-1, 120, mask=mask)
        atoms.set_angle(CnearC2, Cnear, total_atoms-1, 120, mask=mask)
    atoms.set_distance(total_atoms-1, Cnear, 1.2290, fix=1, mask=mask)
    for i in range(2):
        Hnear = near(Cnear, atoms, 'H', limit = False)
        del atoms[Hnear]
    return atoms

def CNH(atoms, Cnear, state=[1, 1]):
    atoms.append('N')
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
    CnearC1 = near(Cnear, atoms, 'C')
    CnearC2 = near(Cnear, atoms, 'C', ignore=[CnearC1], reverse=True)
    #print(CnearC1, CnearC2)
    mask = [0]*(total_atoms-2) + [2]
    mask[Cnear] = 1
    for k in range(10):
        if state[0] == 1:
            atoms.set_distance(Cnear, CnearC1, 1.5220, fix=1, mask=mask)
        elif state[0] == 2:
            atoms.set_distance(Cnear, CnearC1, 1.5220, fix=1, mask=mask)
        elif state[0] == 6:
            atoms.set_distance(Cnear, CnearC1, 1.5220, fix=1, mask=mask)
        elif state[0] == 7:
            atoms.set_distance(Cnear, CnearC1, 1.3350, fix=1, mask=mask)
        elif state[0] == 8:
            atoms.set_distance(Cnear, CnearC1, 1.3640, fix=1, mask=mask)
        elif state[0] == 82:
            atoms.set_distance(Cnear, CnearC1, 1.5100, fix=1, mask=mask)
        if state[1] == 1:
            atoms.set_distance(Cnear, CnearC2, 1.5220, fix=1, mask=mask)
        elif state[1] == 2:
            atoms.set_distance(Cnear, CnearC2, 1.5220, fix=1, mask=mask)
        elif state[1] == 6:
            atoms.set_distance(Cnear, CnearC2, 1.5220, fix=1, mask=mask)
        elif state[1] == 7:
            atoms.set_distance(Cnear, CnearC2, 1.3350, fix=1, mask=mask)
        elif state[1] == 8:
            atoms.set_distance(Cnear, CnearC2, 1.3640, fix=1, mask=mask)
        elif state[1] == 82:
            atoms.set_distance(Cnear, CnearC2, 1.5100, fix=1, mask=mask)
    mask[Cnear] = 0
    for k in range(10):
        atoms.set_angle(CnearC1, Cnear, total_atoms-2, 120, mask=mask)
        atoms.set_angle(CnearC2, Cnear, total_atoms-2, 120, mask=mask)
    mask[-2] = 0
    atoms.set_angle(Cnear, total_atoms-2, total_atoms-1, 120.5, mask=mask)
    atoms.set_distance(total_atoms-1, total_atoms-2, 1.0100, fix=1, mask=mask)
    atoms.set_dihedral(Cnear1, Cnear, total_atoms-2,\
                       total_atoms-1, 179.5, mask=mask)
    for i in range(2):
        Hnear = near(Cnear, atoms, 'H', limit = False)
        del atoms[Hnear]
    return atoms

def O(atoms, Cnear, state=[1, 1]):
    atoms.append('O')
    total_atoms = len(atoms)
    x = atoms.get_positions()
    x[-1][0] = x[-1][0] - 50
    x[-1][1] = x[-1][1] - 50
    x[-1][2] = x[-1][2] - 50
    atoms.set_positions(x)
    CnearC1 = near(Cnear, atoms, 'C')
    CnearC2 = near(Cnear, atoms, 'C', ignore=[CnearC1], reverse=True)
    print(CnearC1, CnearC2)
    mask = [0]*(total_atoms-1) + [1]
    atoms.set_distance(total_atoms-1, Cnear, 0, fix=1)
    for k in range(10):
        if state[0] == 1:
            atoms.set_distance(total_atoms-1, CnearC1, 1.410, fix=1)
        elif state[0] == 2 or state[0] == 6 or state[0] == 82:
            atoms.set_distance(total_atoms-1, CnearC1, 1.364, fix=1)
        elif state[0] == 4:
            atoms.set_distance(total_atoms-1, CnearC1, 1.327, fix=1)
        elif state[0] == 164:
            atoms.set_distance(total_atoms-1, CnearC1, 1.670, fix=1)
        if state[1] == 1:
            atoms.set_distance(total_atoms-1, CnearC2, 1.410, fix=1)
        elif state[1] == 2 or state[1] == 6 or state[1] == 82:
            atoms.set_distance(total_atoms-1, CnearC2, 1.364, fix=1)
        elif state[1] == 4:
            atoms.set_distance(total_atoms-1, CnearC2, 1.327, fix=1)
        elif state[1] == 164:
            atoms.set_distance(total_atoms-1, CnearC2, 1.670, fix=1)
    del atoms[Cnear]
    for i in range(2):
        total_atoms = len(atoms)
        Hnear = near(total_atoms-1, atoms, 'H', limit=False)
        del atoms[Hnear]
    return atoms

def S(atoms, Cnear, state=[1, 1]):
    atoms.append('S')
    total_atoms = len(atoms)
    x = atoms.get_positions()
    x[-1][0] = x[-1][0] - 50
    x[-1][1] = x[-1][1] - 50
    x[-1][2] = x[-1][2] - 50
    atoms.set_positions(x)
    CnearC1 = near(Cnear, atoms, 'C')
    CnearC2 = near(Cnear, atoms, 'C', ignore=[CnearC1], reverse=True)
    print(CnearC1, CnearC2)
    mask = [0]*(total_atoms-1) + [1]
    atoms.set_distance(total_atoms-1, Cnear, 0, fix=1)
    for k in range(10):
        if state[0] == 1:
            atoms.set_distance(total_atoms-1, CnearC1, 1.810, fix=1)
        elif state[0] == 2:
            atoms.set_distance(total_atoms-1, CnearC1, 1.770, fix=1)
        if state[1] == 1:
            atoms.set_distance(total_atoms-1, CnearC2, 1.810, fix=1)
        elif state[1] == 2:
            atoms.set_distance(total_atoms-1, CnearC2, 1.770, fix=1)
    Hnear = near(total_atoms-1, atoms, 'C', limit=False)
    del atoms[Hnear]
    for i in range(2):
        total_atoms = len(atoms)
        Hnear = near(total_atoms-1, atoms, 'H', limit=False)
        del atoms[Hnear]
    return atoms

def N(atoms, Cnear, state=[1, 1]):
    atoms.append('N')
    total_atoms = len(atoms)
    x = atoms.get_positions()
    x[-1][0] = x[-1][0] - 50
    x[-1][1] = x[-1][1] - 50
    x[-1][2] = x[-1][2] - 50
    atoms.set_positions(x)
    CnearC1 = near(Cnear, atoms, 'C')
    CnearC2 = near(Cnear, atoms, 'C', ignore=[CnearC1], reverse=True)
    print(CnearC1, CnearC2)
    mask = [0]*(total_atoms-1) + [1]
    atoms.set_distance(total_atoms-1, Cnear, 0, fix=1)
    for k in range(10):
        if state[0] == 1:
            atoms.set_distance(total_atoms-1, CnearC1, 1.4480, fix=1)
        elif state[0] == 2 or state[0] == 6 or state[0] == 82:
            atoms.set_distance(total_atoms-1, CnearC1, 1.3400, fix=1)
        elif state[0] == 4:
            atoms.set_distance(total_atoms-1, CnearC1, 1.4900, fix=1)
        if state[1] == 1:
            atoms.set_distance(total_atoms-1, CnearC2, 1.4480, fix=1)
        elif state[1] == 2 or state[1] == 6 or state[1] == 82:
            atoms.set_distance(total_atoms-1, CnearC2, 1.3400, fix=1)
        elif state[1] == 4:
            atoms.set_distance(total_atoms-1, CnearC2, 1.4900, fix=1)
    CnearC1 = numjudge(CnearC1, Cnear)
    CnearC2 = numjudge(CnearC2, Cnear)
    del atoms[Cnear]
    for i in range(2):
        total_atoms = len(atoms)
        Hnear = near(total_atoms-1, atoms, 'H', limit=False)
        del atoms[Hnear]
        CnearC1 = numjudge(CnearC1, Hnear)
        CnearC2 = numjudge(CnearC2, Hnear)
    atoms.append('H')
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-1) + [1]
    CnearC3 = near(CnearC1, atoms, 'C', ignore = [CnearC1, CnearC2])
    for k in range(10):
        atoms.set_angle(CnearC1, total_atoms-2, total_atoms-1, 109.5, mask=mask)
        atoms.set_angle(CnearC2, total_atoms-2, total_atoms-1, 109.5, mask=mask)
        atoms.set_dihedral(CnearC3, CnearC1, total_atoms-2,\
                           total_atoms-1, 150, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-2, 1.0100, fix=1)
    return atoms