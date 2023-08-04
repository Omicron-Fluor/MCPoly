import math as m
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.visualize import view

def bind_selection(ring_size, atoms, bond, state=[1, 1], place=0, strict={1:'X'}, side=0, shuffle=0, less5=False, stable=False):
    '''This function will bind the ring onto a bond. 5 means adding a 5-membered-ring bound. 
    Adding a 6 in the front will make sure the conjugation of the ring.''' 
    
    if ring_size in (2, 'bind2'):
        return bind2(atoms, bond, state)
    elif ring_size in (3, 'bind3'):
        return bind3(atoms, bond, state, strict)
    elif ring_size in (5, 'bind5'):
        return bind5(atoms, bond, state, strict, side, shuffle, less5, stable)
    elif ring_size in (6, 'bind6'):
        return bind6(atoms, bond, state, strict, side, shuffle, less5)
    elif ring_size in (65, 'bindconj5'):
        return bindconj5(atoms, bond, state, strict)
    elif ring_size in (66, 'bindconj6'):
        return bindconj6(atoms, bond, state, strict)

def bind2(atoms, bond, state=[1, 1]):
    i1 = bond[0]
    i2 = bond[1]
    elements = atoms.get_chemical_symbols()
    element1 = elements[i1]
    element2 = elements[i2]
    Cnear1 = near(i1, atoms, 'Not H', ignore=[i1, i2])
    Cnear2 = near(i2, atoms, 'Not H', ignore=[i1, i2], reverse=True)
    total_atoms = len(atoms)
    mask = [0] * total_atoms
    mask[i1] = 1
    mask[i2] = 1
    bondjudge2([Cnear1, Cnear2], atoms, [i1, i2], [element1, element2],\
               state, mask)
    for j in range(2):
        CnearH = near(i1, atoms, 'H')
        del atoms[CnearH]
        i1 = numjudge(i1, CnearH)
        i2 = numjudge(i2, CnearH)
    for j in range(2):
        CnearH = near(i2, atoms, 'H')
        del atoms[CnearH]
        i1 = numjudge(i1, CnearH)
        i2 = numjudge(i2, CnearH)
    atoms.append('H')
    atoms.append('H')
    x = atoms.get_positions()
    x[-2][0] = x[-2][0] - 100
    x[-2][1] = x[-2][1] - 50
    x[-2][2] = x[-2][2] - 50
    x[-1][0] = x[-1][0] + 50
    x[-1][1] = x[-1][1] + 50
    x[-1][2] = x[-1][2] + 50
    atoms.set_positions(x)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-2) + [1]*2
    for k in range(10):
        atoms.set_angle(Cnear1, i1, total_atoms-2, 118, mask=mask)
        atoms.set_angle(i2, i1, total_atoms-2, 118, mask=mask)
    atoms.set_distance(total_atoms-2, i1, 1.0800, fix=1, mask=mask)
    mask[-2] = 0
    for k in range(10):
        atoms.set_angle(Cnear2, i2, total_atoms-1, 118, mask=mask)
        atoms.set_angle(i1, i2, total_atoms-1, 118, mask=mask)
    atoms.set_distance(total_atoms-1, i2, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-2, i1, i2, total_atoms-1, 0.5, mask=mask)
    if element1 == 'N':
        del atoms[-2]
    if element2 == 'N':
        del atoms[-1]
    return atoms

def bind3(atoms, bond, state=[1, 1], strict={1:'X'}):
    i1 = bond[0]
    i2 = bond[1]
    Cnear = near(i1, atoms, 'C')
    total_atoms = len(atoms)
    if strict[1] == 'O':
        atoms.append('O')
        mask = [0]*total_atoms + [1]
        if state == [1, 1]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.4100, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.4100, fix=1, mask=mask)
        elif state == [1, 2]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.4100, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.3640, fix=1, mask=mask)
        elif state == [2, 1]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.3640, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.4100, fix=1, mask=mask)
        elif state == [2, 2] or state == [6, 6]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.3640, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.3640, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i1, i2, total_atoms, 75, mask=mask)
        for j in range(2):
            CnearH = near(total_atoms - j, atoms, 'H')
            del atoms[CnearH]
    elif strict[1] == 'S':
        atoms.append('S')
        mask = [0]*total_atoms + [1]
        if state == [1, 1]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.8100, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.8100, fix=1, mask=mask)
        elif state == [1, 2]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.8100, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.7700, fix=1, mask=mask)
        elif state == [2, 1]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.7700, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.8100, fix=1, mask=mask)
        elif state == [2, 2] or state == [6, 6]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.7700, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.7700, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i1, i2, total_atoms, 75, mask=mask)
        for j in range(2):
            CnearH = near(total_atoms, atoms, 'H')
            del atoms[CnearH]
            total_atoms = numjudge(total_atoms, CnearH)
            i1 = numjudge(i1, CnearH)
            i2 = numjudge(i2, CnearH)
    elif strict[1] == 'N':
        atoms.append('N')
        atoms.append('H')
        mask = [0]*total_atoms + [1]*2
        if state == [1, 1]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.4480, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.4480, fix=1, mask=mask)
        elif state == [1, 2]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.4480, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.3400, fix=1, mask=mask)
        elif state == [2, 1]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.3400, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.4480, fix=1, mask=mask)
        elif state == [2, 2] or state == [6, 6]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.3400, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.3400, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i1, i2, total_atoms, 75, mask=mask)
        mask[-2] = 0
        for j in range(2):
            CnearH = near(total_atoms, atoms, 'H')
            del atoms[CnearH]
            total_atoms = numjudge(total_atoms, CnearH)
            i1 = numjudge(i1, CnearH)
            i2 = numjudge(i2, CnearH)
        x = atoms.get_positions()
        x[-1][0] = x[-1][0] + 50
        x[-1][1] = x[-1][1] + 50
        x[-1][2] = x[-1][2] + 50
        atoms.set_positions(x)
        total_atoms = len(atoms)
        mask = [0]*(total_atoms-1) + [1]
        atoms.set_distance(total_atoms-1, total_atoms-2, 0, fix=1, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-2, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-1, total_atoms-2, i1, i2, -90.5, mask=mask)
        for k in range(10):
            atoms.set_angle(i1, total_atoms-2, total_atoms-1, 111.5, mask=mask)
            atoms.set_angle(i2, total_atoms-2, total_atoms-1, 111.5, mask=mask)
    else:
        atoms.append('C')
        atoms.append('H')
        atoms.append('H')
        mask = [0]*total_atoms + [1]*3
        if state == [1, 1]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.5290, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.5290, fix=1, mask=mask)
        elif state == [1, 2]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.5290, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.5095, fix=1, mask=mask)
        elif state == [2, 1]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.5095, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.5290, fix=1, mask=mask)
        elif state == [2, 2]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.5095, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.5095, fix=1, mask=mask)
        elif state == [6, 6]:
            for k in range(10):
                atoms.set_distance(total_atoms, i1, 1.5105, fix=1, mask=mask)
                atoms.set_distance(total_atoms, i2, 1.5105, fix=1, mask=mask)
        atoms.set_dihedral(Cnear, i1, i2, total_atoms, 75, mask=mask)
        mask[-3] = 0
        for j in range(2):
            CnearH = near(total_atoms, atoms, 'H')
            del atoms[CnearH]
            total_atoms = numjudge(total_atoms, CnearH)
            i1 = numjudge(i1, CnearH)
            i2 = numjudge(i2, CnearH)
        x = atoms.get_positions()
        x[-2][0] = x[-2][0] + 50
        x[-2][1] = x[-2][1] + 50
        x[-2][2] = x[-2][2] + 50
        x[-1][0] = x[-1][0] - 50
        x[-1][1] = x[-1][1] - 50
        x[-1][2] = x[-1][2] - 50
        atoms.set_positions(x)
        total_atoms = len(atoms)
        mask = [0]*(total_atoms-2) + [1]*2
        atoms.set_distance(total_atoms-2, total_atoms-3, 0, fix=1, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-3, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-2, total_atoms-3, i1, Cnear, 110, mask=mask)
        for k in range(10):
            atoms.set_angle(i1, total_atoms-3, total_atoms-2, 111.5, mask=mask)
            atoms.set_angle(i2, total_atoms-3, total_atoms-2, 111.5, mask=mask)
        mask[-2] = 0
        atoms.set_distance(total_atoms-1, total_atoms-3, 0, fix=1, mask=mask)
        atoms.set_distance(total_atoms-1, total_atoms-3, 1.0900, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-1, total_atoms-3, i1, Cnear, -70.5, mask=mask)
        for k in range(10):
            atoms.set_angle(i1, total_atoms-3, total_atoms-1, 111.5, mask=mask)
            atoms.set_angle(i2, total_atoms-3, total_atoms-1, 111.5, mask=mask)
        atoms.set_dihedral(total_atoms-1, total_atoms-3, i1, i2, -90.5, mask=mask)
    return atoms

def bind5(atoms, bond, state=[1, 1], strict={1:'X'}, side=0, shuffle=0, less5=False, stable=False):
    i1 = bond[0]
    i2 = bond[1]
    if side == 0:
        Hnear1 = near(i1, atoms, 'H')
        Hnear2 = near(i2, atoms, 'H')
    elif side == 1:
        Hnear1 = near(i1, atoms, 'H', reverse=True)
        Hnear2 = near(i2, atoms, 'H', reverse=True)
    elif side == 2:
        Hnear1 = near(i1, atoms, 'H')
        Hnear2 = near(i2, atoms, 'H', reverse=True)
    elif side == 3:
        Hnear1 = near(i1, atoms, 'H', reverse=True)
        Hnear2 = near(i2, atoms, 'H')
    print(Hnear1, Hnear2)
    total_atoms = len(atoms)
    mask = [0] * total_atoms
    mask[Hnear1] = 1
    atoms.set_angle(i2, i1, Hnear1, 100, mask=mask)
    mask[Hnear1] = 0
    mask[Hnear2] = 1
    atoms.set_angle(i1, i2, Hnear2, -100, mask=mask)
    if (less5 == True and side < 2) or (less5 == False and side >= 2):
        atoms.set_dihedral(Hnear1, i1, i2, Hnear2, -35, mask=mask)
        xoy = near(Hnear2, atoms, limit=False)
        if atoms.get_distance(xoy, Hnear2) <= 0.9:
            atoms.set_dihedral(Hnear1, i1, i2, Hnear2, 35, mask=mask)
        xoy = near(i2, atoms, target='Not H')
        if xoy == i1:
            xoy = near(i2, atoms, target='Not H', reverse=True)
        if atoms.get_angle(xoy, i2, Hnear2) <= 90:
            atoms.set_dihedral(Hnear1, i1, i2, Hnear2, 35, mask=mask)
    elif less5 == False and side < 2:
        atoms.set_dihedral(Hnear1, i1, i2, Hnear2, 55, mask=mask)
        xoy = near(Hnear2, atoms, limit=False)
        if atoms.get_distance(xoy, Hnear2) <= 0.9:
            atoms.set_dihedral(Hnear1, i1, i2, Hnear2, -55, mask=mask)
        xoy = near(i2, atoms, target='Not H')
        if xoy == i1:
            xoy = near(i2, atoms, target='Not H', reverse=True)
        if atoms.get_angle(xoy, i2, Hnear2) <= 90:
            atoms.set_dihedral(Hnear1, i1, i2, Hnear2, -55, mask=mask)
    elif less5 == True and side >= 2:
        atoms.set_dihedral(Hnear1, i1, i2, Hnear2, 80, mask=mask)
        xoy = near(Hnear2, atoms, limit=False)
        if atoms.get_distance(xoy, Hnear2) <= 0.9:
            atoms.set_dihedral(Hnear1, i1, i2, Hnear2, -80, mask=mask)
        xoy = near(i2, atoms, target='Not H')
        if xoy == i1:
            xoy = near(i2, atoms, target='Not H', reverse=True)
        if atoms.get_angle(xoy, i2, Hnear2) <= 90:
            atoms.set_dihedral(Hnear1, i1, i2, Hnear2, -80, mask=mask)
        mask[Hnear2] = 0
        mask[Hnear1] = 1
        if atoms.get_dihedral(Hnear1, i1, i2, Hnear2)<85:
            atoms.set_dihedral(Hnear2, i2, i1, Hnear1, 50, mask=mask)
        else:
            atoms.set_dihedral(Hnear2, i2, i1, Hnear1, 360 - 50, mask=mask)
    if less5 == True and side >= 2:
        elements = atoms.get_chemical_symbols()
        elements[Hnear1] = 'C'
        elements[Hnear2] = 'C'
        atoms.set_chemical_symbols(elements)
        atoms.set_distance(Hnear1, i1, 1.0400, fix=1)
        atoms.set_distance(Hnear2, i2, 1.0400, fix=1)
        Hnear3 = near(i1, atoms, 'H')
        Hnear4 = near(i2, atoms, 'H')
        atoms.set_distance(Hnear1, i1, 1.0900, fix=1)
        atoms.set_distance(Hnear2, i2, 1.0900, fix=1)
        elements[Hnear1] = 'H'
        elements[Hnear2] = 'H'
        atoms.set_chemical_symbols(elements)
    #    x = atoms.get_positions()
    #    x[Hnear3][0] = x[Hnear3][0] - 50
    #    x[Hnear3][1] = x[Hnear3][1] - 50
    #    x[Hnear3][2] = x[Hnear3][2] - 50
    #    x[Hnear4][0] = x[Hnear4][0] + 50
    #    x[Hnear4][1] = x[Hnear4][1] + 50
    #    x[Hnear4][2] = x[Hnear4][2] + 50
    #    atoms.set_positions(x)
    if less5 == True and side >= 2:
        Cnear1 = near(i1, atoms, 'Not H')
        Cnear2 = near(i2, atoms, 'Not H')
        mask[Hnear1] = 0
        mask[Hnear2] = 0
        mask[Hnear3] = 1
        mask[Hnear4] = 1
        if Cnear1 == i2:
            Cnear1 = near(i1, atoms, 'Not H', reverse=True)
        if Cnear2 == i1:
            Cnear2 = near(i2, atoms, 'Not H', reverse=True)
        for k in range(10):
            atoms.set_angle(i2, i1, Hnear3, 105.0, mask=mask)
            atoms.set_angle(Hnear1, i1, Hnear3, 105.0, mask=mask)
        atoms.set_distance(Hnear3, i1, 1.0900, fix=1, mask=mask)
        mask[Hnear3] = 0
        for k in range(10):
            atoms.set_angle(i1, i2, Hnear4, 105.0, mask=mask)
            atoms.set_angle(Hnear2, i2, Hnear4, 105.0, mask=mask)
        atoms.set_distance(Hnear4, i2, 1.0900, fix=1, mask=mask)
    if shuffle == 0:
        atoms = CH3(atoms, Hnear1, state[0])
        atoms = CH3(atoms, Hnear2, state[1], ring=True)
    elif shuffle == 1:
        atoms = CH3(atoms, Hnear1, state[0], ring=True)
        atoms = CH3(atoms, Hnear2, state[1])
    elif shuffle == 2:
        atoms = CH3(atoms, Hnear1, state[0])
        atoms = CH3(atoms, Hnear2, state[1])
    elif shuffle == 3:
        atoms = CH3(atoms, Hnear1, state[0], ring=True)
        atoms = CH3(atoms, Hnear2, state[1], ring=True)
    atoms.append('C')
    atoms.append('H')
    atoms.append('H')
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
    atoms.set_positions(x)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-3)+[1]*3
    for k in range(10):
        if atoms.get_distance(total_atoms-3, i1) < 1.6 or stable == True:
            atoms.set_distance(total_atoms-3, i1, 2.4, fix=1, mask=mask)
        if atoms.get_distance(total_atoms-3, i2) < 1.6 or stable == True:
            atoms.set_distance(total_atoms-3, i2, 2.4, fix=1, mask=mask)
        atoms.set_distance(total_atoms-3, Hnear1, 1.5290, fix=1, mask=mask)
        atoms.set_distance(total_atoms-3, Hnear2, 1.5290, fix=1, mask=mask)
    mask[-3] = 0
    for j in range(2):
        total_atoms = len(atoms)
        CnearH = near(total_atoms-3, atoms, 'H', limit=False)
        del atoms[CnearH]
        i1 = numjudge(i1, CnearH)
        i2 = numjudge(i2, CnearH)
        Hnear1 = numjudge(Hnear1, CnearH)
        Hnear2 = numjudge(Hnear2, CnearH)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-2) + [1]*2
    if side < 2:
        atoms.set_distance(total_atoms-2, total_atoms-3, 0, fix=1, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-3, 1.0900, fix=1, mask=mask)
        while(1):
            if atoms.get_distance(total_atoms-2, Hnear1) < 1.3:
                atoms.set_distance(total_atoms-2, Hnear1, 1.8, fix=1, mask=mask)
            if atoms.get_distance(total_atoms-2, Hnear2) < 1.3:
                atoms.set_distance(total_atoms-2, Hnear2, 1.8, fix=1, mask=mask)
            xoy = near(total_atoms-1, atoms, limit=False)
            #if atoms.get_distance(xoy, total_atoms-2) <= 0.85:
            #    atoms.set_distance(xoy, total_atoms-2, 1, mask=mask)
            atoms.set_angle(Hnear1, total_atoms-3, total_atoms-2, 111.0, mask=mask)
            atoms.set_angle(Hnear2, total_atoms-3, total_atoms-2, 111.0, mask=mask)
            atoms.set_distance(total_atoms-2, total_atoms-3, 1.0900, fix=1, mask=mask)
            if (atoms.get_distance(total_atoms-2, Hnear1) > 1.4 and 
                atoms.get_distance(total_atoms-2, Hnear2)) > 1.4:
                if (abs(atoms.get_angle(Hnear1, total_atoms-3, total_atoms-2)
                        - atoms.get_angle(Hnear2, total_atoms-3, total_atoms-2)) < 1):
                    #print(atoms.get_angle(Hnear1, total_atoms-3, total_atoms-2)
                    #.     - atoms.get_angle(Hnear2, total_atoms-3, total_atoms-2))
                    break
        #atoms.set_dihedral(total_atoms-2, total_atoms-3, Hnear1, i1, 60, mask=mask)
        mask[-2] = 0
        #atoms.set_distance(total_atoms-1, total_atoms-3, 0, fix=1, mask=mask)
        #atoms.set_distance(total_atoms-1, total_atoms-3, 1.0900, fix=1, mask=mask)
        #while(1):
        #    if atoms.get_distance(total_atoms-1, Hnear1) < 1.3:
        #        atoms.set_distance(total_atoms-1, Hnear1, 1.8, fix=1, mask=mask)
        #    if atoms.get_distance(total_atoms-1, Hnear2) < 1.3:
        #        atoms.set_distance(total_atoms-1, Hnear2, 1.8, fix=1, mask=mask)
        #    atoms.set_angle(Hnear1, total_atoms-3, total_atoms-1, 111.0, mask=mask)
        #    atoms.set_angle(Hnear2, total_atoms-3, total_atoms-1, 111.0, mask=mask)
        #    atoms.set_distance(total_atoms-1, total_atoms-3, 1.0900, fix=1, mask=mask)
        #    if atoms.get_distance(total_atoms-1, Hnear1) > 1.4 and
        #       atoms.get_distance(total_atoms-1, Hnear2) > 1.4:
        #        if abs(atoms.get_angle(Hnear1, total_atoms-3, total_atoms-1)
        #               - atoms.get_angle(Hnear2, total_atoms-3, total_atoms-1)) < 1:
        #            print(atoms.get_angle(Hnear1, total_atoms-3, total_atoms-1)
        #                  - atoms.get_angle(Hnear2, total_atoms-3, total_atoms-1))
        #            break
        #atoms.set_dihedral(total_atoms-1, total_atoms-3, Hnear1, i1, 170, mask=mask)
        #view(atoms)
        #mask[-2] = 0
        #mask[-1] = 1
        atoms.set_distance(total_atoms-1, total_atoms-2, 0, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-1, total_atoms-3, Hnear1, i1,\
                           atoms.get_dihedral(total_atoms-1, total_atoms-3, Hnear1, i1)+110, mask=mask)
        xoy = near(total_atoms-1, atoms, limit=False)
        if atoms.get_distance(xoy, total_atoms-1) <= 0.9:
            atoms.set_dihedral(total_atoms-1, total_atoms-3, Hnear1,\
                               i1, atoms.get_dihedral(total_atoms-1, total_atoms-3, Hnear1, i1)-220, mask=mask)
    else:
        atoms.set_distance(total_atoms-2, total_atoms-3, 0, fix=1, mask=mask)
        atoms.set_distance(total_atoms-2, total_atoms-3, 1.0900, fix=1, mask=mask)
        while(1):
            if atoms.get_distance(total_atoms-2, Hnear1) < 1.3:
                atoms.set_distance(total_atoms-2, Hnear1, 1.8, fix=1, mask=mask)
            if atoms.get_distance(total_atoms-2, Hnear2) < 1.3:
                atoms.set_distance(total_atoms-2, Hnear2, 1.8, fix=1, mask=mask)
            atoms.set_angle(Hnear1, total_atoms-3, total_atoms-2, 111.0, mask=mask)
            atoms.set_angle(Hnear2, total_atoms-3, total_atoms-2, 111.0, mask=mask)
            atoms.set_distance(total_atoms-2, total_atoms-3, 1.0900, fix=1, mask=mask)
            if (atoms.get_distance(total_atoms-2, Hnear1) > 1.4 and 
                atoms.get_distance(total_atoms-2, Hnear2) > 1.4):
                break
        mask[-2] = 0
        atoms.set_distance(total_atoms-1, total_atoms-2, 0, fix=1, mask=mask)
        #atoms.set_distance(total_atoms-1, total_atoms-3, 1.0900, fix=1, mask=mask)
        #while(1):
        #    if atoms.get_distance(total_atoms-1, Hnear1)<1.3:
        #        atoms.set_distance(total_atoms-1, Hnear1, 1.8, fix=1, mask=mask)
        #    if atoms.get_distance(total_atoms-1, Hnear2)<1.3:
        #        atoms.set_distance(total_atoms-1, Hnear2, 1.8, fix=1, mask=mask)
        #    atoms.set_angle(Hnear1, total_atoms-3, total_atoms-1, 111.0, mask=mask)
        #    atoms.set_angle(Hnear2, total_atoms-3, total_atoms-1, 111.0, mask=mask)
        #    atoms.set_distance(total_atoms-1, total_atoms-3, 1.0900, fix=1, mask=mask)
        #    if atoms.get_distance(total_atoms-1, 4)>1.4 and atoms.get_distance(total_atoms-1, 6)>1.4:
        #        break
        atoms.set_dihedral(total_atoms-1, total_atoms-3, Hnear1, i1,\
                           atoms.get_dihedral(total_atoms-1, total_atoms-3, Hnear1, i1)+110, mask=mask)
        xoy = near(total_atoms-1, atoms, limit=False)
        if atoms.get_distance(xoy, total_atoms-1) <= 0.9:
            atoms.set_dihedral(total_atoms-1, total_atoms-3, Hnear1,\
                               i1, atoms.get_dihedral(total_atoms-1, total_atoms-3, Hnear1, i1)-220, mask=mask)
    needelements = []
    places = []
    for j in range(3):
        try:
            needelements.append(strict[j+2])
            places.append(j+2)
        except:
            pass
    for i in range(len(needelements)):
        total_atoms = len(atoms)
        mask = [0]*total_atoms
        if places[i] == 2:
            mask[Hnear1] = 1
            elements = atoms.get_chemical_symbols()
            elements[Hnear1] = needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i] == 'N':
                for k in range(10):
                    if state[0] == 1:
                        atoms.set_distance(Hnear1, i1, 1.4480, fix=1, mask=mask)
                    elif state[0] == 2 or state[0] == 6:
                        atoms.set_distance(Hnear1, i1, 1.3400, fix=1, mask=mask)
                    atoms.set_distance(Hnear1, total_atoms-3, 1.4480, fix=1, mask=mask)
                mask[Hnear1] = 0
                total_atoms = len(atoms)
                CnearH = near(Hnear1, atoms, 'H', limit=False)
                del atoms[CnearH]
                Hnear1 = numjudge(Hnear1, CnearH)
                i1 = numjudge(i1, CnearH)
                i2 = numjudge(i2, CnearH)
                CnearH = near(Hnear1, atoms, 'H', limit=False)
                mask[CnearH] = 1
                for k in range(10):
                    atoms.set_angle(i1, Hnear1, CnearH, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-3, Hnear1, CnearH, 109.5, mask=mask)
                atoms.set_dihedral(i2, i1, Hnear1, CnearH, -150, mask=mask)
                atoms.set_distance(CnearH, Hnear1, 1.0100, fix=1)
            elif needelements[i] == 'O':
                for k in range(10):
                    if state[0] == 1:
                        atoms.set_distance(Hnear1, i1, 1.4100, fix=1, mask=mask)
                    elif state[0] == 2 or state[0] == 6:
                        atoms.set_distance(Hnear1, i1, 1.3640, fix=1, mask=mask)
                    atoms.set_distance(Hnear1, total_atoms-3, 1.4100, fix=1, mask=mask)
                for j in range(2):
                    total_atoms = len(atoms)
                    CnearH = near(Hnear1, atoms, 'H', limit=False)
                    del atoms[CnearH]
                    Hnear1 = numjudge(Hnear1, CnearH)
        if places[i] == 3:
            mask[total_atoms-3] = 1
            elements = atoms.get_chemical_symbols()
            elements[total_atoms-3] = needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i] == 'N':
                for k in range(10):
                    atoms.set_distance(total_atoms-3, Hnear1, 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms-3, Hnear2, 1.4480, fix=1, mask=mask)
                mask[total_atoms-3] = 0
                total_atoms = len(atoms)
                CnearH = near(total_atoms-3, atoms, 'H', limit=False)
                del atoms[CnearH]
                Hnear1 = numjudge(Hnear1, CnearH)
                i1 = numjudge(i1, CnearH)
                i2 = numjudge(i2, CnearH)
                CnearH = near(total_atoms-3, atoms, 'H', limit=False)
                mask[CnearH] = 1
                for k in range(10):
                    atoms.set_angle(Hnear1, total_atoms-3, CnearH, 109.5, mask=mask)
                    atoms.set_angle(Hnear2, total_atoms-3, CnearH, 109.5, mask=mask)
                atoms.set_dihedral(i1, Hnear1, total_atoms-3, CnearH, -150, mask=mask)
                atoms.set_distance(CnearH, total_atoms-3, 1.0100, fix=1)
            elif needelements[i] == 'O':
                for k in range(10):
                    atoms.set_distance(total_atoms-3, Hnear1, 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms-3, Hnear2, 1.4100, fix=1, mask=mask)
                for j in range(2):
                    total_atoms = len(atoms)
                    CnearH = near(total_atoms-3+j, atoms, 'H', limit=False)
                    del atoms[CnearH]
        if places[i] == 4:
            mask[Hnear2] = 1
            elements = atoms.get_chemical_symbols()
            elements[Hnear2] = needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i] == 'N':
                for k in range(10):
                    if state[0] == 1:
                        atoms.set_distance(Hnear2, i2, 1.4480, fix=1, mask=mask)
                    elif state[0] == 2 or state[0] == 6:
                        atoms.set_distance(Hnear2, i2, 1.3400, fix=1, mask=mask)
                    atoms.set_distance(Hnear2, total_atoms-3, 1.4480, fix=1, mask=mask)
                mask[Hnear2] = 0
                total_atoms = len(atoms)
                CnearH = near(Hnear2, atoms, 'H', limit=False)
                del atoms[CnearH]
                Hnear1 = numjudge(Hnear1, CnearH)
                i1 = numjudge(i1, CnearH)
                i2 = numjudge(i2, CnearH)
                CnearH = near(Hnear2, atoms, 'H', limit=False)
                mask[CnearH] = 1
                for k in range(10):
                    atoms.set_angle(i2, Hnear2, CnearH, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-3, Hnear2, CnearH, 109.5, mask=mask)
                atoms.set_dihedral(i1, i2, Hnear2, CnearH, -150, mask=mask)
                atoms.set_distance(CnearH, Hnear2, 1.0100, fix=1)
            elif needelements[i] == 'O':
                for k in range(10):
                    if state[0] == 1:
                        atoms.set_distance(Hnear2, i2, 1.4100, fix=1, mask=mask)
                    elif state[0] == 2 or state[0] == 6:
                        atoms.set_distance(Hnear2, i2, 1.3640, fix=1, mask=mask)
                    atoms.set_distance(Hnear2, total_atoms-3, 1.4100, fix=1, mask=mask)
                for j in range(2):
                    total_atoms = len(atoms)
                    CnearH = near(Hnear2, atoms, 'H', limit=False)
                    del atoms[CnearH]
                    Hnear2 = numjudge(Hnear2, CnearH)
    return atoms

def bind6(atoms, bond, state=[1, 1], strict={1:'X'}, side=0, shuffle=0, less5=False):
    i1 = bond[0]
    i2 = bond[1]
    if side == 0:
        Hnear1 = near(i1, atoms, 'H')
        Hnear2 = near(i2, atoms, 'H')
    elif side == 1:
        Hnear1 = near(i1, atoms, 'H', reverse=True)
        Hnear2 = near(i2, atoms, 'H', reverse=True)
    elif side == 2:
        Hnear1 = near(i1, atoms, 'H')
        Hnear2 = near(i2, atoms, 'H', reverse=True)
    elif side == 3:
        Hnear1 = near(i1, atoms, 'H', reverse=True)
        Hnear2 = near(i2, atoms, 'H')
    print(Hnear1, Hnear2)
    if less5 == True and side >= 2:
        total_atoms = len(atoms)
        mask = [0]*total_atoms
        mask[Hnear1] = 1
        mask[Hnear2] = 1
        atoms.set_angle(i2, i1, Hnear1, 107.8, mask=mask)
        mask[Hnear1] = 0
        atoms.set_angle(i1, i2, Hnear2, 107.8, mask=mask)
        atoms.set_dihedral(Hnear1, i1, i2, Hnear2, 55, mask=mask)
        xoy = near(Hnear2, atoms, limit=False)
        if atoms.get_distance(xoy, Hnear2) <= 0.9:
            atoms.set_dihedral(Hnear1, i1, i2, Hnear2, -55, mask=mask)
        xoy = near(i2, atoms, target='Not H')
        if xoy == i1:
            xoy = near(i2, atoms, target='Not H', reverse=True)
        if atoms.get_angle(xoy, i2, Hnear2) <= 90:
            atoms.set_dihedral(Hnear1, i1, i2, Hnear2, -55, mask=mask)
    if shuffle == 0:
        atoms = CH3(atoms, Hnear1, state[0])
        atoms = CH3(atoms, Hnear2, state[1])
    elif shuffle == 1:
        atoms = CH3(atoms, Hnear1, state[0], ring=True)
        atoms = CH3(atoms, Hnear2, state[1], ring=True)
    elif shuffle == 2:
        atoms = CH3(atoms, Hnear1, state[0], ring=True)
        atoms = CH3(atoms, Hnear2, state[1])
    elif shuffle == 3:
        atoms = CH3(atoms, Hnear1, state[0])
        atoms = CH3(atoms, Hnear2, state[1], ring=True)
    atoms.append('C')
    atoms.append('C')
    for j in range(4):
        atoms.append('H')
    x = atoms.get_positions()
    x[-5][0] = x[-5][0] + 50
    x[-5][1] = x[-5][1] + 50
    x[-5][2] = x[-5][2] + 50
    x[-6][0] = x[-6][0] - 50
    x[-6][1] = x[-6][1] - 50
    x[-6][2] = x[-6][2] - 50
    for i in range(-2, 0):
        x[2*i][0] = x[2*i][0] - 50
        x[2*i][1] = x[2*i][1] - 50
        x[2*i][2] = x[2*i][2] - 50
        x[2*i+1][0] = x[2*i+1][0] + 50
        x[2*i+1][1] = x[2*i+1][1] + 50
        x[2*i+1][2] = x[2*i+1][2] + 50
    atoms.set_positions(x)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-6) + [1]*6
    atoms.set_angle(i1, Hnear1, total_atoms-6, 107.8, mask=mask)
    if state[0] == 1:
        atoms.set_distance(total_atoms-6, Hnear1, 1.5290, fix=1, mask=mask)
    elif state[0] == 2:
        atoms.set_distance(total_atoms-6, Hnear1, 1.5095, fix=1, mask=mask)
    elif state[0] == 6:
        atoms.set_distance(total_atoms-6, Hnear1, 1.5105, fix=1, mask=mask)
    if less5 == False:
        atoms.set_dihedral(i2, i1, Hnear1, total_atoms-6, 47.9, mask=mask)
        if atoms.get_distance(total_atoms-6, Hnear2) >= 2.7:
            atoms.set_dihedral(i2, i1, Hnear1, total_atoms-6, -47.9, mask=mask)
    else:
        atoms.set_dihedral(i2, i1, Hnear1, total_atoms-6, -47.9, mask=mask)
        if atoms.get_distance(total_atoms-6, Hnear2) >= 2.7:
            atoms.set_dihedral(i2, i1, Hnear1, total_atoms-6, 47.9, mask=mask)
            
    mask[-6] = 0
    atoms.set_angle(i2, Hnear2, total_atoms-5, 107.8, mask=mask)
    for k in range(10):
        #print(atoms.get_dihedral(i1, i2, Hnear2, total_atoms-5))
        if (less5 == True and atoms.get_dihedral(Hnear2, i2, i1, Hnear1) > 10 and
            atoms.get_dihedral(Hnear2, i2, i1, Hnear1) < 350):
            if (atoms.get_dihedral(i1, i2, Hnear2, total_atoms-5) <= 30
                and atoms.get_dihedral(i1, i2, Hnear2, total_atoms-5) >= 0):
                atoms.set_dihedral(i1, i2, Hnear2, total_atoms-5, 47.9, mask=mask)
            if (atoms.get_dihedral(i1, i2, Hnear2, total_atoms-5) >= 330
                and atoms.get_dihedral(i1, i2, Hnear2, total_atoms-5) <= 360):
                atoms.set_dihedral(i1, i2, Hnear2, total_atoms-5, -47.9, mask=mask)
        else:
            if (atoms.get_dihedral(i1, i2, Hnear2, total_atoms-5) <= 30 and
                atoms.get_dihedral(i1, i2, Hnear2, total_atoms-5) >= 0):
                atoms.set_dihedral(i1, i2, Hnear2, total_atoms-5, -47.9, mask=mask)
            if (atoms.get_dihedral(i1, i2, Hnear2, total_atoms-5) >= 330 and
                atoms.get_dihedral(i1, i2, Hnear2, total_atoms-5) <= 360):
                atoms.set_dihedral(i1, i2, Hnear2, total_atoms-5, 47.9, mask=mask)
        atoms.set_distance(total_atoms-5, total_atoms-6, 1.5290, fix=1, mask=mask)
        if state[1] == 1:
            atoms.set_distance(total_atoms-5, Hnear2, 1.5290, fix=1, mask=mask)
        elif state[1] == 2:
            atoms.set_distance(total_atoms-5, Hnear2, 1.5095, fix=1, mask=mask)
        elif state[1] == 6:
            atoms.set_distance(total_atoms-5, Hnear2, 1.5105, fix=1, mask=mask)
    mask[-5] = 0
    for j in range(2):
        total_atoms = len(atoms)
        CnearH = near(total_atoms-6+j, atoms, 'H', limit=False)
        del atoms[CnearH]
        i1 = numjudge(i1, CnearH)
        i2 = numjudge(i2, CnearH)
        Hnear1 = numjudge(Hnear1, CnearH)
        Hnear2 = numjudge(Hnear2, CnearH)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-4)+[1]*4
    atoms.set_distance(total_atoms-4, total_atoms-6, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-4, total_atoms-6, Hnear1, i1, -110, mask=mask)
    for k in range(10):
        atoms.set_angle(Hnear1, total_atoms-6, total_atoms-4, 111.5, mask=mask)
        atoms.set_angle(total_atoms-5, total_atoms-6,\
                        total_atoms-4, 111.5, mask=mask)
    mask[-4] = 0
    atoms.set_distance(total_atoms-3, total_atoms-4, 0, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-3, total_atoms-6, Hnear1, i1,\
                       atoms.get_dihedral(total_atoms-3,\
                                          total_atoms-6, Hnear1, i1)-110, mask=mask)
    xoy = near(total_atoms-3, atoms, limit=False)
    if atoms.get_distance(xoy, total_atoms-3) <= 0.9:
        atoms.set_dihedral(total_atoms-3, total_atoms-6, Hnear1, i1,\
                           atoms.get_dihedral(total_atoms-3,\
                                              total_atoms-6, Hnear1, i1)+220, mask=mask)
    mask[-3] = 0
    atoms.set_distance(total_atoms-2, total_atoms-5, 0, fix=1, mask=mask)
    atoms.set_distance(total_atoms-2, total_atoms-5, 1.0900, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-2, total_atoms-5, Hnear2, i2, 110, mask=mask)
    for k in range(10):
        atoms.set_angle(Hnear2, total_atoms-5, total_atoms-2, 111.5, mask=mask)
        atoms.set_angle(total_atoms-6, total_atoms-5,\
                        total_atoms-2, 111.5, mask=mask)
    mask[-2] = 0
    atoms.set_distance(total_atoms-1, total_atoms-2, 0, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-1, total_atoms-5, Hnear2, i2,\
                       atoms.get_dihedral(total_atoms-1,\
                                          total_atoms-5, Hnear2, i2)-110, mask=mask)
    xoy = near(total_atoms-1, atoms, limit=False)
    if atoms.get_distance(xoy, total_atoms-1) <= 0.9:
        atoms.set_dihedral(total_atoms-1, total_atoms-5, Hnear2, i2,\
                           atoms.get_dihedral(total_atoms-1,\
                                              total_atoms-5, Hnear2, i2)+220, mask=mask)

#    for k in range(10):
#        atoms.set_angle(Hnear2, total_atoms-5, total_atoms-1, 111.5, mask=mask)
#        atoms.set_angle(total_atoms-6, total_atoms-5, total_atoms-1, 111.5, mask=mask)
        
    mask[-1] = 0
    needelements = []
    places = []
    for j in range(4):
        try:
            needelements.append(strict[j+2])
            places.append(j+2)
        except:
            pass
    error = 0
    for i in range(len(needelements)):
        total_atoms = len(atoms)
        mask = [0]*total_atoms
        if places[i] == 2:
            mask[Hnear1] = 1
            elements = atoms.get_chemical_symbols()
            elements[Hnear1] = needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i] == 'N':
                for k in range(10):
                    if state[0] == 1:
                        atoms.set_distance(Hnear1, i1, 1.4480, fix=1, mask=mask)
                    elif state[0] == 2 or state[0] == 6:
                        atoms.set_distance(Hnear1, i1, 1.3400, fix=1, mask=mask)
                    atoms.set_distance(Hnear1, total_atoms-6, 1.4480, fix=1, mask=mask)
                mask[Hnear1] = 0
                total_atoms = len(atoms)
                CnearH = near(Hnear1, atoms, 'H', limit=False)
                del atoms[CnearH]
                Hnear1 = numjudge(Hnear1, CnearH)
                i1 = numjudge(i1, CnearH)
                i2 = numjudge(i2, CnearH)
                CnearH = near(Hnear1, atoms, 'H', limit=False)
                mask[CnearH] = 1
                for k in range(10):
                    atoms.set_angle(i1, Hnear1, CnearH, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-6, Hnear1, CnearH, 109.5, mask=mask)
                atoms.set_dihedral(i2, i1, Hnear1, CnearH, -150, mask=mask)
                atoms.set_distance(CnearH, Hnear1, 1.0100, fix=1)
            elif needelements[i] == 'O':
                for k in range(10):
                    if state[0] == 1:
                        atoms.set_distance(Hnear1, i1, 1.4100, fix=1, mask=mask)
                    elif state[0] == 2 or state[0] == 6:
                        atoms.set_distance(Hnear1, i1, 1.3640, fix=1, mask=mask)
                    atoms.set_distance(Hnear1, total_atoms-6, 1.4100, fix=1, mask=mask)
                for j in range(2):
                    total_atoms = len(atoms)
                    CnearH = near(Hnear1, atoms, 'H', limit=False)
                    del atoms[CnearH]
                    Hnear1 = numjudge(Hnear1, CnearH)
        if places[i] == 3:
            mask[total_atoms-6+error] = 1
            elements = atoms.get_chemical_symbols()
            elements[total_atoms-6+error] = needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i] == 'N':
                for k in range(10):
                    atoms.set_distance(total_atoms-6+error, Hnear1, 1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms-6+error, total_atoms-5+error,\
                                       1.4480, fix=1, mask=mask)
                mask[total_atoms-6+error] = 0
                total_atoms = len(atoms)
                CnearH = near(total_atoms-6+error, atoms, 'H', limit=False)
                del atoms[CnearH]
                Hnear1 = numjudge(Hnear1, CnearH)
                i1 = numjudge(i1, CnearH)
                i2 = numjudge(i2, CnearH)
                CnearH = near(total_atoms-6+error, atoms, 'H', limit=False)
                mask[CnearH] = 1
                atoms.set_dihedral(i1, Hnear1, total_atoms-6+error, CnearH, -150, mask=mask)
                for k in range(10):
                    atoms.set_angle(Hnear1, total_atoms-6+error, CnearH, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-5+error, total_atoms-6+error, CnearH,\
                                    109.5, mask=mask)
                atoms.set_distance(CnearH, total_atoms-6+error, 1.0100, fix=1)
                error = error+1
            elif needelements[i] == 'O':
                for k in range(10):
                    atoms.set_distance(total_atoms-6+error, Hnear1, 1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms-6+error, total_atoms-5-error, 1.4100,\
                                       fix=1, mask=mask)
                for j in range(2):
                    total_atoms = len(atoms)
                    CnearH = near(total_atoms-6+error+j, atoms, 'H', limit=False)
                    del atoms[CnearH]
                error = error+2
        if places[i] == 4:
            mask[total_atoms-5+error] = 1
            elements = atoms.get_chemical_symbols()
            elements[total_atoms-5-error] = needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i] == 'N':
                for k in range(10):
                    atoms.set_distance(total_atoms-5+error, total_atoms-6+error,\
                                       1.4480, fix=1, mask=mask)
                    atoms.set_distance(total_atoms-5+error, Hnear2, 1.4480, fix=1, mask=mask)
                mask[total_atoms-5+error] = 0
                total_atoms = len(atoms)
                CnearH = near(total_atoms-5+error, atoms, 'H', limit=False)
                del atoms[CnearH]
                Hnear1 = numjudge(Hnear1, CnearH)
                i1 = numjudge(i1, CnearH)
                i2 = numjudge(i2, CnearH)
                CnearH = near(total_atoms-5+error, atoms, 'H', limit=False)
                mask[CnearH] = 1
                atoms.set_dihedral(i2, Hnear2, total_atoms-5+error, CnearH, -150, mask=mask)
                for k in range(10):
                    atoms.set_angle(Hnear2, total_atoms-5+error, CnearH, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-6+error, total_atoms-5+error,\
                                    CnearH, 109.5, mask=mask)
                atoms.set_distance(CnearH, total_atoms-5+error, 1.0100, fix=1)
                error = error+1
            elif needelements[i] == 'O':
                for k in range(10):
                    atoms.set_distance(total_atoms-5+error, total_atoms-6+error,\
                                       1.4100, fix=1, mask=mask)
                    atoms.set_distance(total_atoms-5+error, Hnear2, 1.4100, fix=1, mask=mask)
                for j in range(2):
                    total_atoms = len(atoms)
                    CnearH = near(total_atoms-5-error+j, atoms, 'H', limit=False)
                    del atoms[CnearH]
                error = error+2
        if places[i] == 5:
            mask[Hnear2] = 1
            elements = atoms.get_chemical_symbols()
            elements[Hnear2] = needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i] == 'N':
                for k in range(10):
                    if state[0] == 1:
                        atoms.set_distance(Hnear2, i2, 1.4480, fix=1, mask=mask)
                    elif state[0] == 2 or state[0] == 6:
                        atoms.set_distance(Hnear2, i2, 1.3400, fix=1, mask=mask)
                    atoms.set_distance(Hnear2, total_atoms-5+error, 1.4480, fix=1, mask=mask)
                mask[Hnear2] = 0
                total_atoms = len(atoms)
                CnearH = near(Hnear2, atoms, 'H', limit=False)
                del atoms[CnearH]
                Hnear1 = numjudge(Hnear1, CnearH)
                i1 = numjudge(i1, CnearH)
                i2 = numjudge(i2, CnearH)
                CnearH = near(Hnear2, atoms, 'H', limit=False)
                mask[CnearH] = 1
                for k in range(10):
                    atoms.set_angle(i2, Hnear2, CnearH, 109.5, mask=mask)
                    atoms.set_angle(total_atoms-5+error, Hnear2, CnearH, 109.5, mask=mask)
                atoms.set_dihedral(i1, i2, Hnear2, CnearH, -150, mask=mask)
                atoms.set_distance(CnearH, Hnear2, 1.0100, fix=1)
            elif needelements[i] == 'O':
                for k in range(10):
                    if state[0] == 1:
                        atoms.set_distance(Hnear2, i2, 1.4100, fix=1, mask=mask)
                    elif state[0] == 2 or state[0] == 6:
                        atoms.set_distance(Hnear2, i2, 1.3640, fix=1, mask=mask)
                    atoms.set_distance(Hnear2, total_atoms-5+error, 1.4100, fix=1, mask=mask)
                for j in range(2):
                    total_atoms = len(atoms)
                    CnearH = near(Hnear2, atoms, 'H', limit=False)
                    del atoms[CnearH]
                    Hnear2 = numjudge(Hnear2, CnearH)
    return atoms

def bindconj5(atoms, bond, state=[1, 1], strict={1:'X'}):
    try:
        if strict[1] == 'X':
            strict[3] = 'N'
    except:
        pass
    i1 = bond[0]
    i2 = bond[1]
    Cnear1 = near(i1, atoms, 'Not H')
    Cnear2 = near(i2, atoms, 'Not H', reverse=True)
    total_atoms = len(atoms)
    mask = [0] * total_atoms
    mask[i1] = 1
    mask[i2] = 1
    for j in range(3):
        try:
            element = strict[j+2]
            place = j + 2
        except:
            pass
    bondjudge5([Cnear1, Cnear2], atoms, [i1, i2], state, mask, place)
    for j in range(2):
        CnearH = near(i1, atoms, 'H')
        del atoms[CnearH]
        i1 = numjudge(i1, CnearH)
        i2 = numjudge(i2, CnearH)
    for j in range(2):
        CnearH = near(i2, atoms, 'H')
        del atoms[CnearH]
        i1 = numjudge(i1, CnearH)
        i2 = numjudge(i2, CnearH)
    if place == 2:
        atoms.append(element)
        atoms.append('C')
        atoms.append('C')
    elif place == 3:
        atoms.append('C')
        atoms.append('C')
        atoms.append(element)
    elif place == 4:
        atoms.append('C')
        atoms.append(element)
        atoms.append('C')
    for j in range(3):
        atoms.append('H')
    x = atoms.get_positions()
    x[-6][0] = x[-6][0] - 100
    x[-6][1] = x[-6][1] - 50
    x[-6][2] = x[-6][2] - 50
    x[-5][0] = x[-5][0] + 50
    x[-5][1] = x[-5][1] + 50
    x[-5][2] = x[-5][2] + 50
    for j in range(-4, 0):
        x[j][0] = x[j][0] + 50
        x[j][1] = x[j][1] + 50
        x[j][2] = x[j][2] + 50
    atoms.set_positions(x)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-6)+[1]*6
    if place == 2:
        if element == 'N':
            atoms.set_angle(i2, i1, total_atoms-6, 107.7, mask=mask)
            if state[0] == 1:
                atoms.set_distance(total_atoms-6, i1, 1.3810, fix=1, mask=mask)
            elif state[0] == 2:
                atoms.set_distance(total_atoms-6, i1, 1.3740, fix=1, mask=mask)
        elif element == 'O':
            atoms.set_angle(i2, i1, total_atoms-6, 110.6, mask=mask)
            atoms.set_distance(total_atoms-6, i1, 1.3600, fix=1, mask=mask)
        mask[-6] = 0
        atoms.set_angle(i1, i2, total_atoms-5, 107.3, mask=mask)
        atoms.set_distance(total_atoms-5, i2, 1.4240, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-6, i1, i2, total_atoms-5, 0.5, mask=mask)
        mask[-5] = 0
        if element == 'N':
            for k in range(10):
                atoms.set_distance(total_atoms-4, total_atoms-6, 1.3810, fix=1, mask=mask)
                atoms.set_distance(total_atoms-4, total_atoms-5, 1.3670, fix=1, mask=mask)
            mask[-4] = 0
            #atoms.set_distance(total_atoms-1, total_atoms-6, 0, fix=1, mask=mask)
            atoms.set_dihedral(total_atoms-5, total_atoms-4, total_atoms-6,\
                               total_atoms-1, 0.5, mask=mask)
            for k in range(10):
                atoms.set_angle(total_atoms-4, total_atoms-6, total_atoms-1, 124.5, mask=mask)
                atoms.set_angle(i1, total_atoms-6, total_atoms-1, 124.5, mask=mask)
            atoms.set_distance(total_atoms-1, total_atoms-6, 1.0100, fix=1, mask=mask)
            mask[-1] = 0
        elif element == 'O':
            for k in range(10):
                atoms.set_distance(total_atoms-4, total_atoms-6, 1.3600, fix=1, mask=mask)
                atoms.set_distance(total_atoms-4, total_atoms-5, 1.3670, fix=1, mask=mask)
            mask[-4] = 0
        atoms.set_distance(total_atoms-3, total_atoms-4, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-3, total_atoms-4, total_atoms-5,\
                           i2, 0.5, mask=mask)
        for k in range(10):
            atoms.set_angle(total_atoms-6, total_atoms-4, total_atoms-3, 126.5, mask=mask)
            atoms.set_angle(total_atoms-5, total_atoms-4, total_atoms-3, 126.5, mask=mask)
        mask[-3] = 0
        atoms.set_distance(total_atoms-2, total_atoms-5, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-2, total_atoms-5, i2, i1, 0.5, mask=mask)
        for k in range(10):
            atoms.set_angle(i2, total_atoms-5, total_atoms-2, 126.5, mask=mask)
            atoms.set_angle(total_atoms-4, total_atoms-5, total_atoms-2, 126.5, mask=mask)
        mask[-2] = 0
    elif place == 3:
        atoms.set_angle(i2, i1, total_atoms-6, 107.3, mask=mask)
        atoms.set_distance(total_atoms-6, i1, 1.3670, fix=1, mask=mask)
        mask[-6] = 0
        atoms.set_angle(i1, i2, total_atoms-5, 107.3, mask=mask)
        atoms.set_distance(total_atoms-5, i2, 1.3670, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-6, i1, i2, total_atoms-5, 0.5, mask=mask)
        mask[-5] = 0
        if element == 'N':
            for k in range(10):
                atoms.set_distance(total_atoms-4, total_atoms-6, 1.3810, fix=1, mask=mask)
                atoms.set_distance(total_atoms-4, total_atoms-5, 1.3810, fix=1, mask=mask)
            mask[-4] = 0
            atoms.set_distance(total_atoms-1, total_atoms-4, 0, fix=1, mask=mask)
            atoms.set_distance(total_atoms-1, total_atoms-4, 1.0100, fix=1, mask=mask)
            atoms.set_dihedral(i1, total_atoms-6, total_atoms-4,\
                               total_atoms-1, 110, mask=mask)
            for k in range(10):
                atoms.set_angle(total_atoms-6, total_atoms-4, total_atoms-1, 124.5, mask=mask)
                atoms.set_angle(total_atoms-5, total_atoms-4, total_atoms-1, 124.5, mask=mask)
            mask[-1] = 0
        elif element == 'O':
            for k in range(10):
                atoms.set_distance(total_atoms-4, total_atoms-6, 1.3600, fix=1, mask=mask)
                atoms.set_distance(total_atoms-4, total_atoms-5, 1.3600, fix=1, mask=mask)
            mask[-4] = 0
        atoms.set_distance(total_atoms-3, total_atoms-6, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-3, total_atoms-6, i1, i2, 0.5, mask=mask)
        for k in range(10):
            atoms.set_angle(i1, total_atoms-6, total_atoms-3, 126.5, mask=mask)
            atoms.set_angle(total_atoms-4, total_atoms-6, total_atoms-3, 126.5, mask=mask)
        mask[-3] = 0
        atoms.set_distance(total_atoms-2, total_atoms-5, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-2, total_atoms-5, i2, i1, 0.5, mask=mask)
        for k in range(10):
            atoms.set_angle(i2, total_atoms-5, total_atoms-2, 126.5, mask=mask)
            atoms.set_angle(total_atoms-4, total_atoms-5, total_atoms-2, 126.5, mask=mask)
        mask[-2] = 0
    elif place == 4:
        atoms.set_angle(i2, i1, total_atoms-6, 107.3, mask=mask)
        atoms.set_distance(total_atoms-6, i1, 1.4240, fix=1, mask=mask)
        mask[-6] = 0
        if element == 'N':
            atoms.set_angle(i1, i2, total_atoms-5, 107.7, mask=mask)
            if state[0] == 1:
                atoms.set_distance(total_atoms-5, i2, 1.3810, fix=1, mask=mask)
            if state[0] == 2:
                atoms.set_distance(total_atoms-5, i2, 1.3740, fix=1, mask=mask)
        elif element == 'O':
            atoms.set_angle(i1, i2, total_atoms-5, 110.6, mask=mask)
            atoms.set_distance(total_atoms-5, i2, 1.3600, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-6, i1, i2, total_atoms-5, 0.5, mask=mask)
        mask[-5] = 0
        if element == 'N':
            for k in range(10):
                atoms.set_distance(total_atoms-4, total_atoms-5, 1.3810, fix=1, mask=mask)
                atoms.set_distance(total_atoms-4, total_atoms-6, 1.3670, fix=1, mask=mask)
            mask[-4] = 0
            #atoms.set_distance(total_atoms-1, total_atoms-6, 0, fix=1, mask=mask)
            atoms.set_dihedral(total_atoms-6, total_atoms-4, total_atoms-5,\
                               total_atoms-1, 0.5, mask=mask)
            for k in range(10):
                atoms.set_angle(total_atoms-4, total_atoms-5, total_atoms-1, 124.5, mask=mask)
                atoms.set_angle(i2, total_atoms-5, total_atoms-1, 124.5, mask=mask)
            atoms.set_distance(total_atoms-1, total_atoms-5, 1.0100, fix=1, mask=mask)
            mask[-1] = 0
        elif element == 'O':
            for k in range(10):
                atoms.set_distance(total_atoms-4, total_atoms-5, 1.3600, fix=1, mask=mask)
                atoms.set_distance(total_atoms-4, total_atoms-6, 1.3670, fix=1, mask=mask)
            mask[-4] = 0
        atoms.set_distance(total_atoms-3, total_atoms-4, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-3, total_atoms-4,\
                           total_atoms-5, i2, 0.5, mask=mask)
        for k in range(10):
            atoms.set_angle(total_atoms-6, total_atoms-4, total_atoms-3, 126.5, mask=mask)
            atoms.set_angle(total_atoms-5, total_atoms-4, total_atoms-3, 126.5, mask=mask)
        mask[-3] = 0
        atoms.set_distance(total_atoms-2, total_atoms-6, 1.0100, fix=1, mask=mask)
        atoms.set_dihedral(total_atoms-2, total_atoms-6,\
                           total_atoms-4, total_atoms-5, 0.5, mask=mask)
        for k in range(10):
            atoms.set_angle(i1, total_atoms-6, total_atoms-2, 126.5, mask=mask)
            atoms.set_angle(total_atoms-4, total_atoms-6, total_atoms-2, 126.5, mask=mask)
        mask[-2] = 0
    if element != 'N':
        del atoms[-1]
    return atoms

def bindconj6(atoms, bond, state=[1, 1], strict={1:'X'}):
    i1 = bond[0]
    i2 = bond[1]
    Cnear1 = near(i1, atoms, 'Not H')
    Cnear2 = near(i2, atoms, 'Not H', reverse=True)
    total_atoms = len(atoms)
    mask = [0]*total_atoms
    mask[i1] = 1
    mask[i2] = 1
    for j in range(4):
        try:
            element = strict[j+2]
            place = j+2
            break
        except:
            place = 1
    bondjudge6([Cnear1, Cnear2], atoms, [i1, i2], state, mask, place)
    for j in range(2):
        CnearH = near(i1, atoms, 'H')
        del atoms[CnearH]
        i1 = numjudge(i1, CnearH)
        i2 = numjudge(i2, CnearH)
    for j in range(2):
        CnearH = near(i2, atoms, 'H')
        del atoms[CnearH]
        i1 = numjudge(i1, CnearH)
        i2 = numjudge(i2, CnearH)    
    if place == 2:
        atoms.append(element)
        atoms.append('C')
        atoms.append('C')
        atoms.append('C')
    elif place == 3:
        atoms.append('C')
        atoms.append('C')
        atoms.append(element)
        atoms.append('C')
    elif place == 4:
        atoms.append('C')
        atoms.append('C')
        atoms.append('C')
        atoms.append(element)
    elif place == 5:
        atoms.append('C')
        atoms.append(element)
        atoms.append('C')
        atoms.append('C')
    else:
        for j in range(4):
            atoms.append('C')
    for j in range(4):
        atoms.append('H')
    x = atoms.get_positions()
    x[-8][0] = x[-8][0] - 100
    x[-8][1] = x[-8][1] - 50
    x[-8][2] = x[-8][2] - 50
    x[-7][0] = x[-7][0] + 50
    x[-7][1] = x[-7][1] + 50
    x[-7][2] = x[-7][2] + 50
    x[-6][0] = x[-6][0] - 50
    x[-6][1] = x[-6][1] - 50
    x[-6][2] = x[-6][2] - 50
    x[-5][0] = x[-5][0] + 50
    x[-5][1] = x[-5][1] + 50
    x[-5][2] = x[-5][2] + 50
    for j in range(-2, 0):
        x[2*j][0] = x[2*j][0] + 50
        x[2*j][1] = x[2*j][1] + 50
        x[2*j][2] = x[2*j][2] + 50
        x[2*j+1][0] = x[2*j+1][0] - 50
        x[2*j+1][1] = x[2*j+1][1] - 50
        x[2*j+1][2] = x[2*j+1][2] - 50
    atoms.set_positions(x)
    total_atoms = len(atoms)
    mask = [0]*(total_atoms-8) + [1]*8
    atoms.set_angle(i2, i1, total_atoms-8, 120, mask=mask)
    if state[0] == 2:
        atoms.set_distance(total_atoms-8, i1, 1.4040, fix=1, mask=mask)
    elif state[0] == 82:
        atoms.set_distance(total_atoms-8, i1, 1.4330, fix=1, mask=mask)
    else:
        atoms.set_distance(total_atoms-8, i1, 1.4000, fix=1, mask=mask)
    mask[-8] = 0
    atoms.set_angle(i1, i2, total_atoms-7, 120, mask=mask)
    if state[0] == 2:
        atoms.set_distance(total_atoms-7, i2, 1.4040, fix=1, mask=mask)
    elif state[0] == 82:
        atoms.set_distance(total_atoms-7, i2, 1.4330, fix=1, mask=mask)
    else:
        atoms.set_distance(total_atoms-7, i2, 1.4000, fix=1, mask=mask)
    if state[0] == 2 or state[0] == 82:
        atoms.set_dihedral(total_atoms-8, i1, i2, total_atoms-7, 5, mask=mask)
    else:
        atoms.set_dihedral(total_atoms-8, i1, i2, total_atoms-7, 0.5, mask=mask)
    mask[-7] = 0
    atoms.set_angle(i1, total_atoms-8, total_atoms-6, 120, mask=mask)
    if state[0] == 82:
        atoms.set_distance(total_atoms-6, total_atoms-8, 1.4330, fix=1, mask=mask)
    else:
        atoms.set_distance(total_atoms-6, total_atoms-8, 1.4000, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-6, total_atoms-8, i1, i2, -0.5, mask=mask)
    mask[-6] = 0
    atoms.set_angle(i2, total_atoms-7, total_atoms-5, 120, mask=mask)
    if state[0] == 82:
        atoms.set_distance(total_atoms-5, total_atoms-7, 1.4330, fix=1, mask=mask)
    else:
        atoms.set_distance(total_atoms-5, total_atoms-7, 1.4000, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-5, total_atoms-7, i2, i1, -0.5, mask=mask)
    mask[-5] = 0
    atoms.set_angle(total_atoms-4, total_atoms-8, i1, 120, mask=mask)
    atoms.set_distance(total_atoms-4, total_atoms-8, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(i2, i1, total_atoms-8, total_atoms-4, -179.5, mask=mask)
    mask[-4] = 0
    atoms.set_angle(total_atoms-3, total_atoms-7, i2, 120, mask=mask)
    atoms.set_distance(total_atoms-3, total_atoms-7, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(i1, i2, total_atoms-7, total_atoms-3, -179.5, mask=mask)
    mask[-3] = 0
    atoms.set_angle(total_atoms-8, total_atoms-6, total_atoms-2, 120, mask=mask)
    atoms.set_dihedral(total_atoms-2, total_atoms-6,\
                       total_atoms-8, i1, 179.5, mask=mask)
    atoms.set_distance(total_atoms-2, total_atoms-6, 1.0800, fix=1, mask=mask)
    mask[-2] = 0
    atoms.set_angle(total_atoms-7, total_atoms-5, total_atoms-1, 120, mask=mask)
    atoms.set_distance(total_atoms-1, total_atoms-5, 1.0800, fix=1, mask=mask)
    atoms.set_dihedral(total_atoms-8, total_atoms-6,\
                       total_atoms-7, total_atoms-1, 179.5, mask=mask)
    if place == 2:
        elements = atoms.get_chemical_symbols()
        elements[-8] = 'N'
        atoms.set_chemical_symbols(elements)
        for k in range(10):
            atoms.set_distance(total_atoms-8, total_atoms-6, 1.3390, fix=1)
            if state[0] == 2:
                atoms.set_distance(total_atoms-8, i1, 1.3540, fix=1)
            else:
                atoms.set_distance(total_atoms-8, i1, 1.3390, fix=1)
        CnearH = near(total_atoms-8, atoms, 'H')
        del atoms[CnearH]
    elif place == 3:
        elements = atoms.get_chemical_symbols()
        elements[-6] = 'N'
        atoms.set_chemical_symbols(elements)
        for k in range(10):
            atoms.set_distance(total_atoms-6, total_atoms-8, 1.3390, fix=1)
            atoms.set_distance(total_atoms-6, total_atoms-5, 1.3390, fix=1)
        CnearH = near(total_atoms-6, atoms, 'H')
        del atoms[CnearH]
    elif place == 4:
        elements = atoms.get_chemical_symbols()
        elements[-5] = 'N'
        atoms.set_chemical_symbols(elements)
        for k in range(10):
            atoms.set_distance(total_atoms-5, total_atoms-6, 1.3390, fix=1)
            atoms.set_distance(total_atoms-5, total_atoms-7, 1.3390, fix=1)
        CnearH = near(total_atoms-5, atoms, 'H')
        del atoms[CnearH]
    elif place == 5:
        elements = atoms.get_chemical_symbols()
        elements[-7] = 'N'
        atoms.set_chemical_symbols(elements)
        for k in range(10):
            if state[0] == 2:
                atoms.set_distance(total_atoms-7, i2, 1.3540, fix=1)
            else:
                atoms.set_distance(total_atoms-7, i2, 1.3390, fix=1)
            atoms.set_distance(total_atoms-7, total_atoms-5, 1.3390, fix=1)
        CnearH = near(total_atoms-7, atoms, 'H')
        del atoms[CnearH]
    return atoms
        
def numjudge(i, deleted):
    if deleted<i:
        return i-1
    else:
        return i

def bondjudge2(Cnears, atoms, iis, elements, state, mask, times=10):
    print(Cnears, iis, elements)
    if 'N' in elements:
        for k in range(times):
            if state[0] == state[1] and state[0] == 2:
                atoms.set_distance(iis[0], Cnears[0], 1.5095, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.3400, fix=1, mask=mask)
                atoms.set_distance(iis[0], iis[1], 1.2900, fix=1, mask=mask)
        return 7
    for k in range(times):
        if state[0] == state[1] and state[0] == 2:
            atoms.set_distance(iis[0], Cnears[0], 1.4000, fix=1, mask=mask)
            atoms.set_distance(iis[1], Cnears[1], 1.4000, fix=1, mask=mask)
            atoms.set_distance(iis[0], iis[1], 1.4000, fix=1, mask=mask)
        else:
            for i in range(2):
                if state[i] == 1:
                    atoms.set_distance(iis[i], Cnears[i], 1.5095, fix=1, mask=mask)
                elif state[i] == 2:
                    atoms.set_distance(iis[i], Cnears[i], 1.4600, fix=1, mask=mask)
                elif state[i] == 7:
                    atoms.set_distance(iis[i], Cnears[i], 1.3400, fix=1, mask=mask)
                elif state[i] == 8:
                    atoms.set_distance(iis[i], Cnears[i], 1.3640, fix=1, mask=mask)
                elif state[i] == 16:
                    atoms.set_distance(iis[i], Cnears[i], 1.7700, fix=1, mask=mask)
                elif state[i] == 82:
                    atoms.set_distance(iis[i], Cnears[i], 1.4440, fix=1, mask=mask)
            atoms.set_distance(iis[0], iis[1], 1.3400, fix=1, mask=mask)

def bondjudge5(Cnears, atoms, iis, state, mask, place, times=10):
    if place == 2:
        for k in range(times):
            if state[0] == state[1] and state[0] == 1:
                atoms.set_distance(iis[0], Cnears[0], 1.5040, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.4950, fix=1, mask=mask)
                atoms.set_distance(i1, i2, 1.3670, fix=1, mask=mask)
            elif state[0] == state[1] and state[0] == 2:
                atoms.set_distance(iis[0], Cnears[0], 1.4040, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.4000, fix=1, mask=mask)
                atoms.set_distance(i1, i2, 1.4240, fix=1, mask=mask)
    elif place == 3:
        for k in range(times):
            if state[0] == state[1] and state[0] == 2:
                atoms.set_distance(iis[0], Cnears[0], 1.4000, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.4000, fix=1, mask=mask)
            else:
                for i in range(2):
                    if state[i] == 1:
                        atoms.set_distance(iis[i], Cnears[i], 1.4950, fix=1, mask=mask)
                    elif state[i] == 2:
                        atoms.set_distance(iis[i], Cnears[i], 1.4330, fix=1, mask=mask)
                    elif state[i] == 7:
                        atoms.set_distance(iis[i], Cnears[i], 1.4950, fix=1, mask=mask)
                    elif state[i] == 8:
                        atoms.set_distance(iis[i], Cnears[i], 1.3640, fix=1, mask=mask)
                    elif state[i] == 16:
                        atoms.set_distance(iis[i], Cnears[i], 1.7700, fix=1, mask=mask)
                    #elif state[i] == 82:
                    #    atoms.set_distance(iis[i], Cnears[i], 1.4440?, fix=1, mask=mask)
            atoms.set_distance(i1, i2, 1.4240, fix=1, mask=mask)
    elif place == 4:
        for k in range(times):
            if state[0] == state[1] and state[0] == 1:
                atoms.set_distance(iis[0], Cnears[0], 1.4950, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.5040, fix=1, mask=mask)
            elif state[0] == state[1] and state[0] == 2:
                atoms.set_distance(iis[0], Cnears[0], 1.4000, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.4040, fix=1, mask=mask)
            atoms.set_distance(i1, i2, 1.4240, fix=1, mask=mask)

def bondjudge6(Cnears, atoms, iis, state, mask, place, times=10):
    if place == 2 or place == 5:
        for k in range(times):
            if state[0] == state[1] and state[0] == 1:
                atoms.set_distance(iis[0], Cnears[0], 1.5100, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.5100, fix=1, mask=mask)
                atoms.set_distance(i1, i2, 1.4000, fix=1, mask=mask)
            elif state[0] == state[1] and state[0] == 2:
                atoms.set_distance(iis[0], Cnears[0], 1.4040, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.4040, fix=1, mask=mask)
                atoms.set_distance(i1, i2, 1.3700, fix=1, mask=mask)
    elif place == 1 or place == 3 or place == 4:
        for k in range(times):
            if state[0] == state[1] and state[0] == 2:
                atoms.set_distance(iis[0], Cnears[0], 1.4040, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.4040, fix=1, mask=mask)
                atoms.set_distance(i1, i2, 1.3700, fix=1, mask=mask)
            elif state[0] == state[1] and state[0] == 82:
                atoms.set_distance(iis[0], Cnears[0], 1.4190, fix=1, mask=mask)
                atoms.set_distance(iis[1], Cnears[1], 1.4190, fix=1, mask=mask)
                atoms.set_distance(i1, i2, 1.3700, fix=1, mask=mask)
            else:
                for i in range(2):
                    if state[i] == 1:
                        atoms.set_distance(iis[i], Cnears[i], 1.5105, fix=1, mask=mask)
                    elif state[i] == 2:
                        atoms.set_distance(iis[i], Cnears[i], 1.4330, fix=1, mask=mask)
                    elif state[i] == 7:
                        atoms.set_distance(iis[i], Cnears[i], 1.3400, fix=1, mask=mask)
                    elif state[i] == 8:
                        atoms.set_distance(iis[i], Cnears[i], 1.3640, fix=1, mask=mask)
                    elif state[i] == 16:
                        atoms.set_distance(iis[i], Cnears[i], 1.7700, fix=1, mask=mask)
                    elif state[i] == 82:
                        atoms.set_distance(iis[i], Cnears[i], 1.4190, fix=1, mask=mask)
                    if state[i] == 2:
                        atoms.set_distance(i1, i2, 1.4040, fix=1, mask=mask)
                    else:
                        atoms.set_distance(i1, i2, 1.4000, fix=1, mask=mask)


def near(i, atoms, target='ALL', ignore=[], limit=True, reverse=False):
    elements = atoms.get_chemical_symbols()
    total_atoms = len(atoms)
    j = 0
    targeted_distances = []
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
        if total_atoms+_near in ignore:
            continue
        if reverse == True:
            if total_atoms+_near in ignore:
                _near = _near - 1
                continue
            if abs(min(targeted_distances)-atoms.get_distances(i,\
                                                               total_atoms+_near)) <= 0.001:
                #print(_near, atoms.get_distances(i, _near))
                return total_atoms+_near
            _near = _near - 1
        else:
            if _near in ignore:
                _near = _near + 1
                continue
            if abs(min(targeted_distances)-atoms.get_distances(i, _near)) <= 0.001:
                #print(_near, atoms.get_distances(i, _near))
                return _near
            #print(min(targeted_distances), atoms.get_distances(i, _near))
            _near = _near + 1