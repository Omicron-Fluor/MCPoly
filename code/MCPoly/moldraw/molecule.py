import math as m
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.visualize import view
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit import Chem
from rdkit.ML.Cluster import Butina
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolToFile
import py3Dmol
from rdkit.Chem import rdDetermineBonds
import os
import re

def gen_conformers(mol, numConfs=100, maxAttempts=1000, pruneRmsThresh=0.1, useExpTorsionAnglePrefs=True, useBasicKnowledge=True, enforceChirality=True):
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, maxAttempts=maxAttempts, pruneRmsThresh=pruneRmsThresh, useExpTorsionAnglePrefs=useExpTorsionAnglePrefs, useBasicKnowledge=useBasicKnowledge, enforceChirality=enforceChirality, numThreads=0)
    return list(ids)

def calc_energy(mol, conformerId, minimizeIts):
    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conformerId)
    ff.Initialize()
    ff.CalcEnergy()
    results = []
    if minimizeIts > 0:
        results.append(ff.Minimize(maxIts=minimizeIts))
    results.append(ff.CalcEnergy())
    return results

def cluster_conformers(mol, mode="RMSD", threshold=2.0):
    if mode == "TFD":
        dmat = TorsionFingerprints.GetTFDMatrix(mol)
    else:
        dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
    rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True)
    return rms_clusters

def cleanup_qm9_xyz(fname):
    ind = open(fname).readlines()
    nAts = int(ind[0])
    # There are two smiles in the data: the one from GDB and the one assigned from the
    # 3D coordinates in the QM9 paper using OpenBabel (I think).
    gdb_smi,relax_smi = ind[-2].split()[:2]
    ind[1] = '\n'
    ind = ind[:nAts+2]
    for i in range(2,nAts+2):
        l = ind[i]
        l = l.split('\t')
        l.pop(-1)
        ind[i] = '\t'.join(l)+'\n'
    ind = ''.join(ind)
    return ind,gdb_smi,relax_smi

def draw_with_spheres(mol):
    v = py3Dmol.view(width=600,height=400)
    IPythonConsole.addMolToView(mol,v)
    v.zoomTo()
    v.setStyle({'sphere':{'radius':0.4},'stick':{'radius':0.1}});
    v.show()
    return v

def createmol(filename):
    ind = open('{0}.xyz'.format(filename),'r+')
    mains = ind.readlines()
    main=''
    for i,line in enumerate(mains):
        if i==1:
            main=main+'\n'
        else:
            main=main+mains[i]
    raw_mol = Chem.MolFromXYZBlock(main)
    conn_mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineConnectivity(conn_mol)
    print(Chem.MolToMolBlock(conn_mol),file=open('{0}.mol'.format(filename),'w+'))

def sub_selection(substitution,atoms,i,state=1,place=0,strict={1:'X'},shuffle=0,angle=-1):
    '''This function will change the hydrogen into other substitution.''' 
    
    if substitution=='H' or substitution=='F' or substitution=='Cl' or substitution=='Br' or substitution=='I':
        return Halogen(substitution,atoms,i,state)
    elif substitution=='OH':
        return OH(atoms,i,state,angle)
    elif substitution=='OMe' or substitution=='OCH3':
        return OMe(atoms,i,state,angle)
    elif substitution=='SH':
        return SH(atoms,i,state,angle)
    elif substitution=='SMe' or substitution=='SCH3':
        return SMe(atoms,i,state,angle)
    elif substitution=='NH2':
        return NH2(atoms,i,state)
    elif substitution=='NMe2' or substitution=='N(CH3)2':
        return NMe2(atoms,i,state)
    elif substitution=='CH3' or substitution=='Me':
        return CH3(atoms,i,state)
    elif substitution=='C2H3':
        return C2H3(atoms,i,state,shuffle)
    elif substitution=='C2H5' or substitution=='Et' or substitution=='i-C3H7' or substitution=='i-Pr' or substitution=='iPr' or substitution=='t-C4H9' or substitution=='t-Bu' or substitution=='tBu':
        if substitution=='C2H5' or substitution=='Et':
            return CChainR(atoms,i,1,state)
        if substitution=='i-C3H7' or substitution=='i-Pr' or substitution=='iPr':
            return CChainR(atoms,i,2,state)
        if substitution=='t-C4H9' or substitution=='t-Bu' or substitution=='tBu':
            return CChainR(atoms,i,3,state)
    elif substitution=='C3H7' or substitution=='n-C3H7' or substitution=='Pr' or substitution=='n-Pr' or substitution=='nPr' or substitution=='C4H9' or substitution=='n-C4H9' or substitution=='nBu' or substitution=='Bu' or substitution=='n-Bu': 
        if substitution=='C3H7' or substitution=='n-C3H7' or substitution=='n-Pr' or substitution=='nPr' or substitution=='Pr':
            return CChainL(atoms,i,2,state)
        if substitution=='C4H9' or substitution=='n-C4H9' or substitution=='n-Bu' or substitution=='nBu' or substitution=='Bu':
            return CChainL(atoms,i,3,state)
    elif substitution=='CF3' or substitution=='CCl3' or substitution=='CBr3' or substitution=='CI3':
        return trihalo(substitution,atoms,i,state)
    elif substitution=='NO2':
        return NO2(atoms,i,state)
    elif substitution=='CHO':
        return CHO(atoms,i,state)
    elif substitution=='C=NH':
        return C2N(atoms,i,state,shuffle)
    elif substitution=='COCH3' or substitution=='COMe':
        return COMe(atoms,i,state)
    elif substitution=='COOH':
        return COOH(atoms,i,state)
    elif substitution=='COOMe' or substitution=='COOCH3':
        return COOMe(atoms,i,state)
    elif substitution=='CN':
        return CN(atoms,i,state)
    elif substitution=='Ph':
        return Ph(atoms,i,state,shuffle)
    elif substitution=='Py' or substitution=='pyridine' or substitution=='Pyridine':
        return Py(atoms,i,state,place,shuffle)
    elif substitution=='PL' or substitution=='pyrrole' or substitution=='Pyrrole':
        return Pyrrole(atoms,i,state,place,shuffle)
    elif substitution=='cycloC3H5' or substitution=='cycloPr' or  substitution=='cyclo3':
        return ring3(atoms,i,state,strict)
    elif substitution=='cycloC5H9' or substitution=='cycloPe' or substitution=='cyclo5':
        return ring5(atoms,i,state,strict,shuffle)
    elif substitution=='cycloC6H11' or substitution=='cycloHe' or substitution=='cyclo6':
        return ring6(atoms,i,state,strict,shuffle)
    elif substitution=='SOMe' or substitution=='SOCH3':
        return SOMe(atoms,i,state)
    elif substitution=='SO2Me' or substitution=='SO2CH3':
        return SO2_('CH3',atoms,i,state)
    elif substitution=='SO3H':
        return SO2_('OH',atoms,i,state)
    elif substitution=='double5' or substitution=='d5' or substitution=='D5' or substitution=='bicyclo[2.2.1]':
        return D(5,atoms,i,state,shuffle)
    elif substitution=='double6' or substitution=='d6' or substitution=='D6' or substitution=='bicyclo[2.2.2]':
        return D(6,atoms,i,state,shuffle)
    elif substitution=='Bpin':
        return Bpin(atoms,i,state,shuffle)
    else:
        assert 'Error! No such substitution!'
    
def Halogen(substitution,atoms,i,state=1):
    elements=atoms.get_chemical_symbols()
    elements[i]=substitution
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    if state==1:
        if substitution=='F':
            atoms.set_distance(i,Cnear,1.3605,fix=1)
        elif substitution=='Cl':
            atoms.set_distance(i,Cnear,1.781,fix=1)
        elif substitution=='Br':
            atoms.set_distance(i,Cnear,1.945,fix=1)
        elif substitution=='I':
            atoms.set_distance(i,Cnear,2.210,fix=1)
        return atoms
    elif state==2 or state==6:
        if substitution=='F':
            if state==6:
                atoms.set_distance(i,Cnear,1.354,fix=1)
            if state==2:
                atoms.set_distance(i,Cnear,1.340,fix=1)
        elif substitution=='Cl':
            atoms.set_distance(i,Cnear,1.725,fix=1)
        elif substitution=='Br':
            atoms.set_distance(i,Cnear,1.870,fix=1)
        elif substitution=='I':
            atoms.set_distance(i,Cnear,2.080,fix=1)
        return atoms

def OH(atoms,i,state=1,angle=-1):
    atoms.append('H')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-1)+[1]
    x=atoms.get_positions()
    x[-1][0]=x[-1][0]+50
    x[-1][1]=x[-1][1]+50
    x[-1][2]=x[-1][2]+50
    atoms.set_positions(x)
    elements=atoms.get_chemical_symbols()
    elements[i]='O'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    if state==1:
        atoms.set_distance(i,Cnear,1.410,fix=1)
    elif state==2 or state==6:
        atoms.set_distance(i,Cnear,1.364,fix=1)
    elif state==4:
        atoms.set_distance(i,Cnear,1.327,fix=1)
    elif state==164:
        atoms.set_distance(i,Cnear,1.670,fix=1)
    else:
        atoms.set_distance(i,Cnear,1.410,fix=1)
    atoms.set_distance(total_atoms-1,i,0.945,fix=1)
    atoms.set_angle(Cnear,i,total_atoms-1,103.5)
    if angle>0:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-1,angle,mask=mask)
    if state==2 or state==6:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-1,5)
    return atoms

def OMe(atoms,i,state=1,angle=-1):
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    total_atoms=len(atoms)
    atoms=OH(atoms,i,state,angle)
    total_atoms=len(atoms)
    atoms=CH3(atoms,total_atoms-1,8)
    #print(Cnear,CnearC)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-4)+[1]*4
    angle=0
    atoms.set_distance(i,total_atoms-1,2,mask=mask)
    atoms.set_dihedral(CnearC,Cnear,i,total_atoms-4,105,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-4,total_atoms-1,88,mask=mask)
    mask[-3:]=[0,0,0]
    atoms.set_distance(i,total_atoms-4,1.43)
    return atoms

def SH(atoms,i,state=1,angle=-1):
    atoms.append('H')
    total_atoms=len(atoms)
    x=atoms.get_positions()
    x[-1][0]=x[-1][0]+50
    x[-1][1]=x[-1][1]+50
    x[-1][2]=x[-1][2]+50
    atoms.set_positions(x)
    elements=atoms.get_chemical_symbols()
    elements[i]='S'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    if state==16:
        atoms.set_distance(i,Cnear,2.057,fix=1)
    else:
        atoms.set_distance(i,Cnear,1.810,fix=1)
    atoms.set_angle(Cnear,i,total_atoms-1,101.5)
    if angle>0:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-1,angle)
    atoms.set_distance(total_atoms-1,i,1.336,fix=1)
    return atoms

def SMe(atoms,i,state=1,angle=-1):
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    atoms=SH(atoms,i,state,angle)
    total_atoms=len(atoms)
    atoms=CH3(atoms,total_atoms-1,16)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-4)+[1]*4
    angle=0
    atoms.set_dihedral(CnearC,Cnear,i,total_atoms-4,105,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-4,total_atoms-1,88,mask=mask)
    mask[-3:]=[0,0,0]
    atoms.set_distance(i,total_atoms-4,1.82)
    return atoms

def NH2(atoms,i,state=1):
    atoms.append('H')
    atoms.append('H')
    total_atoms=len(atoms)
    x=atoms.get_positions()
    x[-2][0]=x[-2][0]+50
    x[-2][1]=x[-2][1]+50
    x[-2][2]=x[-2][2]+50
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    atoms.set_positions(x)
    elements=atoms.get_chemical_symbols()
    elements[i]='N'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    mask=[0]*(total_atoms-2)+[1]*2
    if state==1:
        atoms.set_distance(i,Cnear,1.448,fix=1)
    elif state==2 or state==6:
        atoms.set_distance(i,Cnear,1.340,fix=1)
    else:
        atoms.set_distance(i,Cnear,1.448,fix=1)
    atoms.set_distance(total_atoms-2,i,1.0100,fix=1)
    atoms.set_angle(Cnear,i,total_atoms-2,110.5)
    if state==1:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-2,70,mask=mask)
    elif state==2 or state==6:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-2,-175,mask=mask)
    mask[-2]=0
    atoms.set_angle(total_atoms-2,i,total_atoms-1,110.5)
    atoms.set_distance(total_atoms-1,i,1.0100,fix=1)
    atoms.set_angle(Cnear,i,total_atoms-1,110.5)
    if state==1:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-1,179,mask=mask)
    elif state==2 or state==6:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-1,70,mask=mask)
    #print(atoms.get_positions())
    return atoms

def NMe2(atoms,i):
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    atoms=NH2(atoms,i)
    total_atoms=len(atoms)
    atomA=total_atoms-2
    atomB=total_atoms-1
    mask=[0]*total_atoms
    mask[atomA]=1
    atoms.set_distance(i,atomA,1.448,fix=1,mask=mask)
    mask[atomA]=0
    atoms=CH3(atoms,atomA,7)
    total_atoms=len(atoms)
    mask=[0]*total_atoms
    mask[atomB]=1
    atoms.set_distance(i,atomB,1.448,fix=1,mask=mask)
    mask[atomB]=0
    atoms=CH3(atoms,atomB,7)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-8)+[1]*8
    atoms.set_dihedral(CnearC,Cnear,i,total_atoms-8,95,mask=mask)
    return atoms

def C2H3(atoms,i,state=1,shuffle=0):
    atoms.append('H')
    atoms.append('C')
    atoms.append('H')
    atoms.append('H')
    Cnear=near(i, atoms)
    total_atoms=len(atoms)
    elements=atoms.get_chemical_symbols()
    elements[i]='C'
    if state==1:
        atoms.set_distance(i,Cnear,1.5095,fix=1)
    elif state==2:
        atoms.set_distance(i,Cnear,1.4600,fix=1)
    elif state==6:
        atoms.set_distance(i,Cnear,1.4330,fix=1)
    elif state==82:
        atoms.set_distance(i,Cnear,1.4440,fix=1)
    elif state==162:
        atoms.set_distance(i,Cnear,1.7700,fix=1)
    elif state==164:
        atoms.set_distance(i,Cnear,1.7600,fix=1)
    else:
        atoms.set_distance(i,Cnear,1.5095,fix=1)
    CnearH=near(Cnear,atoms,'H',ignore=[Cnear,i])
    #print(CnearH)
    atoms.set_chemical_symbols(elements)
    mask=[0]*(total_atoms-4)+[1]*4
    x=atoms.get_positions()
    x[-4][0]=x[-4][0]+50
    x[-4][1]=x[-4][1]+50
    x[-4][2]=x[-4][2]+50
    for j in range(-3,0):
        x[j][0]=x[j][0]-50
        x[j][1]=x[j][1]-50
        x[j][2]=x[j][2]-50
    atoms.set_positions(x)
    atoms.set_angle(Cnear,i,total_atoms-4,120,mask=mask)
    atoms.set_distance(total_atoms-4,i,1.0800,fix=1,mask=mask)
    if shuffle==1:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-4,5.0,mask=mask)
    else:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-4,-175.0,mask=mask)
    mask[-4]=0
    atoms.set_angle(Cnear,i,total_atoms-3,120,mask=mask)
    atoms.set_distance(total_atoms-3,i,1.3400,fix=1,mask=mask)
    if shuffle==1:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-3,-175.0,mask=mask)
    else:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-3,5.0,mask=mask)
    mask[-3]=0
    atoms.set_angle(i,total_atoms-3,total_atoms-2,120,mask=mask)
    atoms.set_distance(total_atoms-2,total_atoms-3,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-3,total_atoms-2,0.5,mask=mask)
    mask[-2]=0
    atoms.set_angle(i,total_atoms-3,total_atoms-1,120,mask=mask)
    atoms.set_distance(total_atoms-1,total_atoms-3,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-3,total_atoms-1,-179.5,mask=mask)
    return atoms

def CH3(atoms,i,state=1,ring=False):
    atoms.append('H')
    atoms.append('H')
    atoms.append('H')
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    total_atoms=len(atoms)
    elements=atoms.get_chemical_symbols()
    elements[i]='C'
    if state==1:
        atoms.set_distance(i,Cnear,1.5290,fix=1)
    elif state==2:
        atoms.set_distance(i,Cnear,1.5095,fix=1)
    elif state==3:
        atoms.set_distance(i,Cnear,1.4700,fix=1)
    elif state==6:
        atoms.set_distance(i,Cnear,1.5105,fix=1)
    elif state==7:
        atoms.set_distance(i,Cnear,1.4480,fix=1)
    elif state==8:
        atoms.set_distance(i,Cnear,1.4100,fix=1)
    elif state==82:
        atoms.set_distance(i,Cnear,1.5220,fix=1)
    elif state==16:
        atoms.set_distance(i,Cnear,1.8100,fix=1)
    elif state==162:
        atoms.set_distance(i,Cnear,1.7900,fix=1)
    elif state==164:
        atoms.set_distance(i,Cnear,1.7700,fix=1)
    else:
        atoms.set_distance(i,Cnear,1.5290,fix=1)
    atoms.set_chemical_symbols(elements)
    mask=[0]*(total_atoms-3)+[1]*3
    x=atoms.get_positions()
    x[-3][0]=x[-3][0]-50
    x[-3][1]=x[-3][1]-50
    x[-3][2]=x[-3][2]-50
    x[-2][0]=x[-2][0]-50
    x[-2][1]=x[-2][1]-50
    x[-2][2]=x[-2][2]-50
    atoms.set_positions(x)
    atoms.set_angle(Cnear,i,total_atoms-3,109.5,mask=mask)
    atoms.set_distance(total_atoms-3,i,1.0900,fix=1,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-3,109.5,mask=mask)
    if ring==True:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-3,120,mask=mask)
    else:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-3,179,mask=mask)
    mask[-3]=0
    atoms.set_angle(total_atoms-3,i,total_atoms-2,104.5,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-2,109.5,mask=mask)
    atoms.set_distance(total_atoms-2,i,1.0900,fix=1,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-2,109.5,mask=mask)
    if ring==True:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-2,0.5,mask=mask)
    else:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-2,59,mask=mask)
    mask[-2]=0
    x=atoms.get_positions()
    x[-1][0]=x[-1][0]+50
    x[-1][1]=x[-1][1]+50
    x[-1][2]=x[-1][2]+50
    atoms.set_positions(x)
    atoms.set_angle(Cnear,i,total_atoms-1,109.5,mask=mask)
    atoms.set_angle(total_atoms-3,i,total_atoms-1,109.5,mask=mask)
    atoms.set_distance(total_atoms-1,i,1.0900,fix=1,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-1,109.5,mask=mask)
    if ring==True:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-1,-120,mask=mask)
    else:
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-1,-60,mask=mask)
    return atoms

def CChainL(atoms,i,chain,state=1):
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    atoms=CH3(atoms,i,state)
    for i in range(chain):
        total_atoms=len(atoms)
        mask=[0]*total_atoms
        total_atoms=len(atoms)
        atom=total_atoms-2
        mask[atom]=1
        atoms.set_distance(i,atom,1.5290,fix=1,mask=mask)
        mask[atom]=0
        atoms=CH3(atoms,atom)
    return atoms

def CChainR(atoms,i,chain,state=1):
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    atoms=CH3(atoms,i,state)
    #view(atoms)
    total_atoms=len(atoms)
    atomA=total_atoms-3
    atomB=total_atoms-2
    atomC=total_atoms-1
    mask=[0]*total_atoms
    if chain>=1:
        mask[atomB]=1
        atoms.set_distance(i,atomB,1.5290,fix=1,mask=mask)
        mask[atomB]=0
        atoms=CH3(atoms,atomB)
        if chain>=2:
            total_atoms=len(atoms)
            mask=[0]*total_atoms
            mask[atomC]=1
            atoms.set_distance(i,atomC,1.5290,fix=1,mask=mask)
            mask[atomC]=0
            atoms=CH3(atoms,atomC)
            if chain==3:
                total_atoms=len(atoms)
                mask=[0]*total_atoms
                mask[atomA]=1
                atoms.set_distance(i,atomA,1.5290,fix=1,mask=mask)
                mask[atomA]=0
                atoms=CH3(atoms,atomA)
    total_atoms=len(atoms)
    mask=[0]*total_atoms
    if chain==1:
        mask[-6:]=[1]*6
        atoms.set_dihedral(CnearC,Cnear,i,atomB,110,mask=mask)
    if chain==2:
        mask[-9:]=[1]*9
        atoms.set_dihedral(CnearC,Cnear,i,atomB,95,mask=mask)
    return atoms

def trihalo(substitution,atoms,i,state=1):
    import re
    atoms=CH3(atoms,i,state)
    halo=re.search(r'[A-Z]+[a-z]*',substitution)
    for i in range(-3,0):
        Halogen(halo.group(0)[1:],atoms,i)
    return atoms

def NO2(atoms,i,state=1):
    atoms.append('O')
    atoms.append('O')
    total_atoms=len(atoms)
    elements=atoms.get_chemical_symbols()
    elements[i]='N'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    CnearH=near(Cnear,atoms,'H')
    if state==1:
        atoms.set_distance(i,Cnear,1.4900,fix=1)
    elif state==6:
        atoms.set_distance(i,Cnear,1.4600,fix=1)
    atoms.set_dihedral(CnearH,Cnear,i,total_atoms-2,-175.0,mask=mask)
    atoms.set_distance(total_atoms-2,i,1.2250,fix=1)
    atoms.set_angle(Cnear,i,total_atoms-2,119.5)
    atoms.set_angle(total_atoms-2,i,total_atoms-1,119.5)
    atoms.set_distance(total_atoms-1,i,1.2250,fix=1)
    atoms.set_angle(Cnear,i,total_atoms-1,120.5)
    #print(atoms.get_positions())
    return atoms

def CHO(atoms,i,state=1):
    atoms.append('O')
    atoms.append('H')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-2)+[1]*2
    elements=atoms.get_chemical_symbols()
    elements[i]='C'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    CnearH=near(Cnear,atoms,'H')
    if state==1:
        atoms.set_distance(i,Cnear,1.5220,fix=1)
    elif state==2:
        atoms.set_distance(i,Cnear,1.4440,fix=1)
    elif state==6:
        atoms.set_distance(i,Cnear,1.4900,fix=1)
    elif state==7:
        atoms.set_distance(i,Cnear,1.3350,fix=1)
    elif state==8:
        atoms.set_distance(i,Cnear,1.3640,fix=1)
    elif state==82:
        atoms.set_distance(i,Cnear,1.5100,fix=1)
    x=atoms.get_positions()
    #Separate into a side
    x[-2][0]=x[-2][0]-50
    x[-2][1]=x[-2][1]-50
    x[-2][2]=x[-2][2]-50
    atoms.set_positions(x)
    atoms.set_angle(Cnear,i,total_atoms-2,120.5,mask=mask)
    atoms.set_distance(total_atoms-2,i,1.2290,fix=1,mask=mask)
    atoms.set_dihedral(CnearH,Cnear,i,total_atoms-2,55,mask=mask)#
    mask[-2]=0
    x=atoms.get_positions()
    #Separate into another side
    x[-1][0]=x[-1][0]+50
    x[-1][1]=x[-1][1]+50
    x[-1][2]=x[-1][2]+50
    atoms.set_positions(x)
    for k in range(10):
        atoms.set_angle(total_atoms-2,i,total_atoms-1,120,mask=mask)
        atoms.set_angle(Cnear,i,total_atoms-1,120,mask=mask)
    atoms.set_distance(total_atoms-1,i,1.0900,fix=1,mask=mask)
    return atoms

def C2N(atoms,i,state=1,shuffle=0):
    atoms.append('N')
    atoms.append('H')
    atoms.append('H')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-3)+[1]*3
    elements=atoms.get_chemical_symbols()
    elements[i]='C'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    CnearH=near(Cnear,atoms,'H',ignore=[Cnear,i])
    if state==1:
        atoms.set_distance(i,Cnear,1.5220,fix=1)
    elif state==2:
        atoms.set_distance(i,Cnear,1.4440,fix=1)
    elif state==6:
        atoms.set_distance(i,Cnear,1.5220,fix=1)
    elif state==7:
        atoms.set_distance(i,Cnear,1.3350,fix=1)
    elif state==8:
        atoms.set_distance(i,Cnear,1.3640,fix=1)
    elif state==82:
        atoms.set_distance(i,Cnear,1.5100,fix=1)
    x=atoms.get_positions()
    #Separate into a side
    x[-3][0]=x[-3][0]-50
    x[-3][1]=x[-3][1]-50
    x[-3][2]=x[-3][2]-50
    atoms.set_positions(x)
    atoms.set_angle(Cnear,i,total_atoms-3,120.5,mask=mask)
    atoms.set_distance(total_atoms-3,i,1.3580,fix=1,mask=mask)
    if shuffle==0:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-3,55,mask=mask)
    elif shuffle==1:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-3,179.5,mask=mask)
    mask[-3]=0
    x=atoms.get_positions()
    x[-2][0]=x[-2][0]+50
    x[-2][1]=x[-2][1]+50
    x[-2][2]=x[-2][2]+50
    atoms.set_positions(x)
    atoms.set_angle(i,total_atoms-3,total_atoms-2,120.5,mask=mask)
    atoms.set_distance(total_atoms-2,total_atoms-3,1.0100,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-3,total_atoms-2,179.5,mask=mask)
    mask[-2]=0
    x=atoms.get_positions()
    #Separate into another side
    x[-1][0]=x[-1][0]+50
    x[-1][1]=x[-1][1]+50
    x[-1][2]=x[-1][2]+50
    atoms.set_positions(x)
    if shuffle==1:
        atoms.set_dihedral(Cnear,i,total_atoms-3,total_atoms-1,55,mask=mask)
    for k in range(10):
        atoms.set_angle(total_atoms-3,i,total_atoms-1,120,mask=mask)
        atoms.set_angle(Cnear,i,total_atoms-1,120,mask=mask)
    atoms.set_distance(total_atoms-1,i,1.0100,fix=1,mask=mask)
    return atoms

def COMe(atoms,i,state=1):
    atoms=CHO(atoms,i,state)
    total_atoms=len(atoms)
    atoms=CH3(atoms,total_atoms-1,82)
    total_atoms=len(atoms)
    return atoms
    
def COOH(atoms,i,state=1):
    atoms.append('O')
    atoms.append('O')
    atoms.append('H')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-3)+[1]*3
    elements=atoms.get_chemical_symbols()
    elements[i]='C'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    CnearH=near(Cnear,atoms,'H')
    if state==1:
        atoms.set_distance(i,Cnear,1.5220,fix=1)
    elif state==2:
        atoms.set_distance(i,Cnear,1.4440,fix=1)
    elif state==6:
        atoms.set_distance(i,Cnear,1.4900,fix=1)
    elif state==7:
        atoms.set_distance(i,Cnear,1.3350,fix=1)
    elif state==8:
        atoms.set_distance(i,Cnear,1.3640,fix=1)
    elif state==82:
        atoms.set_distance(i,Cnear,1.5095,fix=1)
    x=atoms.get_positions()
    #Separate into a side
    x[-3][0]=x[-3][0]-50
    x[-3][1]=x[-3][1]-50
    x[-3][2]=x[-3][2]-50
    atoms.set_positions(x)
    atoms.set_angle(Cnear,i,total_atoms-3,120.5,mask=mask)
    atoms.set_distance(total_atoms-3,i,1.3640,fix=1,mask=mask)
    atoms.set_dihedral(CnearH,Cnear,i,total_atoms-3,10.4,mask=mask)#
    mask[-3]=0
    x=atoms.get_positions()
    #Separate into another side
    x[-2][0]=x[-2][0]+50
    x[-2][1]=x[-2][1]+50
    x[-2][2]=x[-2][2]+50
    atoms.set_positions(x)
    for k in range(10):
        atoms.set_angle(total_atoms-3,i,total_atoms-2,120,mask=mask)
        atoms.set_angle(Cnear,i,total_atoms-2,120,mask=mask)
    atoms.set_distance(total_atoms-2,i,1.2290,fix=1,mask=mask)
    mask[-2]=0
    atoms.set_distance(total_atoms-1,total_atoms-3,0.945,fix=1)
    atoms.set_angle(i,total_atoms-3,total_atoms-1,120.5)
    atoms.set_dihedral(total_atoms-2,i,total_atoms-3,total_atoms-1,5)
    #print(atoms.get_positions())
    #print(atoms.get_chemical_symbols())
    return atoms

def COOMe(atoms,i,state=1):
    atoms=COOH(atoms,i,state)
    total_atoms=len(atoms)
    atoms=CH3(atoms,total_atoms-1,8)
    return atoms

def CN(atoms,i,state=1):
    atoms.append('N')
    total_atoms=len(atoms)
    elements=atoms.get_chemical_symbols()
    elements[i]='C'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    if state==1:
        atoms.set_distance(i,Cnear,1.4700,fix=1)
    elif state==6:
        atoms.set_distance(i,Cnear,1.4510,fix=1)
    atoms.set_distance(total_atoms-1,i,1.157,fix=1)
    atoms.set_angle(Cnear,i,total_atoms-1,179.5)
    #print(atoms.get_positions())
    return atoms

def Ph(atoms,i,state=1,shuffle=0):
    for j in range(5):
        atoms.append('C')
    for j in range(5):
        atoms.append('H')
    Cnear=near(i, atoms)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-10)+[1]*10
    elements=atoms.get_chemical_symbols()
    elements[i]='C'
    atoms.set_chemical_symbols(elements)
    x=atoms.get_positions()
    #Separate into a side
    for j in range(-10,-8):
        x[j][0]=x[j][0]-50
        x[j][1]=x[j][1]-50
        x[j][2]=x[j][2]-50
    for j in range(-5,-3):
        x[j][0]=x[j][0]-50
        x[j][1]=x[j][1]-50
        x[j][2]=x[j][2]-50
    for j in range(-8,-5):
        x[j][0]=x[j][0]+50
        x[j][1]=x[j][1]+50
        x[j][2]=x[j][2]+50
    for j in range(-3,0):
        x[j][0]=x[j][0]+50
        x[j][1]=x[j][1]+50
        x[j][2]=x[j][2]+50
    atoms.set_positions(x)
    #if shuffle==0:
    CnearH=near(Cnear,atoms,'H',ignore=[Cnear,i])
    #elif shuffle==1:
    #    CnearH=near(Cnear,atoms,'C')
    #print(Cnear,CnearH)
    if state==1:
        atoms.set_distance(i,Cnear,1.5105,fix=1)
    elif state==2:
        atoms.set_distance(i,Cnear,1.4330,fix=1)
    elif state==6:
        atoms.set_distance(i,Cnear,1.5100,fix=1)
    elif state==7:
        atoms.set_distance(i,Cnear,1.3400,fix=1)
    elif state==8:
        atoms.set_distance(i,Cnear,1.3640,fix=1)
    elif state==16:
        atoms.set_distance(i,Cnear,1.7700,fix=1)
    elif state==82:
        atoms.set_distance(i,Cnear,1.4900,fix=1)
    atoms.set_angle(Cnear,i,total_atoms-10,120,mask=mask)
    atoms.set_distance(total_atoms-10,i,1.4,fix=1,mask=mask)
    if shuffle==0:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-10,5.4,mask=mask)
    elif shuffle==1:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-10,95.4,mask=mask)
    mask[-10]=0
    atoms.set_angle(Cnear,i,total_atoms-6,120,mask=mask)
    atoms.set_distance(total_atoms-6,i,1.4,fix=1,mask=mask)
    if shuffle==0:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-6,-175.0,mask=mask)
    elif shuffle==1:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-6,-85.0,mask=mask)
    mask[-6]=0
    atoms.set_angle(i,total_atoms-10,total_atoms-9,120,mask=mask)
    atoms.set_distance(total_atoms-9,total_atoms-10,1.4,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-10,total_atoms-9,179.5,mask=mask)
    mask[-9]=0
    atoms.set_angle(i,total_atoms-6,total_atoms-7,120,mask=mask)
    atoms.set_distance(total_atoms-7,total_atoms-6,1.4,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-6,total_atoms-7,179.5,mask=mask)
    mask[-7]=0
    atoms.set_angle(total_atoms-6,total_atoms-7,total_atoms-8,120,mask=mask)
    atoms.set_distance(total_atoms-8,total_atoms-7,1.4,fix=1,mask=mask)
    atoms.set_angle(total_atoms-7,total_atoms-8,total_atoms-9,120,mask=mask)
    atoms.set_dihedral(i,total_atoms-6,total_atoms-7,total_atoms-8,0.5,mask=mask)
    mask[-8]=0
    atoms.set_angle(i,total_atoms-10,total_atoms-5,120,mask=mask)
    atoms.set_distance(total_atoms-5,total_atoms-10,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-10,total_atoms-5,0.5,mask=mask)
    mask[-5]=0
    atoms.set_angle(i,total_atoms-6,total_atoms-1,120,mask=mask)
    atoms.set_distance(total_atoms-1,total_atoms-6,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-6,total_atoms-1,0.5,mask=mask)
    mask[-1]=0
    atoms.set_angle(total_atoms-10,total_atoms-9,total_atoms-4,120,mask=mask)
    atoms.set_distance(total_atoms-4,total_atoms-9,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-5,total_atoms-10,total_atoms-9,total_atoms-4,0.5,mask=mask)
    mask[-4]=0
    atoms.set_angle(total_atoms-6,total_atoms-7,total_atoms-2,120,mask=mask)
    atoms.set_distance(total_atoms-2,total_atoms-7,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-1,total_atoms-6,total_atoms-7,total_atoms-2,0.5,mask=mask)
    mask[-2]=0
    atoms.set_angle(total_atoms-9,total_atoms-8,total_atoms-3,120,mask=mask)
    atoms.set_distance(total_atoms-3,total_atoms-8,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-4,total_atoms-9,total_atoms-8,total_atoms-3,0.5,mask=mask)
    return atoms

def Py(atoms,i,state=1,place=4,shuffle=0):
    atoms=Ph(atoms,i,state,shuffle)
    total_atoms=len(atoms)
    mask=[0]*total_atoms

def Pyrrole(atoms,i,state=1,place=0,shuffle=0):
    if place==0:
        for j in range(4):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear=near(i, atoms)
        total_atoms=len(atoms)
        mask=[0]*(total_atoms-8)+[1]*8
        elements=atoms.get_chemical_symbols()
        elements[i]='N'
        if shuffle==0:
            CnearH=near(Cnear,atoms,'H')
        elif shuffle==1:
            CnearH=near(Cnear,atoms,'C')
        atoms.set_chemical_symbols(elements)
        if state==1:
            atoms.set_distance(i,Cnear,1.475,fix=1,ignore=i)
        elif state==6:
            atoms.set_distance(i,Cnear,1.440,fix=1,ignore=i)
        x=atoms.get_positions()
        #Separate into a side
        for j in range(-8,-6):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-4,-2):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-6,-4):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        for j in range(-2,0):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        atoms.set_positions(x)
        atoms.set_angle(Cnear,i,total_atoms-8,125.35,mask=mask)
        atoms.set_distance(total_atoms-8,i,1.381,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-8,5.4,mask=mask)
        mask[-8]=0
        atoms.set_angle(Cnear,i,total_atoms-5,125.35,mask=mask)
        atoms.set_distance(total_atoms-5,i,1.381,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-5,-175.0,mask=mask)
        mask[-5]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-7,107.7,mask=mask)
        atoms.set_distance(total_atoms-7,total_atoms-8,1.367,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-7,179.5,mask=mask)
        mask[-7]=0
        atoms.set_angle(i,total_atoms-5,total_atoms-6,107.7,mask=mask)
        atoms.set_distance(total_atoms-6,total_atoms-5,1.367,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-6,-179.5,mask=mask)
        mask[-6]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-4,126,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-8,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-4,0.5,mask=mask)
        mask[-4]=0
        atoms.set_angle(i,total_atoms-5,total_atoms-1,126,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-5,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-1,-0.5,mask=mask)
        mask[-1]=0
        atoms.set_angle(total_atoms-8,total_atoms-7,total_atoms-3,126,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-7,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-8,total_atoms-7,total_atoms-3,179.5,mask=mask)
        mask[-3]=0
        atoms.set_angle(total_atoms-5,total_atoms-6,total_atoms-2,126,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-6,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-5,total_atoms-6,total_atoms-2,-179.5,mask=mask)
        mask[-2]=0
    elif place==1:
        atoms.append('N')
        for j in range(3):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear=near(i, atoms)
        total_atoms=len(atoms)
        mask=[0]*(total_atoms-8)+[1]*8
        elements=atoms.get_chemical_symbols()
        elements[i]='C'
        if shuffle==0:
            CnearH=near(Cnear,atoms,'H',ignore=i)
        elif shuffle==1:
            CnearH=near(Cnear,atoms,'C',ignore=i)
        atoms.set_chemical_symbols(elements)
        if state==1:
            atoms.set_distance(i,Cnear,1.5100,fix=1)
        elif state==6:
            atoms.set_distance(i,Cnear,1.4600,fix=1)
        x=atoms.get_positions()
        #Separate into a side
        for j in range(-8,-6):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-4,-2):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-6,-4):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        for j in range(-2,0):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        atoms.set_positions(x)
        atoms.set_angle(Cnear,i,total_atoms-8,120.8,mask=mask)
        atoms.set_distance(total_atoms-8,i,1.381,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-8,6.0,mask=mask)
        mask[-8]=0
        atoms.set_angle(Cnear,i,total_atoms-5,131.58,mask=mask)
        atoms.set_distance(total_atoms-5,i,1.367,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-5,-175.0,mask=mask)
        mask[-5]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-7,109.8,mask=mask)
        atoms.set_distance(total_atoms-7,total_atoms-8,1.381,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-7,179.5,mask=mask)
        mask[-7]=0
        atoms.set_angle(i,total_atoms-5,total_atoms-6,107.3,mask=mask)
        atoms.set_distance(total_atoms-6,total_atoms-5,1.424,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-6,-179.5,mask=mask)
        mask[-6]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-4,126,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-8,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-4,0.5,mask=mask)
        mask[-4]=0
        atoms.set_angle(i,total_atoms-5,total_atoms-1,126,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-5,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-1,-0.5,mask=mask)
        mask[-1]=0
        atoms.set_angle(total_atoms-8,total_atoms-7,total_atoms-3,126,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-7,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-8,total_atoms-7,total_atoms-3,179.5,mask=mask)
        mask[-3]=0
        atoms.set_angle(total_atoms-5,total_atoms-6,total_atoms-2,126,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-6,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-5,total_atoms-6,total_atoms-2,-179.5,mask=mask)
        mask[-2]=0
    elif place==-1:
        atoms.append('N')
        for j in range(3):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear=near(i, atoms)
        total_atoms=len(atoms)
        mask=[0]*(total_atoms-8)+[1]*8
        elements=atoms.get_chemical_symbols()
        elements[i]='C'
        if shuffle==0:
            CnearH=near(Cnear,atoms,'H',ignore=i)
        elif shuffle==1:
            CnearH=near(Cnear,atoms,'C',ignore=i)
        atoms.set_chemical_symbols(elements)
        if state==1:
            atoms.set_distance(i,Cnear,1.5100,fix=1)
        elif state==6:
            atoms.set_distance(i,Cnear,1.4600,fix=1)
        x=atoms.get_positions()
        #Separate into a side
        for j in range(-8,-6):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-4,-2):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-6,-4):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        for j in range(-2,0):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        atoms.set_positions(x)
        atoms.set_angle(Cnear,i,total_atoms-5,131.58,mask=mask)
        atoms.set_distance(total_atoms-5,i,1.367,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-5,6.0,mask=mask)
        mask[-5]=0
        atoms.set_angle(Cnear,i,total_atoms-8,120.8,mask=mask)
        atoms.set_distance(total_atoms-8,i,1.381,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-8,-175.0,mask=mask)
        mask[-8]=0
        atoms.set_angle(i,total_atoms-5,total_atoms-6,107.3,mask=mask)
        atoms.set_distance(total_atoms-6,total_atoms-5,1.424,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-6,179.5,mask=mask)
        mask[-6]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-7,109.8,mask=mask)
        atoms.set_distance(total_atoms-7,total_atoms-8,1.381,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-7,-179.5,mask=mask)
        mask[-7]=0
        atoms.set_angle(i,total_atoms-7,total_atoms-4,126,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-5,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-4,0.5,mask=mask)
        mask[-4]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-1,126,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-8,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-1,-0.5,mask=mask)
        mask[-1]=0
        atoms.set_angle(total_atoms-5,total_atoms-6,total_atoms-3,126,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-6,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-5,total_atoms-6,total_atoms-3,179.5,mask=mask)
        mask[-3]=0
        atoms.set_angle(total_atoms-8,total_atoms-7,total_atoms-2,126,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-7,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-8,total_atoms-7,total_atoms-2,-179.5,mask=mask)
        mask[-2]=0
    elif place==2:
        atoms.append('C')
        atoms.append('N')
        for j in range(2):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear=near(i, atoms)
        total_atoms=len(atoms)
        mask=[0]*(total_atoms-8)+[1]*8
        elements=atoms.get_chemical_symbols()
        elements[i]='C'
        if shuffle==0:
            CnearH=near(Cnear,atoms,'H',ignore=i)
        elif shuffle==1:
            CnearH=near(Cnear,atoms,'C',ignore=i)
        atoms.set_chemical_symbols(elements)
        if state==1:
            atoms.set_distance(i,Cnear,1.5100,fix=1)
        elif state==6:
            atoms.set_distance(i,Cnear,1.4600,fix=1)
        x=atoms.get_positions()
        #Separate into a side
        for j in range(-8,-6):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-4,-2):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-6,-4):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        for j in range(-2,0):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        atoms.set_positions(x)
        atoms.set_angle(Cnear,i,total_atoms-8,125.4,mask=mask)
        atoms.set_distance(total_atoms-8,i,1.367,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-8,5.4,mask=mask)
        mask[-8]=0
        atoms.set_angle(Cnear,i,total_atoms-5,126.75,mask=mask)
        atoms.set_distance(total_atoms-5,i,1.424,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-5,-175.0,mask=mask)
        mask[-5]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-7,107.7,mask=mask)
        atoms.set_distance(total_atoms-7,total_atoms-8,1.381,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-7,179.5,mask=mask)
        mask[-7]=0
        atoms.set_angle(i,total_atoms-5,total_atoms-6,107.3,mask=mask)
        atoms.set_distance(total_atoms-6,total_atoms-5,1.367,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-6,-179.5,mask=mask)
        mask[-6]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-4,126,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-8,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-4,0.5,mask=mask)
        mask[-4]=0
        atoms.set_angle(i,total_atoms-5,total_atoms-1,126,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-5,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-1,-0.5,mask=mask)
        mask[-1]=0
        atoms.set_angle(total_atoms-8,total_atoms-7,total_atoms-3,126,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-7,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-8,total_atoms-7,total_atoms-3,179.5,mask=mask)
        mask[-3]=0
        atoms.set_angle(total_atoms-5,total_atoms-6,total_atoms-2,126,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-6,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-5,total_atoms-6,total_atoms-2,-179.5,mask=mask)
        mask[-2]=0
    elif place==-2:
        atoms.append('C')
        atoms.append('N')
        for j in range(2):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        Cnear=near(i, atoms)
        total_atoms=len(atoms)
        mask=[0]*(total_atoms-8)+[1]*8
        elements=atoms.get_chemical_symbols()
        elements[i]='C'
        if shuffle==0:
            CnearH=near(Cnear,atoms,'H',ignore=i)
        elif shuffle==1:
            CnearH=near(Cnear,atoms,'C',ignore=i)
        atoms.set_chemical_symbols(elements)
        if state==1:
            atoms.set_distance(i,Cnear,1.5100,fix=1)
        elif state==6:
            atoms.set_distance(i,Cnear,1.4600,fix=1)
        x=atoms.get_positions()
        #Separate into a side
        for j in range(-8,-6):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-4,-2):
            x[j][0]=x[j][0]-50
            x[j][1]=x[j][1]-50
            x[j][2]=x[j][2]-50
        for j in range(-6,-4):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        for j in range(-2,0):
            x[j][0]=x[j][0]+50
            x[j][1]=x[j][1]+50
            x[j][2]=x[j][2]+50
        atoms.set_positions(x)
        atoms.set_angle(Cnear,i,total_atoms-8,125.4,mask=mask)
        atoms.set_distance(total_atoms-8,i,1.367,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-8,175.0,mask=mask)
        mask[-8]=0
        atoms.set_angle(Cnear,i,total_atoms-5,126.75,mask=mask)
        atoms.set_distance(total_atoms-5,i,1.424,fix=1,mask=mask)
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-5,-5.4,mask=mask)
        mask[-5]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-7,107.7,mask=mask)
        atoms.set_distance(total_atoms-7,total_atoms-8,1.381,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-7,-179.5,mask=mask)
        mask[-7]=0
        atoms.set_angle(i,total_atoms-5,total_atoms-6,107.3,mask=mask)
        atoms.set_distance(total_atoms-6,total_atoms-5,1.367,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-6,179.5,mask=mask)
        mask[-6]=0
        atoms.set_angle(i,total_atoms-5,total_atoms-4,126,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-5,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-5,total_atoms-4,0.5,mask=mask)
        mask[-4]=0
        atoms.set_angle(i,total_atoms-8,total_atoms-1,126,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-8,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms-8,total_atoms-1,-0.5,mask=mask)
        mask[-1]=0
        atoms.set_angle(total_atoms-8,total_atoms-7,total_atoms-3,126,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-7,1.0800,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-8,total_atoms-7,total_atoms-3,179.5,mask=mask)
        mask[-3]=0
        atoms.set_angle(total_atoms-5,total_atoms-6,total_atoms-2,126,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-6,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms-5,total_atoms-6,total_atoms-2,-179.5,mask=mask)
        mask[-2]=0
    return atoms

def ring3(atoms,i,state=1,strict={1:'X'}):
    for j in range(2):
        atoms.append('C')
    for j in range(5):
        atoms.append('H')
    total_atoms=len(atoms)
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    mask=[0]*(total_atoms-7)+[1]*7
    elements=atoms.get_chemical_symbols()
    elements[i]='C'
    try:
        if strict[1]=='N':
            elements[i]='N'
        else:
            elements[i]='C'
    except:
        elements[i]='C'
    atoms.set_chemical_symbols(elements)
    if elements[i]=='C':
        if state==1:
            atoms.set_distance(i,Cnear,1.5290,fix=1)
        elif state==2 or state==82:
            atoms.set_distance(i,Cnear,1.5095,fix=1)
        elif state==3:
            atoms.set_distance(i,Cnear,1.4700,fix=1)
        elif state==6:
            atoms.set_distance(i,Cnear,1.5105,fix=1)
        elif state==7:
            atoms.set_distance(i,Cnear,1.4480,fix=1)
        elif state==8:
            atoms.set_distance(i,Cnear,1.4100,fix=1)
        elif state==16:
            atoms.set_distance(i,Cnear,1.8100,fix=1)
        elif state==162:
            atoms.set_distance(i,Cnear,1.7900,fix=1)
        elif state==164:
            atoms.set_distance(i,Cnear,1.7700,fix=1)
    elif elements[i]=='N':
        if state==1:
            atoms.set_distance(i,Cnear,1.448,fix=1)
        elif state==2 or state==6:
            atoms.set_distance(i,Cnear,1.340,fix=1)
        else:
            atoms.set_distance(i,Cnear,1.448,fix=1)
    x=atoms.get_positions()
    x[-5][0]=x[-5][0]-50
    x[-5][1]=x[-5][1]-50
    x[-5][2]=x[-5][2]-50
    x[-7][0]=x[-7][0]-50
    x[-7][1]=x[-7][1]-50
    x[-7][2]=x[-7][2]-50
    x[-6][0]=x[-6][0]+50
    x[-6][1]=x[-6][1]+50
    x[-6][2]=x[-6][2]+50
    for j in range(-4,-2):
        x[j][0]=x[j][0]-50
        x[j][1]=x[j][1]-50
        x[j][2]=x[j][2]-50
    for j in range(-2,0):
        x[j][0]=x[j][0]+50
        x[j][1]=x[j][1]+50
        x[j][2]=x[j][2]+50
    atoms.set_positions(x)
    if elements[i]=='C':
        atoms.set_angle(Cnear,i,total_atoms-5,109.5,mask=mask)
        atoms.set_distance(total_atoms-5,i,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-5,-150,mask=mask)
        mask[-5]=0
    else:
        mask[-5]=0
    if elements[i]=='N':
        atoms.set_angle(Cnear,i,total_atoms-7,122.5,mask=mask)
        atoms.set_distance(total_atoms-7,i,1.4480,fix=1,mask=mask)
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-7,0.5,mask=mask)
        mask[-7]=0
        atoms.set_angle(Cnear,i,total_atoms-6,122.5,mask=mask)
        atoms.set_distance(total_atoms-6,i,1.4480,fix=1,mask=mask)
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-6,73.5,mask=mask)
        mask[-6]=0
    elif elements[i]=='C':
        atoms.set_angle(Cnear,i,total_atoms-7,122.5,mask=mask)
        atoms.set_distance(total_atoms-7,i,1.5290,fix=1,mask=mask)
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-7,0.5,mask=mask)
        mask[-7]=0
        atoms.set_angle(Cnear,i,total_atoms-6,122.5,mask=mask)
        atoms.set_distance(total_atoms-6,i,1.5290,fix=1,mask=mask)
        atoms.set_dihedral(CnearC,Cnear,i,total_atoms-6,73.5,mask=mask)
        mask[-6]=0
    atoms.set_angle(i,total_atoms-7,total_atoms-4,122.5,mask=mask)
    atoms.set_distance(total_atoms-4,total_atoms-7,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-7,total_atoms-4,5,mask=mask)
    mask[-4]=0
    atoms.set_angle(i,total_atoms-7,total_atoms-3,122.5,mask=mask)
    atoms.set_distance(total_atoms-3,total_atoms-7,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-7,total_atoms-3,-140,mask=mask)
    mask[-3]=0
    atoms.set_angle(i,total_atoms-6,total_atoms-2,122.5,mask=mask)
    atoms.set_distance(total_atoms-2,total_atoms-6,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-6,total_atoms-2,-5,mask=mask)
    mask[-2]=0
    atoms.set_angle(i,total_atoms-6,total_atoms-1,122.5,mask=mask)
    atoms.set_distance(total_atoms-1,total_atoms-6,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-6,total_atoms-1,140,mask=mask)
    mask[-1]=0
    for j in range(2):
        try:
            ring_element=strict[j+2]
        except: 
            continue
        elements=atoms.get_chemical_symbols()
        if ring_element=='O':
            elements[-7+j]='O'
            atoms.set_chemical_symbols(elements)
            mask=[0]*total_atoms
            mask[-7+j]=1
            if j==0:
                for l in range(10):
                    atoms.set_distance(total_atoms-6,total_atoms-7,1.4100,mask=mask)
                    atoms.set_distance(i,total_atoms-7,1.4100,mask=mask)
            elif j==1:
                for l in range(10):
                    atoms.set_distance(i,total_atoms-6,1.4100,mask=mask)
                    atoms.set_distance(total_atoms-7,total_atoms-6,1.4100,mask=mask)
        elif ring_element=='S':
            elements[-7+j]='S'
            atoms.set_chemical_symbols(elements)
            mask=[0]*total_atoms
            mask[-7+j]=1
            if j==0:
                for l in range(10):
                    atoms.set_distance(total_atoms-6,total_atoms-7,1.8100,mask=mask)
                    atoms.set_distance(i,total_atoms-7,1.8100,mask=mask)
            elif j==1:
                for l in range(10):
                    atoms.set_distance(i,total_atoms-6,1.8100,mask=mask)
                    atoms.set_distance(total_atoms-7,total_atoms-6,1.8100,mask=mask)
        elif ring_element=='N':
            elements[-7+j]='N'
            atoms.set_chemical_symbols(elements)
            mask=[0]*total_atoms
            mask[-7+j]=1
            if j==0:
                for l in range(10):
                    atoms.set_distance(total_atoms-6,total_atoms-7,1.4480,mask=mask)
                    atoms.set_distance(i,total_atoms-7,1.4480,mask=mask)
                mask[-7]=0
                mask[-3]=1
                for l in range(10):
                    atoms.set_angle(i,total_atoms-7,total_atoms-3,109.5,mask=mask)
                    atoms.set_angle(total_atoms-6,total_atoms-7,total_atoms-3,109.5,mask=mask)
                atoms.set_distance(total_atoms-3,total_atoms-7,1.0100,fix=1,mask=mask)
                mask[-3]=0
            elif j==1:
                for l in range(10):
                    atoms.set_distance(i,total_atoms-6,1.4480,mask=mask)
                    atoms.set_distance(total_atoms-7,total_atoms-6,1.4480,mask=mask)
                mask[-6]=0
                mask[-1]=1
                for l in range(10):
                    atoms.set_angle(i,total_atoms-6,total_atoms-1,109.5,mask=mask)
                    atoms.set_angle(total_atoms-7,total_atoms-6,total_atoms-1,109.5,mask=mask)
                atoms.set_distance(total_atoms-1,total_atoms-6,1.0100,fix=1,mask=mask)
                mask[-1]=0
        mask[-7+j]=0
    try:
        if strict[1]=='N':
            del atoms[-5]
        else:
            pass
    except:
        pass
    deleted=[]
    for j in range(2):
        try:
            ring_element=strict[j+2]
        except: 
            continue
        if ring_element=='N':
            if j+2==2:
                deleted.append(total_atoms-4)
            if j+2==3:
                deleted.append(total_atoms-2)
        else:
            if j+2==2:
                deleted.append(total_atoms-3)
                deleted.append(total_atoms-4)
            if j+2==3:
                deleted.append(total_atoms-1)
                deleted.append(total_atoms-2)
    deleted.sort(reverse=True)
    for d in deleted:
        del atoms[d]
    return atoms

def ring5(atoms,i,state=1,strict={1:'X'},shuffle=0):
    Cnear=near(i, atoms)
    try:
        if strict[1]=='N':
            pass
    except:
        strict[1]='X'
    if strict[1]=='N':
        atoms=NH2(atoms,i,state)
    else:
        atoms=CH3(atoms,i,state)
    for j in range(2):
        atoms.append('C')
    for j in range(8):
        atoms.append('H')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-10)+[1]*10
    elements=atoms.get_chemical_symbols()
    if shuffle==0:
        box=[-12,-11]
    elif shuffle==1:
        box=[-13,-12]
    elif shuffle==2:
        box=[-11,-13]
    for j in box:
        elements[total_atoms+j]='C'
        atoms.set_chemical_symbols(elements)
        if strict[1]=='N':
            atoms.set_distance(total_atoms+j,i,1.4480,fix=1)
        else:
            atoms.set_distance(total_atoms+j,i,1.5290,fix=1)
    x=atoms.get_positions()
    x[-10][0]=x[-10][0]-50
    x[-10][1]=x[-10][1]-50
    x[-10][2]=x[-10][2]-50
    x[-9][0]=x[-9][0]+50
    x[-9][1]=x[-9][1]+50
    x[-9][2]=x[-9][2]+50
    for j in range(-8,-4):
        x[j][0]=x[j][0]-50
        x[j][1]=x[j][1]-50
        x[j][2]=x[j][2]-50
    for j in range(-4,0):
        x[j][0]=x[j][0]+50
        x[j][1]=x[j][1]+50
        x[j][2]=x[j][2]+50
    if strict[1]=='N':
        atoms.set_positions(x)
        atoms.set_angle(i,total_atoms+box[0],total_atoms-10,104.98,mask=mask)
        atoms.set_distance(total_atoms-10,total_atoms+box[0],1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-10,164.01,mask=mask)
        mask[-10]=0
        atoms.set_angle(i,total_atoms+box[1],total_atoms-9,104.98,mask=mask)
        atoms.set_distance(total_atoms-9,total_atoms+box[1],1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-9,-164.01,mask=mask)
        mask[-9]=0
    else:
        atoms.set_positions(x)
        atoms.set_angle(i,total_atoms+box[0],total_atoms-10,99.37,mask=mask)
        atoms.set_distance(total_atoms-10,total_atoms+box[0],1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-10,-161,mask=mask)
        mask[-10]=0
        atoms.set_angle(i,total_atoms+box[1],total_atoms-9,99.37,mask=mask)
        atoms.set_distance(total_atoms-9,total_atoms+box[1],1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-9,161,mask=mask)
        mask[-9]=0
    atoms.set_angle(i,total_atoms+box[0],total_atoms-8,109.5,mask=mask)
    atoms.set_distance(total_atoms-8,total_atoms+box[0],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-8,75,mask=mask)
    mask[-8]=0
    atoms.set_angle(i,total_atoms+box[0],total_atoms-7,109.5,mask=mask)
    atoms.set_distance(total_atoms-7,total_atoms+box[0],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-7,-45,mask=mask)
    mask[-7]=0
    atoms.set_angle(i,total_atoms+box[1],total_atoms-6,109.5,mask=mask)
    atoms.set_distance(total_atoms-6,total_atoms+box[1],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-6,-75,mask=mask)
    mask[-6]=0
    atoms.set_angle(i,total_atoms+box[1],total_atoms-5,109.5,mask=mask)
    atoms.set_distance(total_atoms-5,total_atoms+box[1],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-5,45,mask=mask)
    mask[-5]=0
    if strict[1]=='N':
        atoms.set_angle(total_atoms+box[0],total_atoms-10,total_atoms-4,109.5,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-10,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0],total_atoms-10,total_atoms-4,-140.5,mask=mask)
        mask[-4]=0
        atoms.set_angle(total_atoms+box[0],total_atoms-10,total_atoms-3,109.5,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-10,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0],total_atoms-10,total_atoms-3,95,mask=mask)
        mask[-3]=0
        atoms.set_angle(total_atoms+box[1],total_atoms-9,total_atoms-2,109.5,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-9,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1],total_atoms-9,total_atoms-2,140.5,mask=mask)
        mask[-2]=0
        atoms.set_angle(total_atoms+box[1],total_atoms-9,total_atoms-1,109.5,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-9,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1],total_atoms-9,total_atoms-1,-95,mask=mask)
        mask[-1]=0
    else:
        atoms.set_angle(total_atoms+box[0],total_atoms-10,total_atoms-4,109.5,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-10,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0],total_atoms-10,total_atoms-4,-90.5,mask=mask)
        mask[-4]=0
        atoms.set_angle(total_atoms+box[0],total_atoms-10,total_atoms-3,109.5,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-10,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0],total_atoms-10,total_atoms-3,145,mask=mask)
        mask[-3]=0
        atoms.set_angle(total_atoms+box[1],total_atoms-9,total_atoms-2,109.5,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-9,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1],total_atoms-9,total_atoms-2,90.5,mask=mask)
        mask[-2]=0
        atoms.set_angle(total_atoms+box[1],total_atoms-9,total_atoms-1,109.5,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-9,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1],total_atoms-9,total_atoms-1,-145,mask=mask)
        mask[-1]=0
    l=0
    for j in range(4):
        try:
            ring_element=strict[j+2]
        except: 
            continue
        if ring_element=='O':
            mask.append(1)
            atoms.append(ring_element)
            if j+2==2:
                atoms.set_distance(total_atoms+l,total_atoms+box[0],0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-10,1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,i,1.4100,fix=1,mask=mask)
            elif j+2==3:
                atoms.set_distance(total_atoms+l,total_atoms-10,0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[0],1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-9,1.4100,fix=1,mask=mask)
            elif j+2==4:
                atoms.set_distance(total_atoms+l,total_atoms-9,0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[1],1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-10,1.4100,fix=1,mask=mask)
            elif j+2==5:
                atoms.set_distance(total_atoms+l,total_atoms+box[1],0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,i,1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-9,1.4100,fix=1,mask=mask)
            mask[-1]=0
            l=l+1
        elif ring_element=='S':
            mask.append(1)
            atoms.append(ring_element)
            if j+2==4:
                atoms.set_distance(total_atoms+l,total_atoms-9,0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[1],1.8100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-10,1.8100,fix=1,mask=mask)
            elif j+2==3:
                atoms.set_distance(total_atoms+l,total_atoms-10,0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[0],1.8100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-9,1.8100,fix=1,mask=mask)
            elif j+2==5:
                atoms.set_distance(total_atoms+l,total_atoms+box[1],0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,i,1.8100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-9,1.8100,fix=1,mask=mask)
            elif j+2==2:
                atoms.set_distance(total_atoms+l,total_atoms+box[0],0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-10,1.8100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,i,1.8100,fix=1,mask=mask)
            mask[-1]=0
            l=l+1
        elif ring_element=='N':
            mask.append(1)
            mask.append(1)
            atoms.append(ring_element)
            atoms.append('H')
            x=atoms.get_positions()
            x[-1][0]=x[-1][0]-50
            x[-1][1]=x[-1][1]-50
            x[-1][2]=x[-1][2]-50
            atoms.set_positions(x)
            if j+2==4:
                atoms.set_distance(total_atoms+l,total_atoms-9,0,fix=1,mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[1],1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-10,1.4480,fix=1,mask=mask)
                mask[-2]=0
                for k in range(10):
                    atoms.set_angle(total_atoms-10,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                    atoms.set_angle(total_atoms+box[1],total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                atoms.set_distance(total_atoms+l+1,total_atoms+l,1.0100,fix=1,mask=mask)
            elif j+2==3:
                atoms.set_distance(total_atoms+l,total_atoms-10,0,fix=1,mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[0],1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-9,1.4480,fix=1,mask=mask)
                mask[-2]=0
                for k in range(10):
                    atoms.set_angle(total_atoms+box[0],total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                    atoms.set_angle(total_atoms-9,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                atoms.set_distance(total_atoms+l+1,total_atoms+l,1.0100,fix=1,mask=mask)
            elif j+2==5:
                atoms.set_distance(total_atoms+l,total_atoms+box[1],0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,i,1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-9,1.4480,fix=1,mask=mask)
                mask[-2]=0
                for k in range(10):
                    atoms.set_angle(i,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                    atoms.set_angle(total_atoms-9,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                atoms.set_distance(total_atoms+l+1,total_atoms+l,1.0100,fix=1,mask=mask)
            elif j+2==2:
                atoms.set_distance(total_atoms+l,total_atoms+box[0],0,fix=1,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-10,1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,i,1.4480,fix=1,mask=mask)
                mask[-2]=0
                for k in range(10):
                    atoms.set_angle(i,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                    atoms.set_angle(total_atoms-10,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                atoms.set_distance(total_atoms+l+1,total_atoms+l,1.0100,fix=1,mask=mask)
            mask[-1]=0
            l=l+2
    deleted=[]
    for j in range(5):
        try:
            ring_element=strict[j+2]
        except: 
            continue
        if j+2==4:
            deleted.append(total_atoms-1)
            deleted.append(total_atoms-2)
            deleted.append(total_atoms-9)
        if j+2==3:
            deleted.append(total_atoms-3)
            deleted.append(total_atoms-4)
            deleted.append(total_atoms-10)
        if j+2==5:
            deleted.append(total_atoms-5)
            deleted.append(total_atoms-6)
            deleted.append(total_atoms+box[1])
        if j+2==2:
            deleted.append(total_atoms-7)
            deleted.append(total_atoms-8)
            deleted.append(total_atoms+box[0])
    deleted.sort(reverse=True)
    for d in deleted:
        del atoms[d]
    return atoms

def ring6(atoms,i,state=1,strict={1:'X'},shuffle=0):
    Cnear=near(i, atoms)
    try:
        if strict[1]=='N':
            pass
    except:
        strict[1]='X'
    if strict[1]=='N':
        atoms=NH2(atoms,i,state)
    else:
        atoms=CH3(atoms,i,state)
    for j in range(3):
        atoms.append('C')
    for j in range(10):
        atoms.append('H')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-13)+[1]*13
    elements=atoms.get_chemical_symbols()
    if shuffle==0:
        box=[-15,-14]
    elif shuffle==1:
        box=[-16,-15]
    elif shuffle==2:
        box=[-14,-16]
    for j in box:
        elements[total_atoms+j]='C'
        atoms.set_chemical_symbols(elements)
        if strict[1]=='N':
            atoms.set_distance(total_atoms+j,i,1.4480,fix=1)
        else:
            atoms.set_distance(total_atoms+j,i,1.5290,fix=1)
    x=atoms.get_positions()
    x[-13][0]=x[-13][0]-50
    x[-13][1]=x[-13][1]-50
    x[-13][2]=x[-13][2]-50
    x[-12][0]=x[-12][0]-50
    x[-12][1]=x[-12][1]-50
    x[-12][2]=x[-12][2]-50
    x[-11][0]=x[-11][0]+50
    x[-11][1]=x[-11][1]+50
    x[-11][2]=x[-11][2]+50
    for j in range(-10,-5):
        x[j][0]=x[j][0]-50
        x[j][1]=x[j][1]-50
        x[j][2]=x[j][2]-50
    for j in range(-5,0):
        x[j][0]=x[j][0]+50
        x[j][1]=x[j][1]+50
        x[j][2]=x[j][2]+50
    atoms.set_positions(x)
    if strict[1]=='N':
        atoms.set_angle(i,total_atoms+box[0],total_atoms-13,115.7,mask=mask)
        atoms.set_distance(total_atoms-13,total_atoms+box[0],1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-13,170,mask=mask)
        mask[-13]=0
        atoms.set_angle(i,total_atoms+box[1],total_atoms-12,115.7,mask=mask)
        atoms.set_distance(total_atoms-12,total_atoms+box[1],1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-12,-170,mask=mask)
        mask[-12]=0
        atoms.set_angle(total_atoms+box[0],total_atoms-13,total_atoms-11,113.7,mask=mask)
        for k in range(10):
            atoms.set_distance(total_atoms-11,total_atoms-12,1.5290,fix=1,mask=mask)
            atoms.set_distance(total_atoms-11,total_atoms-13,1.5290,fix=1,mask=mask)
        mask[-11]=0
    else:
        atoms.set_angle(i,total_atoms+box[0],total_atoms-13,107.4,mask=mask)
        atoms.set_distance(total_atoms-13,total_atoms+box[0],1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-13,170,mask=mask)
        mask[-13]=0
        atoms.set_angle(i,total_atoms+box[1],total_atoms-12,107.4,mask=mask)
        atoms.set_distance(total_atoms-12,total_atoms+box[1],1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-12,-170,mask=mask)
        mask[-12]=0
        atoms.set_angle(total_atoms+box[0],total_atoms-13,total_atoms-11,107.4,mask=mask)
        atoms.set_distance(total_atoms-11,total_atoms-13,1.5290,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0],total_atoms-13,total_atoms-11,55,mask=mask)
        mask[-11]=0
    atoms.set_angle(i,total_atoms+box[0],total_atoms-10,109.5,mask=mask)
    atoms.set_distance(total_atoms-10,total_atoms+box[0],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-10,-60,mask=mask)
    mask[-10]=0
    atoms.set_angle(i,total_atoms+box[0],total_atoms-9,109.5,mask=mask)
    atoms.set_distance(total_atoms-9,total_atoms+box[0],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-9,60,mask=mask)
    mask[-9]=0
    atoms.set_angle(i,total_atoms+box[1],total_atoms-8,109.5,mask=mask)
    atoms.set_distance(total_atoms-8,total_atoms+box[1],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-8,-60,mask=mask)
    mask[-8]=0
    atoms.set_angle(i,total_atoms+box[1],total_atoms-7,109.5,mask=mask)
    atoms.set_distance(total_atoms-7,total_atoms+box[1],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-7,60,mask=mask)
    mask[-7]=0
    if strict[1]=='N':
        atoms.set_angle(total_atoms+box[0],total_atoms-13,total_atoms-6,109.5,mask=mask)
        atoms.set_distance(total_atoms-6,total_atoms-13,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0],total_atoms-13,total_atoms-6,165,mask=mask)
        mask[-6]=0
        atoms.set_angle(total_atoms+box[0],total_atoms-13,total_atoms-5,109.5,mask=mask)
        atoms.set_distance(total_atoms-5,total_atoms-13,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0],total_atoms-13,total_atoms-5,59.5,mask=mask)
        mask[-5]=0
        atoms.set_angle(total_atoms+box[1],total_atoms-12,total_atoms-4,109.5,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-12,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1],total_atoms-12,total_atoms-4,-165,mask=mask)
        mask[-4]=0
        atoms.set_angle(total_atoms+box[1],total_atoms-12,total_atoms-3,109.5,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-12,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1],total_atoms-12,total_atoms-3,-59.5,mask=mask)
        mask[-3]=0
    else:
        atoms.set_angle(total_atoms+box[0],total_atoms-13,total_atoms-6,109.5,mask=mask)
        atoms.set_distance(total_atoms-6,total_atoms-13,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0],total_atoms-13,total_atoms-6,-75,mask=mask)
        mask[-6]=0
        atoms.set_angle(total_atoms+box[0],total_atoms-13,total_atoms-5,109.5,mask=mask)
        atoms.set_distance(total_atoms-5,total_atoms-13,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0],total_atoms-13,total_atoms-5,179.5,mask=mask)
        mask[-5]=0
        atoms.set_angle(total_atoms+box[1],total_atoms-12,total_atoms-4,109.5,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-12,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1],total_atoms-12,total_atoms-4,75,mask=mask)
        mask[-4]=0
        atoms.set_angle(total_atoms+box[1],total_atoms-12,total_atoms-3,109.5,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-12,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1],total_atoms-12,total_atoms-3,-179.5,mask=mask)
        mask[-3]=0
    atoms.set_angle(total_atoms-13,total_atoms-11,total_atoms-2,109.5,mask=mask)
    atoms.set_distance(total_atoms-2,total_atoms-11,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms+box[0],total_atoms-13,total_atoms-11,total_atoms-2,-179.5,mask=mask)
    mask[-2]=0
    atoms.set_angle(total_atoms-13,total_atoms-11,total_atoms-1,109.5,mask=mask)
    atoms.set_distance(total_atoms-1,total_atoms-11,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms+box[0],total_atoms-12,total_atoms-11,total_atoms-1,300,mask=mask)
    mask[-1]=0
    l=0
    for j in range(5):
        try:
            ring_element=strict[j+2]
        except: 
            continue
        if ring_element=='O':
            mask.append(1)
            atoms.append(ring_element)
            if j+2==4:
                atoms.set_distance(total_atoms+l,total_atoms-11,0.1,fix=1,mask=mask)
                atoms.set_angle(total_atoms+box[0],total_atoms-13,total_atoms+l,100,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms-13,1.4100,fix=1,mask=mask)
                atoms.set_dihedral(i,total_atoms+box[0],total_atoms-13,total_atoms+l,47.9,mask=mask)
            elif j+2==3:
                atoms.set_distance(total_atoms+l,total_atoms-13,0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms+box[1],5,mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[0],1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-11,1.4100,fix=1,mask=mask)
            elif j+2==5:
                atoms.set_distance(total_atoms+l,total_atoms-12,0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms+box[0],5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[1],1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-11,1.4100,fix=1,mask=mask)
            elif j+2==2:
                atoms.set_distance(total_atoms+l,total_atoms+box[0],0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms-12,5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-13,1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,i,1.4100,fix=1,mask=mask)
            elif j+2==6:
                atoms.set_distance(total_atoms+l,total_atoms+box[1],0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms-13,5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-12,1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,i,1.4100,fix=1,mask=mask)
            l=l+1
        elif ring_element=='S':
            mask.append(1)
            atoms.append(ring_element)
            if j+2==4:
                atoms.set_distance(total_atoms+l,total_atoms-11,0.1,fix=1,mask=mask)
                atoms.set_angle(total_atoms+box[0],total_atoms-13,total_atoms+l,127.35,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms-13,1.8100,fix=1,mask=mask)
                atoms.set_dihedral(i,total_atoms+box[0],total_atoms-13,total_atoms+l,47.9,mask=mask)
            elif j+2==3:
                atoms.set_distance(total_atoms+l,total_atoms-13,0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms+box[1],5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[0],1.8100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-11,1.8100,fix=1,mask=mask)
            elif j+2==5:
                atoms.set_distance(total_atoms+l,total_atoms-12,0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms+box[0],5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[1],1.8100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-11,1.8100,fix=1,mask=mask)
            elif j+2==2:
                atoms.set_distance(total_atoms+l,total_atoms+box[0],0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms-12,5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-13,1.8100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,i,1.8100,fix=1,mask=mask)
            elif j+2==6:
                atoms.set_distance(total_atoms+l,total_atoms+box[1],0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms-13,5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-12,1.8100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,i,1.8100,fix=1,mask=mask)
            l=l+1
        elif ring_element=='N':
            mask.append(1)
            mask.append(1)
            atoms.append(ring_element)
            atoms.append('H')
            x=atoms.get_positions()
            x[-1][0]=x[-1][0]-50
            x[-1][1]=x[-1][1]-50
            x[-1][2]=x[-1][2]-50
            atoms.set_positions(x)
            if j+2==4:
                atoms.set_distance(total_atoms+l,total_atoms-11,0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,i,5,mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-13,1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-12,1.4480,fix=1,mask=mask)
                mask[-2]=0
                atoms.set_angle(total_atoms-13,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                atoms.set_distance(total_atoms+l+1,total_atoms+l,1.0100,fix=1,mask=mask)
                atoms.set_dihedral(total_atoms+l+1,total_atoms+l,total_atoms-13,total_atoms+box[0],-150,mask=mask)
            elif j+2==3:
                atoms.set_distance(total_atoms+l,total_atoms-13,0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms-14,5,mask=mask)
                # Pseudo-Optimization
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[0],1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-11,1.4480,fix=1,mask=mask)
                mask[-2]=0
                for k in range(10):
                    atoms.set_angle(total_atoms+box[0],total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                    atoms.set_angle(total_atoms-11,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                atoms.set_distance(total_atoms+l+1,total_atoms+l,1.0100,fix=1,mask=mask)
            elif j+2==5:
                atoms.set_distance(total_atoms+l,total_atoms-12,0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms+box[0],5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms+box[1],1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,total_atoms-11,1.4480,fix=1,mask=mask)
                mask[-2]=0
                for k in range(10):
                    atoms.set_angle(total_atoms+box[0],total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                    atoms.set_angle(total_atoms-11,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                atoms.set_distance(total_atoms+l+1,total_atoms+l,1.0100,fix=1,mask=mask)
            elif j+2==2:
                atoms.set_distance(total_atoms+l,total_atoms+box[0],0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms-12,5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-13,1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,i,1.4480,fix=1,mask=mask)
                mask[-2]=0
                for k in range(10):
                    atoms.set_angle(i,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                    atoms.set_angle(total_atoms-13,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                atoms.set_distance(total_atoms+l+1,total_atoms+l,1.0100,fix=1,mask=mask)
            elif j+2==6:
                atoms.set_distance(total_atoms+l,total_atoms+box[1],0,fix=1,mask=mask)
                atoms.set_distance(total_atoms+l,total_atoms-13,5,mask=mask)
                for k in range(10):
                    atoms.set_distance(total_atoms+l,total_atoms-12,1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms+l,i,1.4480,fix=1,mask=mask)
                mask[-2]=0
                for k in range(10):
                    atoms.set_angle(i,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                    atoms.set_angle(total_atoms-12,total_atoms+l,total_atoms+l+1,109.5,mask=mask)
                atoms.set_distance(total_atoms+l+1,total_atoms+l,1.0100,fix=1,mask=mask)
            l=l+2
        mask[-1]=0
    deleted=[]
    for j in range(5):
        try:
            ring_element=strict[j+2]
        except: 
            continue
        if j+2==4:
            deleted.append(total_atoms-1)
            deleted.append(total_atoms-2)
            deleted.append(total_atoms-11)
        if j+2==3:
            deleted.append(total_atoms-5)
            deleted.append(total_atoms-6)
            deleted.append(total_atoms-13)
        if j+2==5:
            deleted.append(total_atoms-3)
            deleted.append(total_atoms-4)
            deleted.append(total_atoms-12)
        if j+2==2:
            deleted.append(total_atoms-9)
            deleted.append(total_atoms-10)
            deleted.append(total_atoms+box[0])
        if j+2==6:
            deleted.append(total_atoms-7)
            deleted.append(total_atoms-8)
            deleted.append(total_atoms+box[1])
    deleted.sort(reverse=True)
    for d in deleted:
        del atoms[d]
    return atoms
    
def SOMe(atoms,i,state=1):
    atoms.append('O')
    atoms.append('H')
    total_atoms=len(atoms)
    elements=atoms.get_chemical_symbols()
    elements[i]='S'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    mask=[0]*(total_atoms-2)+[1]*2
    x=atoms.get_positions()
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    if state==1:
        atoms.set_distance(i,Cnear,1.7900,fix=1)
    elif state==2 or state==6:
        atoms.set_distance(i,Cnear,1.7700,fix=1)
    atoms.set_distance(total_atoms-2,i,1.5300,fix=1,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-2,107.0,mask=mask)
    atoms.set_dihedral(CnearC,Cnear,i,total_atoms-2,65,mask=mask)
    mask[-2]=0
    atoms.set_angle(total_atoms-2,i,total_atoms-1,107.0,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-1,96.0,mask=mask)
    atoms.set_distance(total_atoms-1,i,1.7900,fix=1,mask=mask)
    atoms.set_dihedral(CnearC,Cnear,i,total_atoms-1,175,mask=mask)
    atoms=CH3(atoms,total_atoms-1,state=162)
    return atoms

def SO2_(addition,atoms,i,state=1):
    atoms.append('O')
    atoms.append('O')
    atoms.append('H')
    total_atoms=len(atoms)
    elements=atoms.get_chemical_symbols()
    elements[i]='S'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    CnearC=near(Cnear,atoms,'C')
    mask=[0]*(total_atoms-3)+[1]*3
    x=atoms.get_positions()
    x[-3][0]=x[-3][0]-50
    x[-3][1]=x[-3][1]-50
    x[-3][2]=x[-3][2]-50
    x[-2][0]=x[-2][0]+50
    x[-2][1]=x[-2][1]+50
    x[-2][2]=x[-2][2]+50
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    if state==1:
        atoms.set_distance(i,Cnear,1.7700,fix=1)
    elif state==2 or state==6:
        atoms.set_distance(i,Cnear,1.7600,fix=1)
    atoms.set_distance(total_atoms-3,i,1.5300,fix=1,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-3,108.9,mask=mask)
    atoms.set_dihedral(CnearC,Cnear,i,total_atoms-2,65,mask=mask)
    mask[-3]=0
    atoms.set_distance(total_atoms-2,i,1.5300,fix=1,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-2,108.9,mask=mask)
    atoms.set_dihedral(CnearC,Cnear,i,total_atoms-2,-70,mask=mask)
    mask[-2]=0
    atoms.set_angle(total_atoms-2,i,total_atoms-1,108.9,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-1,97.67,mask=mask)
    atoms.set_distance(total_atoms-1,i,1.7700,fix=1,mask=mask)
    atoms.set_dihedral(CnearC,Cnear,i,total_atoms-1,175,mask=mask)
    if addition=='CH3' or addition=='Me':
        atoms=CH3(atoms,total_atoms-1,state=164)
    if addition=='OH':
        atoms=OH(atoms,total_atoms-1,state=164)
    #print(atoms.get_positions())
    return atoms

def D(member,atoms,i,state=1,shuffle=0):
    Cnear=near(i, atoms)
    atoms=CH3(atoms,i,state)
    for j in range(3):
        atoms.append('C')
    for j in range(8):
        atoms.append('H')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-11)+[1]*11
    elements=atoms.get_chemical_symbols()
    if shuffle==0:
        box=[-13,-12]
    elif shuffle==1:
        box=[-14,-13]
    elif shuffle==2:
        box=[-12,-14]
    for j in box:
        mask[j]=1
        elements[total_atoms+j]='C'
        atoms.set_chemical_symbols(elements)
        atoms.set_distance(total_atoms+j,i,1.5290,fix=1)
        atoms.set_angle(Cnear,i,total_atoms+j,113,mask=mask)
        mask[j]=0
    x=atoms.get_positions()
    x[-11][0]=x[-11][0]-50
    x[-11][1]=x[-11][1]-50
    x[-11][2]=x[-11][2]-50
    x[-10][0]=x[-10][0]+50
    x[-10][1]=x[-10][1]+50
    x[-10][2]=x[-10][2]+50
    x[-9][0]=x[-9][0]-50
    x[-9][1]=x[-9][1]-50
    x[-9][2]=x[-9][2]-50
    for j in range(-8,-4):
        x[j][0]=x[j][0]-50
        x[j][1]=x[j][1]-50
        x[j][2]=x[j][2]-50
    for j in range(-4,0):
        x[j][0]=x[j][0]+50
        x[j][1]=x[j][1]+50
        x[j][2]=x[j][2]+50
    atoms.set_positions(x)
    atoms.set_angle(i,total_atoms+box[0],total_atoms-11,102.3,mask=mask)
    atoms.set_distance(total_atoms-11,total_atoms+box[0],1.5290,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-11,235,mask=mask)
    mask[-11]=0
    atoms.set_angle(i,total_atoms+box[1],total_atoms-10,115.5,mask=mask)
    atoms.set_distance(total_atoms-10,total_atoms+box[1],1.5290,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-10,165,mask=mask)
    mask[-10]=0
    atoms.set_angle(total_atoms+box[0],total_atoms-11,total_atoms-9,110.5,mask=mask)
    atoms.set_dihedral(i,total_atoms+box[0],total_atoms-11,total_atoms-9,290,mask=mask)
    for k in range(10):
        atoms.set_distance(total_atoms-9,total_atoms-11,1.5290,fix=1,mask=mask)
        atoms.set_distance(total_atoms-9,total_atoms-10,1.5290,fix=1,mask=mask)
    mask[-9]=0
    atoms.set_angle(i,total_atoms+box[0],total_atoms-8,109.5,mask=mask)
    atoms.set_distance(total_atoms-8,total_atoms+box[0],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-8,-60,mask=mask)
    mask[-8]=0
    atoms.set_angle(i,total_atoms+box[0],total_atoms-7,109.5,mask=mask)
    atoms.set_distance(total_atoms-7,total_atoms+box[0],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[0],total_atoms-7,60,mask=mask)
    mask[-7]=0
    atoms.set_angle(i,total_atoms+box[1],total_atoms-6,109.5,mask=mask)
    atoms.set_distance(total_atoms-6,total_atoms+box[1],1.0900,fix=1,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms+box[1],total_atoms-6,-60,mask=mask)
    mask[-6]=0
    atoms.set_angle(total_atoms+box[0],total_atoms-11,total_atoms-5,109.5,mask=mask)
    atoms.set_distance(total_atoms-5,total_atoms-11,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(i,total_atoms+box[0],total_atoms-11,total_atoms-5,179.5,mask=mask)
    mask[-5]=0
    atoms.set_angle(total_atoms+box[1],total_atoms-10,total_atoms-4,109.5,mask=mask)
    atoms.set_distance(total_atoms-4,total_atoms-10,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(i,total_atoms+box[1],total_atoms-10,total_atoms-4,75,mask=mask)
    mask[-4]=0
    atoms.set_angle(total_atoms+box[1],total_atoms-10,total_atoms-3,109.5,mask=mask)
    atoms.set_distance(total_atoms-3,total_atoms-10,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(i,total_atoms+box[1],total_atoms-10,total_atoms-3,-179.5,mask=mask)
    mask[-3]=0
    atoms.set_angle(total_atoms-11,total_atoms-9,total_atoms-2,109.5,mask=mask)
    atoms.set_distance(total_atoms-2,total_atoms-9,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms+box[0],total_atoms-11,total_atoms-9,total_atoms-2,-179.5,mask=mask)
    mask[-2]=0
    atoms.set_angle(total_atoms-11,total_atoms-9,total_atoms-1,109.5,mask=mask)
    atoms.set_distance(total_atoms-1,total_atoms-9,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms+box[0],total_atoms-11,total_atoms-9,total_atoms-1,300,mask=mask)
    mask[-1]=0
    if member==5:
        for j in range(1):
            atoms.append('C')
        for j in range(2):
            atoms.append('H')
        total_atoms=len(atoms)
        mask=[0]*(total_atoms-3)+[1]*3
        atoms.set_angle(i,total_atoms+box[1]-3,total_atoms-3,109.5,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms+box[1]-3,1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[1]-3,total_atoms-3,60,mask=mask)
        for k in range(10):
            atoms.set_distance(total_atoms-3,total_atoms+box[1]-3,1.5290,fix=1,mask=mask)
            atoms.set_distance(total_atoms-3,total_atoms-14,1.5290,fix=1,mask=mask)
        mask[-3]=0
        atoms.set_angle(total_atoms+box[1]-3,total_atoms-3,total_atoms-2,113.2,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-3,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1]-3,total_atoms-3,total_atoms-2,185,mask=mask)
        mask[-2]=0
        atoms.set_angle(total_atoms+box[1]-3,total_atoms-3,total_atoms-1,113.2,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-3,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1]-3,total_atoms-3,total_atoms-1,300,mask=mask)
    elif member==6:
        for j in range(2):
            atoms.append('C')
        for j in range(4):
            atoms.append('H')
        total_atoms=len(atoms)
        mask=[0]*(total_atoms-6)+[1]*6
        atoms.set_angle(i,total_atoms+box[1]-6,total_atoms-6,110.5,mask=mask)
        atoms.set_distance(total_atoms-6,total_atoms+box[1]-6,1.5290,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i,total_atoms+box[1]-6,total_atoms-6,67.5,mask=mask)
        mask[-6]=0
        atoms.set_angle(total_atoms+box[0]-6,total_atoms-17,total_atoms-5,110.5,mask=mask)
        atoms.set_distance(total_atoms-5,total_atoms-17,1.5290,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[0]-6,total_atoms-17,total_atoms-5,67.5,mask=mask)
        for k in range(10):
            atoms.set_distance(total_atoms-5,total_atoms-6,1.5290,fix=1,mask=mask)
            atoms.set_distance(total_atoms-5,total_atoms-17,1.5290,fix=1,mask=mask)
        mask[-5]=0
        atoms.set_angle(total_atoms+box[1]-6,total_atoms-6,total_atoms-4,110.2,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-6,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1]-6,total_atoms-6,total_atoms-4,185,mask=mask)
        mask[-4]=0
        atoms.set_angle(total_atoms+box[1]-6,total_atoms-6,total_atoms-3,110.2,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-6,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(i,total_atoms+box[1]-6,total_atoms-6,total_atoms-3,300,mask=mask)
        mask[-3]=0
        atoms.set_angle(total_atoms-17,total_atoms-5,total_atoms-2,110.2,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-5,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms+box[0]-6,total_atoms-17,total_atoms-5,total_atoms-2,185,mask=mask)
        mask[-2]=0
        atoms.set_angle(total_atoms-17,total_atoms-5,total_atoms-1,110.2,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-5,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms+box[0]-6,total_atoms-17,total_atoms-5,total_atoms-1,60,mask=mask)
    return atoms
    
def Bpin(atoms,i,state=1,shuffle=0):
    atoms.append('O')
    atoms.append('O')
    total_atoms=len(atoms)
    elements=atoms.get_chemical_symbols()
    elements[i]='B'
    atoms.set_chemical_symbols(elements)
    Cnear=near(i, atoms)
    mask=[0]*(total_atoms-2)+[1]*2
    CnearH=near(Cnear,atoms,'H')
    x=atoms.get_positions()
    x[-2][0]=x[-2][0]+50
    x[-2][1]=x[-2][1]+50
    x[-2][2]=x[-2][2]+50
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    if state==1:
        atoms.set_distance(i,Cnear,1.5950,fix=1)
    atoms.set_distance(total_atoms-2,i,1.4700,fix=1,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-2,120.0,mask=mask)
    if shuffle==0:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-2,20,mask=mask)
    elif shuffle==1:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-2,-20,mask=mask)
    elif shuffle==5:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-2,0.5,mask=mask)
    mask[-2]=0
    atoms.set_distance(total_atoms-1,i,1.4700,fix=1,mask=mask)
    atoms.set_angle(Cnear,i,total_atoms-1,120.0,mask=mask)
    if shuffle==0:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-1,200,mask=mask)
    elif shuffle==1:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-1,160,mask=mask)
    elif shuffle==5:
        atoms.set_dihedral(CnearH,Cnear,i,total_atoms-1,180.5,mask=mask)
    atoms.append('C')
    atoms.append('C')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-2)+[1]*2
    x=atoms.get_positions()
    x[-2][0]=x[-2][0]+50
    x[-2][1]=x[-2][1]+50
    x[-2][2]=x[-2][2]+50
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    atoms.set_distance(total_atoms-2,total_atoms-4,1.4100,fix=1,mask=mask)
    atoms.set_angle(i,total_atoms-4,total_atoms-2,95.0,mask=mask)
    atoms.set_dihedral(Cnear,i,total_atoms-4,total_atoms-2,125,mask=mask)
    mask[-2]=0
    atoms.set_dihedral(Cnear,i,total_atoms-3,total_atoms-1,235,mask=mask)
    for k in range(10):
        atoms.set_distance(total_atoms-1,total_atoms-3,1.4100,fix=1,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-2,1.5290,fix=1,mask=mask)
    for i in range(4):
        atoms.append('H')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-4)+[1]*4
    x=atoms.get_positions()
    x[-4][0]=x[-4][0]-50
    x[-4][1]=x[-4][1]-50
    x[-4][2]=x[-4][2]-50
    x[-3][0]=x[-3][0]+50
    x[-3][1]=x[-3][1]+50
    x[-3][2]=x[-3][2]+50
    x[-2][0]=x[-2][0]-50
    x[-2][1]=x[-2][1]-50
    x[-2][2]=x[-2][2]-50
    x[-1][0]=x[-1][0]+50
    x[-1][1]=x[-1][1]+50
    x[-1][2]=x[-1][2]+50
    atoms.set_dihedral(i,total_atoms-8,total_atoms-6,total_atoms-4,120,mask=mask)
    atoms.set_angle(total_atoms-8,total_atoms-6,total_atoms-4,111.5,mask=mask)
    atoms.set_distance(total_atoms-4,total_atoms-6,1.0900,fix=1)
    mask[-4]=0
    atoms.set_dihedral(i,total_atoms-8,total_atoms-6,total_atoms-3,240,mask=mask)
    atoms.set_angle(total_atoms-8,total_atoms-6,total_atoms-3,111.5,mask=mask)
    atoms.set_distance(total_atoms-3,total_atoms-6,1.0900,fix=1)
    mask[-3]=0
    atoms.set_dihedral(i,total_atoms-7,total_atoms-5,total_atoms-2,120,mask=mask)
    atoms.set_angle(total_atoms-7,total_atoms-5,total_atoms-2,111.5,mask=mask)
    atoms.set_distance(total_atoms-2,total_atoms-5,1.0900,fix=1)
    mask[-2]=0
    atoms.set_dihedral(i,total_atoms-7,total_atoms-5,total_atoms-1,240,mask=mask)
    atoms.set_angle(total_atoms-7,total_atoms-5,total_atoms-1,111.5,mask=mask)
    atoms.set_distance(total_atoms-1,total_atoms-5,1.0900,fix=1)
    if atoms.get_distance(total_atoms-4,total_atoms-5)<=1.3 or atoms.get_distance(total_atoms-3,total_atoms-5)<=1.3 or atoms.get_distance(total_atoms-2,total_atoms-6)<=1.3 or atoms.get_distance(total_atoms-1,total_atoms-6)<=1.3:
        mask=[0]*(total_atoms-4)+[1]*4
        x=atoms.get_positions()
        x[-4][0]=x[-4][0]-50
        x[-4][1]=x[-4][1]-50
        x[-4][2]=x[-4][2]-50
        x[-3][0]=x[-3][0]+50
        x[-3][1]=x[-3][1]+50
        x[-3][2]=x[-3][2]+50
        x[-2][0]=x[-2][0]-50
        x[-2][1]=x[-2][1]-50
        x[-2][2]=x[-2][2]-50
        x[-1][0]=x[-1][0]+50
        x[-1][1]=x[-1][1]+50
        x[-1][2]=x[-1][2]+50
        atoms.set_dihedral(i,total_atoms-8,total_atoms-6,total_atoms-4,300,mask=mask)
        atoms.set_angle(total_atoms-8,total_atoms-6,total_atoms-4,111.5,mask=mask)
        atoms.set_distance(total_atoms-4,total_atoms-6,1.0900,fix=1)
        mask[-4]=0
        atoms.set_dihedral(i,total_atoms-8,total_atoms-6,total_atoms-3,179.5,mask=mask)
        atoms.set_angle(total_atoms-8,total_atoms-6,total_atoms-3,111.5,mask=mask)
        atoms.set_distance(total_atoms-3,total_atoms-6,1.0900,fix=1)
        mask[-3]=0
        atoms.set_dihedral(i,total_atoms-7,total_atoms-5,total_atoms-2,60,mask=mask)
        atoms.set_angle(total_atoms-7,total_atoms-5,total_atoms-2,111.5,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-5,1.0900,fix=1)
        mask[-2]=0
        atoms.set_dihedral(i,total_atoms-7,total_atoms-5,total_atoms-1,179.5,mask=mask)
        atoms.set_angle(total_atoms-7,total_atoms-5,total_atoms-1,111.5,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-5,1.0900,fix=1)
    return atoms
    
def near(i,atoms,target='ALL',ignore=[],limit=True,reverse=False):
    elements=atoms.get_chemical_symbols()
    total_atoms=len(atoms)
    j=0
    targeted_distances=[]
    for element in elements:
        if j in ignore:
            j=j+1
            continue
        if element==target:
            d=atoms.get_distances(i,j)
            if limit==True:
                if d[0]>=0.5:
                    targeted_distances.append(d[0])
            else:
                targeted_distances.append(d[0])
        if target=='Not H' or target=='notH' or target=='not H':
            if element!='H':
                d=atoms.get_distances(i,j)
                if limit==True:
                    if d[0]>=0.5:
                        targeted_distances.append(d[0])
                else:
                    targeted_distances.append(d[0])
        j=j+1
    if reverse==True:
        _near=-1
    else:
        _near=0
    if target=='all' or target=='All' or target=='ALL':
        targeted_distances=atoms.get_distances(i,range(len(atoms)))
        targeted_distances=np.delete(targeted_distances,np.where(targeted_distances==0.0))
    distances=atoms.get_distances(i,range(len(atoms)))
    #print(distances,targeted_distances)
    for distance in distances:
        if total_atoms+_near in ignore:
            continue
        if reverse==True:
            if total_atoms+_near in ignore:
                _near=_near-1
                continue
            if abs(min(targeted_distances)-atoms.get_distances(i,total_atoms+_near))<=0.001:
                #print(_near,atoms.get_distances(i,_near))
                return total_atoms+_near
            _near=_near-1
        else:
            if _near in ignore:
                _near=_near+1
                continue
            if abs(min(targeted_distances)-atoms.get_distances(i,_near))<=0.001:
                #print(_near,atoms.get_distances(i,_near))
                return _near
            #print(min(targeted_distances),atoms.get_distances(i,_near))
            _near=_near+1

def bind_selection(ring_size,atoms,bond,state=[1,1],place=0,strict={1:'X'}, side=0, shuffle=0, less5=False, stable=False):
    '''This function will bind the ring onto a bond. 5 means adding a 5-membered-ring bound. 
    Adding a 6 in the front will make sure the conjugation of the ring.''' 
    
    if ring_size==2 or ring_size=='bind2':
        return bind2(atoms,bond,state)
    elif ring_size==3 or ring_size=='bind3':
        return bind3(atoms,bond,state,strict)
    elif ring_size==5 or ring_size=='bind5':
        return bind5(atoms,bond,state,strict,side,shuffle,less5,stable)
    elif ring_size==6 or ring_size=='bind6':
        return bind6(atoms,bond,state,strict,side,shuffle,less5)
    elif ring_size==65 or ring_size=='bindconj5':
        return bindconj5(atoms,bond,state,strict)
    elif ring_size==66 or ring_size=='bindconj6':
        return bindconj6(atoms,bond,state,strict)

def bind2(atoms,bond,state=[1,1]):
    i1=bond[0]
    i2=bond[1]
    elements=atoms.get_chemical_symbols()
    element1=elements[i1]
    element2=elements[i2]
    Cnear1=near(i1,atoms,'Not H',ignore=[i1,i2])
    Cnear2=near(i2,atoms,'Not H',ignore=[i1,i2],reverse=True)
    total_atoms=len(atoms)
    mask=[0]*total_atoms
    mask[i1]=1
    mask[i2]=1
    bondjudge2([Cnear1,Cnear2],atoms,[i1,i2],[element1,element2],state,mask)
    for j in range(2):
        CnearH=near(i1,atoms,'H')
        del atoms[CnearH]
        i1=numjudge(i1,CnearH)
        i2=numjudge(i2,CnearH)
    for j in range(2):
        CnearH=near(i2,atoms,'H')
        del atoms[CnearH]
        i1=numjudge(i1,CnearH)
        i2=numjudge(i2,CnearH)
    atoms.append('H')
    atoms.append('H')
    x=atoms.get_positions()
    x[-2][0]=x[-2][0]-100
    x[-2][1]=x[-2][1]-50
    x[-2][2]=x[-2][2]-50
    x[-1][0]=x[-1][0]+50
    x[-1][1]=x[-1][1]+50
    x[-1][2]=x[-1][2]+50
    atoms.set_positions(x)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-2)+[1]*2
    for k in range(10):
        atoms.set_angle(Cnear1,i1,total_atoms-2,118,mask=mask)
        atoms.set_angle(i2,i1,total_atoms-2,118,mask=mask)
    atoms.set_distance(total_atoms-2,i1,1.0800,fix=1,mask=mask)
    mask[-2]=0
    for k in range(10):
        atoms.set_angle(Cnear2,i2,total_atoms-1,118,mask=mask)
        atoms.set_angle(i1,i2,total_atoms-1,118,mask=mask)
    atoms.set_distance(total_atoms-1,i2,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-2,i1,i2,total_atoms-1,0.5,mask=mask)
    if element1=='N':
        del atoms[-2]
    if element2=='N':
        del atoms[-1]
    return atoms

def bind3(atoms,bond,state=[1,1],strict={1:'X'}):
    i1=bond[0]
    i2=bond[1]
    Cnear=near(i1,atoms,'C')
    total_atoms=len(atoms)
    if strict[1]=='O':
        atoms.append('O')
        mask=[0]*total_atoms+[1]
        if state==[1,1]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.4100,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.4100,fix=1,mask=mask)
        elif state==[1,2]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.4100,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.3640,fix=1,mask=mask)
        elif state==[2,1]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.3640,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.4100,fix=1,mask=mask)
        elif state==[2,2] or state==[6,6]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.3640,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.3640,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i1,i2,total_atoms,75,mask=mask)
        for j in range(2):
            CnearH=near(total_atoms-j,atoms,'H')
            del atoms[CnearH]
    elif strict[1]=='S':
        atoms.append('S')
        mask=[0]*total_atoms+[1]
        if state==[1,1]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.8100,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.8100,fix=1,mask=mask)
        elif state==[1,2]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.8100,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.7700,fix=1,mask=mask)
        elif state==[2,1]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.7700,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.8100,fix=1,mask=mask)
        elif state==[2,2] or state==[6,6]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.7700,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.7700,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i1,i2,total_atoms,75,mask=mask)
        for j in range(2):
            CnearH=near(total_atoms,atoms,'H')
            del atoms[CnearH]
            total_atoms=numjudge(total_atoms,CnearH)
            i1=numjudge(i1,CnearH)
            i2=numjudge(i2,CnearH)
    elif strict[1]=='N':
        atoms.append('N')
        atoms.append('H')
        mask=[0]*total_atoms+[1]*2
        if state==[1,1]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.4480,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.4480,fix=1,mask=mask)
        elif state==[1,2]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.4480,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.3400,fix=1,mask=mask)
        elif state==[2,1]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.3400,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.4480,fix=1,mask=mask)
        elif state==[2,2] or state==[6,6]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.3400,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.3400,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i1,i2,total_atoms,75,mask=mask)
        mask[-2]=0
        for j in range(2):
            CnearH=near(total_atoms,atoms,'H')
            del atoms[CnearH]
            total_atoms=numjudge(total_atoms,CnearH)
            i1=numjudge(i1,CnearH)
            i2=numjudge(i2,CnearH)
        x=atoms.get_positions()
        x[-1][0]=x[-1][0]+50
        x[-1][1]=x[-1][1]+50
        x[-1][2]=x[-1][2]+50
        atoms.set_positions(x)
        total_atoms=len(atoms)
        mask=[0]*(total_atoms-1)+[1]
        atoms.set_distance(total_atoms-1,total_atoms-2,0,fix=1,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-2,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-1,total_atoms-2,i1,i2,-90.5,mask=mask)
        for k in range(10):
            atoms.set_angle(i1,total_atoms-2,total_atoms-1,111.5,mask=mask)
            atoms.set_angle(i2,total_atoms-2,total_atoms-1,111.5,mask=mask)
    else:
        atoms.append('C')
        atoms.append('H')
        atoms.append('H')
        mask=[0]*total_atoms+[1]*3
        if state==[1,1]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.5290,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.5290,fix=1,mask=mask)
        elif state==[1,2]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.5290,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.5095,fix=1,mask=mask)
        elif state==[2,1]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.5095,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.5290,fix=1,mask=mask)
        elif state==[2,2]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.5095,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.5095,fix=1,mask=mask)
        elif state==[6,6]:
            for k in range(10):
                atoms.set_distance(total_atoms,i1,1.5105,fix=1,mask=mask)
                atoms.set_distance(total_atoms,i2,1.5105,fix=1,mask=mask)
        atoms.set_dihedral(Cnear,i1,i2,total_atoms,75,mask=mask)
        mask[-3]=0
        for j in range(2):
            CnearH=near(total_atoms,atoms,'H')
            del atoms[CnearH]
            total_atoms=numjudge(total_atoms,CnearH)
            i1=numjudge(i1,CnearH)
            i2=numjudge(i2,CnearH)
        x=atoms.get_positions()
        x[-2][0]=x[-2][0]+50
        x[-2][1]=x[-2][1]+50
        x[-2][2]=x[-2][2]+50
        x[-1][0]=x[-1][0]-50
        x[-1][1]=x[-1][1]-50
        x[-1][2]=x[-1][2]-50
        atoms.set_positions(x)
        total_atoms=len(atoms)
        mask=[0]*(total_atoms-2)+[1]*2
        atoms.set_distance(total_atoms-2,total_atoms-3,0,fix=1,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-3,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-2,total_atoms-3,i1,Cnear,110,mask=mask)
        for k in range(10):
            atoms.set_angle(i1,total_atoms-3,total_atoms-2,111.5,mask=mask)
            atoms.set_angle(i2,total_atoms-3,total_atoms-2,111.5,mask=mask)
        mask[-2]=0
        atoms.set_distance(total_atoms-1,total_atoms-3,0,fix=1,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-3,1.0900,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-1,total_atoms-3,i1,Cnear,-70.5,mask=mask)
        for k in range(10):
            atoms.set_angle(i1,total_atoms-3,total_atoms-1,111.5,mask=mask)
            atoms.set_angle(i2,total_atoms-3,total_atoms-1,111.5,mask=mask)
        atoms.set_dihedral(total_atoms-1,total_atoms-3,i1,i2,-90.5,mask=mask)
    return atoms

def bind5(atoms,bond,state=[1,1],strict={1:'X'},side=0,shuffle=0,less5=False,stable=False):
    i1=bond[0]
    i2=bond[1]
    if side==0:
        Hnear1=near(i1,atoms,'H')
        Hnear2=near(i2,atoms,'H')
    elif side==1:
        Hnear1=near(i1,atoms,'H',reverse=True)
        Hnear2=near(i2,atoms,'H',reverse=True)
    elif side==2:
        Hnear1=near(i1,atoms,'H')
        Hnear2=near(i2,atoms,'H',reverse=True)
    elif side==3:
        Hnear1=near(i1,atoms,'H',reverse=True)
        Hnear2=near(i2,atoms,'H')
    print(Hnear1,Hnear2)
    total_atoms=len(atoms)
    mask=[0]*total_atoms
    mask[Hnear1]=1
    atoms.set_angle(i2,i1,Hnear1,100,mask=mask)
    mask[Hnear1]=0
    mask[Hnear2]=1
    atoms.set_angle(i1,i2,Hnear2,-100,mask=mask)
    if (less5==True and side<2) or (less5==False and side>=2):
        atoms.set_dihedral(Hnear1,i1,i2,Hnear2,-35,mask=mask)
        xoy=near(Hnear2,atoms,limit=False)
        if atoms.get_distance(xoy,Hnear2)<=0.9:
            atoms.set_dihedral(Hnear1,i1,i2,Hnear2,35,mask=mask)
        xoy=near(i2,atoms,target='Not H')
        if xoy==i1:
            xoy=near(i2,atoms,target='Not H',reverse=True)
        if atoms.get_angle(xoy,i2,Hnear2)<=90:
            atoms.set_dihedral(Hnear1,i1,i2,Hnear2,35,mask=mask)
    elif less5==False and side<2:
        atoms.set_dihedral(Hnear1,i1,i2,Hnear2,55,mask=mask)
        xoy=near(Hnear2,atoms,limit=False)
        if atoms.get_distance(xoy,Hnear2)<=0.9:
            atoms.set_dihedral(Hnear1,i1,i2,Hnear2,-55,mask=mask)
        xoy=near(i2,atoms,target='Not H')
        if xoy==i1:
            xoy=near(i2,atoms,target='Not H',reverse=True)
        if atoms.get_angle(xoy,i2,Hnear2)<=90:
            atoms.set_dihedral(Hnear1,i1,i2,Hnear2,-55,mask=mask)
    elif less5==True and side>=2:
        atoms.set_dihedral(Hnear1,i1,i2,Hnear2,80,mask=mask)
        xoy=near(Hnear2,atoms,limit=False)
        if atoms.get_distance(xoy,Hnear2)<=0.9:
            atoms.set_dihedral(Hnear1,i1,i2,Hnear2,-80,mask=mask)
        xoy=near(i2,atoms,target='Not H')
        if xoy==i1:
            xoy=near(i2,atoms,target='Not H',reverse=True)
        if atoms.get_angle(xoy,i2,Hnear2)<=90:
            atoms.set_dihedral(Hnear1,i1,i2,Hnear2,-80,mask=mask)
        mask[Hnear2]=0
        mask[Hnear1]=1
        if atoms.get_dihedral(Hnear1,i1,i2,Hnear2)<85:
            atoms.set_dihedral(Hnear2,i2,i1,Hnear1,50,mask=mask)
        else:
            atoms.set_dihedral(Hnear2,i2,i1,Hnear1,360-50,mask=mask)
    if less5==True and side>=2:
        elements=atoms.get_chemical_symbols()
        elements[Hnear1]='C'
        elements[Hnear2]='C'
        atoms.set_chemical_symbols(elements)
        atoms.set_distance(Hnear1,i1,1.0400,fix=1)
        atoms.set_distance(Hnear2,i2,1.0400,fix=1)
        Hnear3=near(i1,atoms,'H')
        Hnear4=near(i2,atoms,'H')
        atoms.set_distance(Hnear1,i1,1.0900,fix=1)
        atoms.set_distance(Hnear2,i2,1.0900,fix=1)
        elements[Hnear1]='H'
        elements[Hnear2]='H'
        atoms.set_chemical_symbols(elements)
    #    x=atoms.get_positions()
    #    x[Hnear3][0]=x[Hnear3][0]-50
    #    x[Hnear3][1]=x[Hnear3][1]-50
    #    x[Hnear3][2]=x[Hnear3][2]-50
    #    x[Hnear4][0]=x[Hnear4][0]+50
    #    x[Hnear4][1]=x[Hnear4][1]+50
    #    x[Hnear4][2]=x[Hnear4][2]+50
    #    atoms.set_positions(x)
    if less5==True and side>=2:
        Cnear1=near(i1,atoms,'Not H')
        Cnear2=near(i2,atoms,'Not H')
        mask[Hnear1]=0
        mask[Hnear2]=0
        mask[Hnear3]=1
        mask[Hnear4]=1
        if Cnear1==i2:
            Cnear1=near(i1,atoms,'Not H',reverse=True)
        if Cnear2==i1:
            Cnear2=near(i2,atoms,'Not H',reverse=True)
        for k in range(10):
            atoms.set_angle(i2,i1,Hnear3,105.0,mask=mask)
            atoms.set_angle(Hnear1,i1,Hnear3,105.0,mask=mask)
        atoms.set_distance(Hnear3,i1,1.0900,fix=1,mask=mask)
        mask[Hnear3]=0
        for k in range(10):
            atoms.set_angle(i1,i2,Hnear4,105.0,mask=mask)
            atoms.set_angle(Hnear2,i2,Hnear4,105.0,mask=mask)
        atoms.set_distance(Hnear4,i2,1.0900,fix=1,mask=mask)
    if shuffle==0:
        atoms=CH3(atoms,Hnear1,state[0])
        atoms=CH3(atoms,Hnear2,state[1],ring=True)
    elif shuffle==1:
        atoms=CH3(atoms,Hnear1,state[0],ring=True)
        atoms=CH3(atoms,Hnear2,state[1])
    elif shuffle==2:
        atoms=CH3(atoms,Hnear1,state[0])
        atoms=CH3(atoms,Hnear2,state[1])
    elif shuffle==3:
        atoms=CH3(atoms,Hnear1,state[0],ring=True)
        atoms=CH3(atoms,Hnear2,state[1],ring=True)
    atoms.append('C')
    atoms.append('H')
    atoms.append('H')
    x=atoms.get_positions()
    x[-3][0]=x[-3][0]-50
    x[-3][1]=x[-3][1]-50
    x[-3][2]=x[-3][2]-50
    x[-2][0]=x[-2][0]+50
    x[-2][1]=x[-2][1]+50
    x[-2][2]=x[-2][2]+50
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    atoms.set_positions(x)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-3)+[1]*3
    for k in range(10):
        if atoms.get_distance(total_atoms-3,i1)<1.6 or stable==True:
            atoms.set_distance(total_atoms-3,i1,2.4,fix=1,mask=mask)
        if atoms.get_distance(total_atoms-3,i2)<1.6 or stable==True:
            atoms.set_distance(total_atoms-3,i2,2.4,fix=1,mask=mask)
        atoms.set_distance(total_atoms-3,Hnear1,1.5290,fix=1,mask=mask)
        atoms.set_distance(total_atoms-3,Hnear2,1.5290,fix=1,mask=mask)
    mask[-3]=0
    for j in range(2):
        total_atoms=len(atoms)
        CnearH=near(total_atoms-3,atoms,'H',limit=False)
        del atoms[CnearH]
        i1=numjudge(i1,CnearH)
        i2=numjudge(i2,CnearH)
        Hnear1=numjudge(Hnear1,CnearH)
        Hnear2=numjudge(Hnear2,CnearH)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-2)+[1]*2
    if side<2:
        atoms.set_distance(total_atoms-2,total_atoms-3,0,fix=1,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-3,1.0900,fix=1,mask=mask)
        while(1):
            if atoms.get_distance(total_atoms-2,Hnear1)<1.3:
                atoms.set_distance(total_atoms-2,Hnear1,1.8,fix=1,mask=mask)
            if atoms.get_distance(total_atoms-2,Hnear2)<1.3:
                atoms.set_distance(total_atoms-2,Hnear2,1.8,fix=1,mask=mask)
            xoy=near(total_atoms-1,atoms,limit=False)
            #if atoms.get_distance(xoy,total_atoms-2)<=0.85:
            #    atoms.set_distance(xoy,total_atoms-2,1,mask=mask)
            atoms.set_angle(Hnear1,total_atoms-3,total_atoms-2,111.0,mask=mask)
            atoms.set_angle(Hnear2,total_atoms-3,total_atoms-2,111.0,mask=mask)
            atoms.set_distance(total_atoms-2,total_atoms-3,1.0900,fix=1,mask=mask)
            if atoms.get_distance(total_atoms-2,Hnear1)>1.4 and atoms.get_distance(total_atoms-2,Hnear2)>1.4:
                if abs(atoms.get_angle(Hnear1,total_atoms-3,total_atoms-2)-atoms.get_angle(Hnear2,total_atoms-3,total_atoms-2))<1:
                    #print(atoms.get_angle(Hnear1,total_atoms-3,total_atoms-2)-atoms.get_angle(Hnear2,total_atoms-3,total_atoms-2))
                    break
        #atoms.set_dihedral(total_atoms-2,total_atoms-3,Hnear1,i1,60,mask=mask)
        mask[-2]=0
        #atoms.set_distance(total_atoms-1,total_atoms-3,0,fix=1,mask=mask)
        #atoms.set_distance(total_atoms-1,total_atoms-3,1.0900,fix=1,mask=mask)
        #while(1):
        #    if atoms.get_distance(total_atoms-1,Hnear1)<1.3:
        #        atoms.set_distance(total_atoms-1,Hnear1,1.8,fix=1,mask=mask)
        #    if atoms.get_distance(total_atoms-1,Hnear2)<1.3:
        #        atoms.set_distance(total_atoms-1,Hnear2,1.8,fix=1,mask=mask)
        #    atoms.set_angle(Hnear1,total_atoms-3,total_atoms-1,111.0,mask=mask)
        #    atoms.set_angle(Hnear2,total_atoms-3,total_atoms-1,111.0,mask=mask)
        #    atoms.set_distance(total_atoms-1,total_atoms-3,1.0900,fix=1,mask=mask)
        #    if atoms.get_distance(total_atoms-1,Hnear1)>1.4 and atoms.get_distance(total_atoms-1,Hnear2)>1.4:
        #        if abs(atoms.get_angle(Hnear1,total_atoms-3,total_atoms-1)-atoms.get_angle(Hnear2,total_atoms-3,total_atoms-1))<1:
        #            print(atoms.get_angle(Hnear1,total_atoms-3,total_atoms-1)-atoms.get_angle(Hnear2,total_atoms-3,total_atoms-1))
        #            break
        #atoms.set_dihedral(total_atoms-1,total_atoms-3,Hnear1,i1,170,mask=mask)
        #view(atoms)
        #mask[-2]=0
        #mask[-1]=1
        atoms.set_distance(total_atoms-1,total_atoms-2,0,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-1,total_atoms-3,Hnear1,i1,atoms.get_dihedral(total_atoms-1,total_atoms-3,Hnear1,i1)+110,mask=mask)
        xoy=near(total_atoms-1,atoms,limit=False)
        if atoms.get_distance(xoy,total_atoms-1)<=0.9:
            atoms.set_dihedral(total_atoms-1,total_atoms-3,Hnear1,i1,atoms.get_dihedral(total_atoms-1,total_atoms-3,Hnear1,i1)-220,mask=mask)
    else:
        atoms.set_distance(total_atoms-2,total_atoms-3,0,fix=1,mask=mask)
        atoms.set_distance(total_atoms-2,total_atoms-3,1.0900,fix=1,mask=mask)
        while(1):
            if atoms.get_distance(total_atoms-2,Hnear1)<1.3:
                atoms.set_distance(total_atoms-2,Hnear1,1.8,fix=1,mask=mask)
            if atoms.get_distance(total_atoms-2,Hnear2)<1.3:
                atoms.set_distance(total_atoms-2,Hnear2,1.8,fix=1,mask=mask)
            atoms.set_angle(Hnear1,total_atoms-3,total_atoms-2,111.0,mask=mask)
            atoms.set_angle(Hnear2,total_atoms-3,total_atoms-2,111.0,mask=mask)
            atoms.set_distance(total_atoms-2,total_atoms-3,1.0900,fix=1,mask=mask)
            if atoms.get_distance(total_atoms-2,Hnear1)>1.4 and atoms.get_distance(total_atoms-2,Hnear2)>1.4:
                break
        mask[-2]=0
        atoms.set_distance(total_atoms-1,total_atoms-2,0,fix=1,mask=mask)
        #atoms.set_distance(total_atoms-1,total_atoms-3,1.0900,fix=1,mask=mask)
        #while(1):
        #    if atoms.get_distance(total_atoms-1,Hnear1)<1.3:
        #        atoms.set_distance(total_atoms-1,Hnear1,1.8,fix=1,mask=mask)
        #    if atoms.get_distance(total_atoms-1,Hnear2)<1.3:
        #        atoms.set_distance(total_atoms-1,Hnear2,1.8,fix=1,mask=mask)
        #    atoms.set_angle(Hnear1,total_atoms-3,total_atoms-1,111.0,mask=mask)
        #    atoms.set_angle(Hnear2,total_atoms-3,total_atoms-1,111.0,mask=mask)
        #    atoms.set_distance(total_atoms-1,total_atoms-3,1.0900,fix=1,mask=mask)
        #    if atoms.get_distance(total_atoms-1,4)>1.4 and atoms.get_distance(total_atoms-1,6)>1.4:
        #        break
        atoms.set_dihedral(total_atoms-1,total_atoms-3,Hnear1,i1,atoms.get_dihedral(total_atoms-1,total_atoms-3,Hnear1,i1)+110,mask=mask)
        xoy=near(total_atoms-1,atoms,limit=False)
        if atoms.get_distance(xoy,total_atoms-1)<=0.9:
            atoms.set_dihedral(total_atoms-1,total_atoms-3,Hnear1,i1,atoms.get_dihedral(total_atoms-1,total_atoms-3,Hnear1,i1)-220,mask=mask)
    needelements=[]
    places=[]
    for j in range(3):
        try:
            needelements.append(strict[j+2])
            places.append(j+2)
        except:
            pass
    for i in range(len(needelements)):
        total_atoms=len(atoms)
        mask=[0]*total_atoms
        if places[i]==2:
            mask[Hnear1]=1
            elements=atoms.get_chemical_symbols()
            elements[Hnear1]=needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i]=='N':
                for k in range(10):
                    if state[0]==1:
                        atoms.set_distance(Hnear1,i1,1.4480,fix=1,mask=mask)
                    elif state[0]==2 or state[0]==6:
                        atoms.set_distance(Hnear1,i1,1.3400,fix=1,mask=mask)
                    atoms.set_distance(Hnear1,total_atoms-3,1.4480,fix=1,mask=mask)
                mask[Hnear1]=0
                total_atoms=len(atoms)
                CnearH=near(Hnear1,atoms,'H',limit=False)
                del atoms[CnearH]
                Hnear1=numjudge(Hnear1,CnearH)
                i1=numjudge(i1,CnearH)
                i2=numjudge(i2,CnearH)
                CnearH=near(Hnear1,atoms,'H',limit=False)
                mask[CnearH]=1
                for k in range(10):
                    atoms.set_angle(i1,Hnear1,CnearH,109.5,mask=mask)
                    atoms.set_angle(total_atoms-3,Hnear1,CnearH,109.5,mask=mask)
                atoms.set_dihedral(i2,i1,Hnear1,CnearH,-150,mask=mask)
                atoms.set_distance(CnearH,Hnear1,1.0100,fix=1)
            elif needelements[i]=='O':
                for k in range(10):
                    if state[0]==1:
                        atoms.set_distance(Hnear1,i1,1.4100,fix=1,mask=mask)
                    elif state[0]==2 or state[0]==6:
                        atoms.set_distance(Hnear1,i1,1.3640,fix=1,mask=mask)
                    atoms.set_distance(Hnear1,total_atoms-3,1.4100,fix=1,mask=mask)
                for j in range(2):
                    total_atoms=len(atoms)
                    CnearH=near(Hnear1,atoms,'H',limit=False)
                    del atoms[CnearH]
                    Hnear1=numjudge(Hnear1,CnearH)
        if places[i]==3:
            mask[total_atoms-3]=1
            elements=atoms.get_chemical_symbols()
            elements[total_atoms-3]=needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i]=='N':
                for k in range(10):
                    atoms.set_distance(total_atoms-3,Hnear1,1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms-3,Hnear2,1.4480,fix=1,mask=mask)
                mask[total_atoms-3]=0
                total_atoms=len(atoms)
                CnearH=near(total_atoms-3,atoms,'H',limit=False)
                del atoms[CnearH]
                Hnear1=numjudge(Hnear1,CnearH)
                i1=numjudge(i1,CnearH)
                i2=numjudge(i2,CnearH)
                CnearH=near(total_atoms-3,atoms,'H',limit=False)
                mask[CnearH]=1
                for k in range(10):
                    atoms.set_angle(Hnear1,total_atoms-3,CnearH,109.5,mask=mask)
                    atoms.set_angle(Hnear2,total_atoms-3,CnearH,109.5,mask=mask)
                atoms.set_dihedral(i1,Hnear1,total_atoms-3,CnearH,-150,mask=mask)
                atoms.set_distance(CnearH,total_atoms-3,1.0100,fix=1)
            elif needelements[i]=='O':
                for k in range(10):
                    atoms.set_distance(total_atoms-3,Hnear1,1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms-3,Hnear2,1.4100,fix=1,mask=mask)
                for j in range(2):
                    total_atoms=len(atoms)
                    CnearH=near(total_atoms-3+j,atoms,'H',limit=False)
                    del atoms[CnearH]
        if places[i]==4:
            mask[Hnear2]=1
            elements=atoms.get_chemical_symbols()
            elements[Hnear2]=needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i]=='N':
                for k in range(10):
                    if state[0]==1:
                        atoms.set_distance(Hnear2,i2,1.4480,fix=1,mask=mask)
                    elif state[0]==2 or state[0]==6:
                        atoms.set_distance(Hnear2,i2,1.3400,fix=1,mask=mask)
                    atoms.set_distance(Hnear2,total_atoms-3,1.4480,fix=1,mask=mask)
                mask[Hnear2]=0
                total_atoms=len(atoms)
                CnearH=near(Hnear2,atoms,'H',limit=False)
                del atoms[CnearH]
                Hnear1=numjudge(Hnear1,CnearH)
                i1=numjudge(i1,CnearH)
                i2=numjudge(i2,CnearH)
                CnearH=near(Hnear2,atoms,'H',limit=False)
                mask[CnearH]=1
                for k in range(10):
                    atoms.set_angle(i2,Hnear2,CnearH,109.5,mask=mask)
                    atoms.set_angle(total_atoms-3,Hnear2,CnearH,109.5,mask=mask)
                atoms.set_dihedral(i1,i2,Hnear2,CnearH,-150,mask=mask)
                atoms.set_distance(CnearH,Hnear2,1.0100,fix=1)
            elif needelements[i]=='O':
                for k in range(10):
                    if state[0]==1:
                        atoms.set_distance(Hnear2,i2,1.4100,fix=1,mask=mask)
                    elif state[0]==2 or state[0]==6:
                        atoms.set_distance(Hnear2,i2,1.3640,fix=1,mask=mask)
                    atoms.set_distance(Hnear2,total_atoms-3,1.4100,fix=1,mask=mask)
                for j in range(2):
                    total_atoms=len(atoms)
                    CnearH=near(Hnear2,atoms,'H',limit=False)
                    del atoms[CnearH]
                    Hnear2=numjudge(Hnear2,CnearH)
    return atoms

def bind6(atoms,bond,state=[1,1],strict={1:'X'},side=0,shuffle=0,less5=False):
    i1=bond[0]
    i2=bond[1]
    if side==0:
        Hnear1=near(i1,atoms,'H')
        Hnear2=near(i2,atoms,'H')
    elif side==1:
        Hnear1=near(i1,atoms,'H',reverse=True)
        Hnear2=near(i2,atoms,'H',reverse=True)
    elif side==2:
        Hnear1=near(i1,atoms,'H')
        Hnear2=near(i2,atoms,'H',reverse=True)
    elif side==3:
        Hnear1=near(i1,atoms,'H',reverse=True)
        Hnear2=near(i2,atoms,'H')
    print(Hnear1,Hnear2)
    if less5==True and side>=2:
        total_atoms=len(atoms)
        mask=[0]*total_atoms
        mask[Hnear1]=1
        mask[Hnear2]=1
        atoms.set_angle(i2,i1,Hnear1,107.8,mask=mask)
        mask[Hnear1]=0
        atoms.set_angle(i1,i2,Hnear2,107.8,mask=mask)
        atoms.set_dihedral(Hnear1,i1,i2,Hnear2,55,mask=mask)
        xoy=near(Hnear2,atoms,limit=False)
        if atoms.get_distance(xoy,Hnear2)<=0.9:
            atoms.set_dihedral(Hnear1,i1,i2,Hnear2,-55,mask=mask)
        xoy=near(i2,atoms,target='Not H')
        if xoy==i1:
            xoy=near(i2,atoms,target='Not H',reverse=True)
        if atoms.get_angle(xoy,i2,Hnear2)<=90:
            atoms.set_dihedral(Hnear1,i1,i2,Hnear2,-55,mask=mask)
    if shuffle==0:
        atoms=CH3(atoms,Hnear1,state[0])
        atoms=CH3(atoms,Hnear2,state[1])
    elif shuffle==1:
        atoms=CH3(atoms,Hnear1,state[0],ring=True)
        atoms=CH3(atoms,Hnear2,state[1],ring=True)
    elif shuffle==2:
        atoms=CH3(atoms,Hnear1,state[0],ring=True)
        atoms=CH3(atoms,Hnear2,state[1])
    elif shuffle==3:
        atoms=CH3(atoms,Hnear1,state[0])
        atoms=CH3(atoms,Hnear2,state[1],ring=True)
    atoms.append('C')
    atoms.append('C')
    for j in range(4):
        atoms.append('H')
    x=atoms.get_positions()
    x[-5][0]=x[-5][0]+50
    x[-5][1]=x[-5][1]+50
    x[-5][2]=x[-5][2]+50
    x[-6][0]=x[-6][0]-50
    x[-6][1]=x[-6][1]-50
    x[-6][2]=x[-6][2]-50
    for i in range(-2,0):
        x[2*i][0]=x[2*i][0]-50
        x[2*i][1]=x[2*i][1]-50
        x[2*i][2]=x[2*i][2]-50
        x[2*i+1][0]=x[2*i+1][0]+50
        x[2*i+1][1]=x[2*i+1][1]+50
        x[2*i+1][2]=x[2*i+1][2]+50
    atoms.set_positions(x)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-6)+[1]*6
    atoms.set_angle(i1,Hnear1,total_atoms-6,107.8,mask=mask)
    if state[0]==1:
        atoms.set_distance(total_atoms-6,Hnear1,1.5290,fix=1,mask=mask)
    elif state[0]==2:
        atoms.set_distance(total_atoms-6,Hnear1,1.5095,fix=1,mask=mask)
    elif state[0]==6:
        atoms.set_distance(total_atoms-6,Hnear1,1.5105,fix=1,mask=mask)
    if less5==False:
        atoms.set_dihedral(i2,i1,Hnear1,total_atoms-6,47.9,mask=mask)
        if atoms.get_distance(total_atoms-6,Hnear2)>=2.7:
            atoms.set_dihedral(i2,i1,Hnear1,total_atoms-6,-47.9,mask=mask)
    else:
        atoms.set_dihedral(i2,i1,Hnear1,total_atoms-6,-47.9,mask=mask)
        if atoms.get_distance(total_atoms-6,Hnear2)>=2.7:
            atoms.set_dihedral(i2,i1,Hnear1,total_atoms-6,47.9,mask=mask)
            
    mask[-6]=0
    atoms.set_angle(i2,Hnear2,total_atoms-5,107.8,mask=mask)
    for k in range(10):
        #print(atoms.get_dihedral(i1,i2,Hnear2,total_atoms-5))
        if less5==True and atoms.get_dihedral(Hnear2,i2,i1,Hnear1)>10 and atoms.get_dihedral(Hnear2,i2,i1,Hnear1)<350:
            if atoms.get_dihedral(i1,i2,Hnear2,total_atoms-5)<=30 and atoms.get_dihedral(i1,i2,Hnear2,total_atoms-5)>=0:
                atoms.set_dihedral(i1,i2,Hnear2,total_atoms-5,47.9,mask=mask)
            if atoms.get_dihedral(i1,i2,Hnear2,total_atoms-5)>=330 and atoms.get_dihedral(i1,i2,Hnear2,total_atoms-5)<=360:
                atoms.set_dihedral(i1,i2,Hnear2,total_atoms-5,-47.9,mask=mask)
        else:
            if atoms.get_dihedral(i1,i2,Hnear2,total_atoms-5)<=30 and atoms.get_dihedral(i1,i2,Hnear2,total_atoms-5)>=0:
                atoms.set_dihedral(i1,i2,Hnear2,total_atoms-5,-47.9,mask=mask)
            if atoms.get_dihedral(i1,i2,Hnear2,total_atoms-5)>=330 and atoms.get_dihedral(i1,i2,Hnear2,total_atoms-5)<=360:
                atoms.set_dihedral(i1,i2,Hnear2,total_atoms-5,47.9,mask=mask)
        atoms.set_distance(total_atoms-5,total_atoms-6,1.5290,fix=1,mask=mask)
        if state[1]==1:
            atoms.set_distance(total_atoms-5,Hnear2,1.5290,fix=1,mask=mask)
        elif state[1]==2:
            atoms.set_distance(total_atoms-5,Hnear2,1.5095,fix=1,mask=mask)
        elif state[1]==6:
            atoms.set_distance(total_atoms-5,Hnear2,1.5105,fix=1,mask=mask)
    mask[-5]=0
    for j in range(2):
        total_atoms=len(atoms)
        CnearH=near(total_atoms-6+j,atoms,'H',limit=False)
        del atoms[CnearH]
        i1=numjudge(i1,CnearH)
        i2=numjudge(i2,CnearH)
        Hnear1=numjudge(Hnear1,CnearH)
        Hnear2=numjudge(Hnear2,CnearH)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-4)+[1]*4
    atoms.set_distance(total_atoms-4,total_atoms-6,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-4,total_atoms-6,Hnear1,i1,-110,mask=mask)
    for k in range(10):
        atoms.set_angle(Hnear1,total_atoms-6,total_atoms-4,111.5,mask=mask)
        atoms.set_angle(total_atoms-5,total_atoms-6,total_atoms-4,111.5,mask=mask)
    mask[-4]=0
    atoms.set_distance(total_atoms-3,total_atoms-4,0,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-3,total_atoms-6,Hnear1,i1,atoms.get_dihedral(total_atoms-3,total_atoms-6,Hnear1,i1)-110,mask=mask)
    xoy=near(total_atoms-3,atoms,limit=False)
    if atoms.get_distance(xoy,total_atoms-3)<=0.9:
        atoms.set_dihedral(total_atoms-3,total_atoms-6,Hnear1,i1,atoms.get_dihedral(total_atoms-3,total_atoms-6,Hnear1,i1)+220,mask=mask)
    mask[-3]=0
    atoms.set_distance(total_atoms-2,total_atoms-5,0,fix=1,mask=mask)
    atoms.set_distance(total_atoms-2,total_atoms-5,1.0900,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-2,total_atoms-5,Hnear2,i2,110,mask=mask)
    for k in range(10):
        atoms.set_angle(Hnear2,total_atoms-5,total_atoms-2,111.5,mask=mask)
        atoms.set_angle(total_atoms-6,total_atoms-5,total_atoms-2,111.5,mask=mask)
    mask[-2]=0
    atoms.set_distance(total_atoms-1,total_atoms-2,0,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-1,total_atoms-5,Hnear2,i2,atoms.get_dihedral(total_atoms-1,total_atoms-5,Hnear2,i2)-110,mask=mask)
    xoy=near(total_atoms-1,atoms,limit=False)
    if atoms.get_distance(xoy,total_atoms-1)<=0.9:
        atoms.set_dihedral(total_atoms-1,total_atoms-5,Hnear2,i2,atoms.get_dihedral(total_atoms-1,total_atoms-5,Hnear2,i2)+220,mask=mask)

#    for k in range(10):
#        atoms.set_angle(Hnear2,total_atoms-5,total_atoms-1,111.5,mask=mask)
#        atoms.set_angle(total_atoms-6,total_atoms-5,total_atoms-1,111.5,mask=mask)
        
    mask[-1]=0
    needelements=[]
    places=[]
    for j in range(4):
        try:
            needelements.append(strict[j+2])
            places.append(j+2)
        except:
            pass
    error=0
    for i in range(len(needelements)):
        total_atoms=len(atoms)
        mask=[0]*total_atoms
        if places[i]==2:
            mask[Hnear1]=1
            elements=atoms.get_chemical_symbols()
            elements[Hnear1]=needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i]=='N':
                for k in range(10):
                    if state[0]==1:
                        atoms.set_distance(Hnear1,i1,1.4480,fix=1,mask=mask)
                    elif state[0]==2 or state[0]==6:
                        atoms.set_distance(Hnear1,i1,1.3400,fix=1,mask=mask)
                    atoms.set_distance(Hnear1,total_atoms-6,1.4480,fix=1,mask=mask)
                mask[Hnear1]=0
                total_atoms=len(atoms)
                CnearH=near(Hnear1,atoms,'H',limit=False)
                del atoms[CnearH]
                Hnear1=numjudge(Hnear1,CnearH)
                i1=numjudge(i1,CnearH)
                i2=numjudge(i2,CnearH)
                CnearH=near(Hnear1,atoms,'H',limit=False)
                mask[CnearH]=1
                for k in range(10):
                    atoms.set_angle(i1,Hnear1,CnearH,109.5,mask=mask)
                    atoms.set_angle(total_atoms-6,Hnear1,CnearH,109.5,mask=mask)
                atoms.set_dihedral(i2,i1,Hnear1,CnearH,-150,mask=mask)
                atoms.set_distance(CnearH,Hnear1,1.0100,fix=1)
            elif needelements[i]=='O':
                for k in range(10):
                    if state[0]==1:
                        atoms.set_distance(Hnear1,i1,1.4100,fix=1,mask=mask)
                    elif state[0]==2 or state[0]==6:
                        atoms.set_distance(Hnear1,i1,1.3640,fix=1,mask=mask)
                    atoms.set_distance(Hnear1,total_atoms-6,1.4100,fix=1,mask=mask)
                for j in range(2):
                    total_atoms=len(atoms)
                    CnearH=near(Hnear1,atoms,'H',limit=False)
                    del atoms[CnearH]
                    Hnear1=numjudge(Hnear1,CnearH)
        if places[i]==3:
            mask[total_atoms-6+error]=1
            elements=atoms.get_chemical_symbols()
            elements[total_atoms-6+error]=needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i]=='N':
                for k in range(10):
                    atoms.set_distance(total_atoms-6+error,Hnear1,1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms-6+error,total_atoms-5+error,1.4480,fix=1,mask=mask)
                mask[total_atoms-6+error]=0
                total_atoms=len(atoms)
                CnearH=near(total_atoms-6+error,atoms,'H',limit=False)
                del atoms[CnearH]
                Hnear1=numjudge(Hnear1,CnearH)
                i1=numjudge(i1,CnearH)
                i2=numjudge(i2,CnearH)
                CnearH=near(total_atoms-6+error,atoms,'H',limit=False)
                mask[CnearH]=1
                atoms.set_dihedral(i1,Hnear1,total_atoms-6+error,CnearH,-150,mask=mask)
                for k in range(10):
                    atoms.set_angle(Hnear1,total_atoms-6+error,CnearH,109.5,mask=mask)
                    atoms.set_angle(total_atoms-5+error,total_atoms-6+error,CnearH,109.5,mask=mask)
                atoms.set_distance(CnearH,total_atoms-6+error,1.0100,fix=1)
                error=error+1
            elif needelements[i]=='O':
                for k in range(10):
                    atoms.set_distance(total_atoms-6+error,Hnear1,1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms-6+error,total_atoms-5-error,1.4100,fix=1,mask=mask)
                for j in range(2):
                    total_atoms=len(atoms)
                    CnearH=near(total_atoms-6+error+j,atoms,'H',limit=False)
                    del atoms[CnearH]
                error=error+2
        if places[i]==4:
            mask[total_atoms-5+error]=1
            elements=atoms.get_chemical_symbols()
            elements[total_atoms-5-error]=needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i]=='N':
                for k in range(10):
                    atoms.set_distance(total_atoms-5+error,total_atoms-6+error,1.4480,fix=1,mask=mask)
                    atoms.set_distance(total_atoms-5+error,Hnear2,1.4480,fix=1,mask=mask)
                mask[total_atoms-5+error]=0
                total_atoms=len(atoms)
                CnearH=near(total_atoms-5+error,atoms,'H',limit=False)
                del atoms[CnearH]
                Hnear1=numjudge(Hnear1,CnearH)
                i1=numjudge(i1,CnearH)
                i2=numjudge(i2,CnearH)
                CnearH=near(total_atoms-5+error,atoms,'H',limit=False)
                mask[CnearH]=1
                atoms.set_dihedral(i2,Hnear2,total_atoms-5+error,CnearH,-150,mask=mask)
                for k in range(10):
                    atoms.set_angle(Hnear2,total_atoms-5+error,CnearH,109.5,mask=mask)
                    atoms.set_angle(total_atoms-6+error,total_atoms-5+error,CnearH,109.5,mask=mask)
                atoms.set_distance(CnearH,total_atoms-5+error,1.0100,fix=1)
                error=error+1
            elif needelements[i]=='O':
                for k in range(10):
                    atoms.set_distance(total_atoms-5+error,total_atoms-6+error,1.4100,fix=1,mask=mask)
                    atoms.set_distance(total_atoms-5+error,Hnear2,1.4100,fix=1,mask=mask)
                for j in range(2):
                    total_atoms=len(atoms)
                    CnearH=near(total_atoms-5-error+j,atoms,'H',limit=False)
                    del atoms[CnearH]
                error=error+2
        if places[i]==5:
            mask[Hnear2]=1
            elements=atoms.get_chemical_symbols()
            elements[Hnear2]=needelements[i]
            atoms.set_chemical_symbols(elements)
            if needelements[i]=='N':
                for k in range(10):
                    if state[0]==1:
                        atoms.set_distance(Hnear2,i2,1.4480,fix=1,mask=mask)
                    elif state[0]==2 or state[0]==6:
                        atoms.set_distance(Hnear2,i2,1.3400,fix=1,mask=mask)
                    atoms.set_distance(Hnear2,total_atoms-5+error,1.4480,fix=1,mask=mask)
                mask[Hnear2]=0
                total_atoms=len(atoms)
                CnearH=near(Hnear2,atoms,'H',limit=False)
                del atoms[CnearH]
                Hnear1=numjudge(Hnear1,CnearH)
                i1=numjudge(i1,CnearH)
                i2=numjudge(i2,CnearH)
                CnearH=near(Hnear2,atoms,'H',limit=False)
                mask[CnearH]=1
                for k in range(10):
                    atoms.set_angle(i2,Hnear2,CnearH,109.5,mask=mask)
                    atoms.set_angle(total_atoms-5+error,Hnear2,CnearH,109.5,mask=mask)
                atoms.set_dihedral(i1,i2,Hnear2,CnearH,-150,mask=mask)
                atoms.set_distance(CnearH,Hnear2,1.0100,fix=1)
            elif needelements[i]=='O':
                for k in range(10):
                    if state[0]==1:
                        atoms.set_distance(Hnear2,i2,1.4100,fix=1,mask=mask)
                    elif state[0]==2 or state[0]==6:
                        atoms.set_distance(Hnear2,i2,1.3640,fix=1,mask=mask)
                    atoms.set_distance(Hnear2,total_atoms-5+error,1.4100,fix=1,mask=mask)
                for j in range(2):
                    total_atoms=len(atoms)
                    CnearH=near(Hnear2,atoms,'H',limit=False)
                    del atoms[CnearH]
                    Hnear2=numjudge(Hnear2,CnearH)
    return atoms

def bindconj5(atoms,bond,state=[1,1],strict={1:'X'}):
    try:
        if strict[1]=='X':
            strict[3]='N'
    except:
        pass
    i1=bond[0]
    i2=bond[1]
    Cnear1=near(i1,atoms,'Not H')
    Cnear2=near(i2,atoms,'Not H',reverse=True)
    total_atoms=len(atoms)
    mask=[0]*total_atoms
    mask[i1]=1
    mask[i2]=1
    for j in range(3):
        try:
            element=strict[j+2]
            place=j+2
        except:
            pass
    bondjudge5([Cnear1,Cnear2],atoms,[i1,i2],state,mask,place)
    for j in range(2):
        CnearH=near(i1,atoms,'H')
        del atoms[CnearH]
        i1=numjudge(i1,CnearH)
        i2=numjudge(i2,CnearH)
    for j in range(2):
        CnearH=near(i2,atoms,'H')
        del atoms[CnearH]
        i1=numjudge(i1,CnearH)
        i2=numjudge(i2,CnearH)
    if place==2:
        atoms.append(element)
        atoms.append('C')
        atoms.append('C')
    elif place==3:
        atoms.append('C')
        atoms.append('C')
        atoms.append(element)
    elif place==4:
        atoms.append('C')
        atoms.append(element)
        atoms.append('C')
    for j in range(3):
        atoms.append('H')
    x=atoms.get_positions()
    x[-6][0]=x[-6][0]-100
    x[-6][1]=x[-6][1]-50
    x[-6][2]=x[-6][2]-50
    x[-5][0]=x[-5][0]+50
    x[-5][1]=x[-5][1]+50
    x[-5][2]=x[-5][2]+50
    for j in range(-4,0):
        x[j][0]=x[j][0]+50
        x[j][1]=x[j][1]+50
        x[j][2]=x[j][2]+50
    atoms.set_positions(x)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-6)+[1]*6
    if place==2:
        if element=='N':
            atoms.set_angle(i2,i1,total_atoms-6,107.7,mask=mask)
            if state[0]==1:
                atoms.set_distance(total_atoms-6,i1,1.3810,fix=1,mask=mask)
            elif state[0]==2:
                atoms.set_distance(total_atoms-6,i1,1.3740,fix=1,mask=mask)
        elif element=='O':
            atoms.set_angle(i2,i1,total_atoms-6,110.6,mask=mask)
            atoms.set_distance(total_atoms-6,i1,1.3600,fix=1,mask=mask)
        mask[-6]=0
        atoms.set_angle(i1,i2,total_atoms-5,107.3,mask=mask)
        atoms.set_distance(total_atoms-5,i2,1.4240,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-6,i1,i2,total_atoms-5,0.5,mask=mask)
        mask[-5]=0
        if element=='N':
            for k in range(10):
                atoms.set_distance(total_atoms-4,total_atoms-6,1.3810,fix=1,mask=mask)
                atoms.set_distance(total_atoms-4,total_atoms-5,1.3670,fix=1,mask=mask)
            mask[-4]=0
            #atoms.set_distance(total_atoms-1,total_atoms-6,0,fix=1,mask=mask)
            atoms.set_dihedral(total_atoms-5,total_atoms-4,total_atoms-6,total_atoms-1,0.5,mask=mask)
            for k in range(10):
                atoms.set_angle(total_atoms-4,total_atoms-6,total_atoms-1,124.5,mask=mask)
                atoms.set_angle(i1,total_atoms-6,total_atoms-1,124.5,mask=mask)
            atoms.set_distance(total_atoms-1,total_atoms-6,1.0100,fix=1,mask=mask)
            mask[-1]=0
        elif element=='O':
            for k in range(10):
                atoms.set_distance(total_atoms-4,total_atoms-6,1.3600,fix=1,mask=mask)
                atoms.set_distance(total_atoms-4,total_atoms-5,1.3670,fix=1,mask=mask)
            mask[-4]=0
        atoms.set_distance(total_atoms-3,total_atoms-4,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-3,total_atoms-4,total_atoms-5,i2,0.5,mask=mask)
        for k in range(10):
            atoms.set_angle(total_atoms-6,total_atoms-4,total_atoms-3,126.5,mask=mask)
            atoms.set_angle(total_atoms-5,total_atoms-4,total_atoms-3,126.5,mask=mask)
        mask[-3]=0
        atoms.set_distance(total_atoms-2,total_atoms-5,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-2,total_atoms-5,i2,i1,0.5,mask=mask)
        for k in range(10):
            atoms.set_angle(i2,total_atoms-5,total_atoms-2,126.5,mask=mask)
            atoms.set_angle(total_atoms-4,total_atoms-5,total_atoms-2,126.5,mask=mask)
        mask[-2]=0
    elif place==3:
        atoms.set_angle(i2,i1,total_atoms-6,107.3,mask=mask)
        atoms.set_distance(total_atoms-6,i1,1.3670,fix=1,mask=mask)
        mask[-6]=0
        atoms.set_angle(i1,i2,total_atoms-5,107.3,mask=mask)
        atoms.set_distance(total_atoms-5,i2,1.3670,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-6,i1,i2,total_atoms-5,0.5,mask=mask)
        mask[-5]=0
        if element=='N':
            for k in range(10):
                atoms.set_distance(total_atoms-4,total_atoms-6,1.3810,fix=1,mask=mask)
                atoms.set_distance(total_atoms-4,total_atoms-5,1.3810,fix=1,mask=mask)
            mask[-4]=0
            atoms.set_distance(total_atoms-1,total_atoms-4,0,fix=1,mask=mask)
            atoms.set_distance(total_atoms-1,total_atoms-4,1.0100,fix=1,mask=mask)
            atoms.set_dihedral(i1,total_atoms-6,total_atoms-4,total_atoms-1,110,mask=mask)
            for k in range(10):
                atoms.set_angle(total_atoms-6,total_atoms-4,total_atoms-1,124.5,mask=mask)
                atoms.set_angle(total_atoms-5,total_atoms-4,total_atoms-1,124.5,mask=mask)
            mask[-1]=0
        elif element=='O':
            for k in range(10):
                atoms.set_distance(total_atoms-4,total_atoms-6,1.3600,fix=1,mask=mask)
                atoms.set_distance(total_atoms-4,total_atoms-5,1.3600,fix=1,mask=mask)
            mask[-4]=0
        atoms.set_distance(total_atoms-3,total_atoms-6,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-3,total_atoms-6,i1,i2,0.5,mask=mask)
        for k in range(10):
            atoms.set_angle(i1,total_atoms-6,total_atoms-3,126.5,mask=mask)
            atoms.set_angle(total_atoms-4,total_atoms-6,total_atoms-3,126.5,mask=mask)
        mask[-3]=0
        atoms.set_distance(total_atoms-2,total_atoms-5,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-2,total_atoms-5,i2,i1,0.5,mask=mask)
        for k in range(10):
            atoms.set_angle(i2,total_atoms-5,total_atoms-2,126.5,mask=mask)
            atoms.set_angle(total_atoms-4,total_atoms-5,total_atoms-2,126.5,mask=mask)
        mask[-2]=0
    elif place==4:
        atoms.set_angle(i2,i1,total_atoms-6,107.3,mask=mask)
        atoms.set_distance(total_atoms-6,i1,1.4240,fix=1,mask=mask)
        mask[-6]=0
        if element=='N':
            atoms.set_angle(i1,i2,total_atoms-5,107.7,mask=mask)
            if state[0]==1:
                atoms.set_distance(total_atoms-5,i2,1.3810,fix=1,mask=mask)
            if state[0]==2:
                atoms.set_distance(total_atoms-5,i2,1.3740,fix=1,mask=mask)
        elif element=='O':
            atoms.set_angle(i1,i2,total_atoms-5,110.6,mask=mask)
            atoms.set_distance(total_atoms-5,i2,1.3600,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-6,i1,i2,total_atoms-5,0.5,mask=mask)
        mask[-5]=0
        if element=='N':
            for k in range(10):
                atoms.set_distance(total_atoms-4,total_atoms-5,1.3810,fix=1,mask=mask)
                atoms.set_distance(total_atoms-4,total_atoms-6,1.3670,fix=1,mask=mask)
            mask[-4]=0
            #atoms.set_distance(total_atoms-1,total_atoms-6,0,fix=1,mask=mask)
            atoms.set_dihedral(total_atoms-6,total_atoms-4,total_atoms-5,total_atoms-1,0.5,mask=mask)
            for k in range(10):
                atoms.set_angle(total_atoms-4,total_atoms-5,total_atoms-1,124.5,mask=mask)
                atoms.set_angle(i2,total_atoms-5,total_atoms-1,124.5,mask=mask)
            atoms.set_distance(total_atoms-1,total_atoms-5,1.0100,fix=1,mask=mask)
            mask[-1]=0
        elif element=='O':
            for k in range(10):
                atoms.set_distance(total_atoms-4,total_atoms-5,1.3600,fix=1,mask=mask)
                atoms.set_distance(total_atoms-4,total_atoms-6,1.3670,fix=1,mask=mask)
            mask[-4]=0
        atoms.set_distance(total_atoms-3,total_atoms-4,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-3,total_atoms-4,total_atoms-5,i2,0.5,mask=mask)
        for k in range(10):
            atoms.set_angle(total_atoms-6,total_atoms-4,total_atoms-3,126.5,mask=mask)
            atoms.set_angle(total_atoms-5,total_atoms-4,total_atoms-3,126.5,mask=mask)
        mask[-3]=0
        atoms.set_distance(total_atoms-2,total_atoms-6,1.0100,fix=1,mask=mask)
        atoms.set_dihedral(total_atoms-2,total_atoms-6,total_atoms-4,total_atoms-5,0.5,mask=mask)
        for k in range(10):
            atoms.set_angle(i1,total_atoms-6,total_atoms-2,126.5,mask=mask)
            atoms.set_angle(total_atoms-4,total_atoms-6,total_atoms-2,126.5,mask=mask)
        mask[-2]=0
    if element!='N':
        del atoms[-1]
    return atoms

def bindconj6(atoms,bond,state=[1,1],strict={1:'X'}):
    i1=bond[0]
    i2=bond[1]
    Cnear1=near(i1,atoms,'Not H')
    Cnear2=near(i2,atoms,'Not H',reverse=True)
    total_atoms=len(atoms)
    mask=[0]*total_atoms
    mask[i1]=1
    mask[i2]=1
    for j in range(4):
        try:
            element=strict[j+2]
            place=j+2
            break
        except:
            place=1
    bondjudge6([Cnear1,Cnear2],atoms,[i1,i2],state,mask,place)
    for j in range(2):
        CnearH=near(i1,atoms,'H')
        del atoms[CnearH]
        i1=numjudge(i1,CnearH)
        i2=numjudge(i2,CnearH)
    for j in range(2):
        CnearH=near(i2,atoms,'H')
        del atoms[CnearH]
        i1=numjudge(i1,CnearH)
        i2=numjudge(i2,CnearH)    
    if place==2:
        atoms.append(element)
        atoms.append('C')
        atoms.append('C')
        atoms.append('C')
    elif place==3:
        atoms.append('C')
        atoms.append('C')
        atoms.append(element)
        atoms.append('C')
    elif place==4:
        atoms.append('C')
        atoms.append('C')
        atoms.append('C')
        atoms.append(element)
    elif place==5:
        atoms.append('C')
        atoms.append(element)
        atoms.append('C')
        atoms.append('C')
    else:
        for j in range(4):
            atoms.append('C')
    for j in range(4):
        atoms.append('H')
    x=atoms.get_positions()
    x[-8][0]=x[-8][0]-100
    x[-8][1]=x[-8][1]-50
    x[-8][2]=x[-8][2]-50
    x[-7][0]=x[-7][0]+50
    x[-7][1]=x[-7][1]+50
    x[-7][2]=x[-7][2]+50
    x[-6][0]=x[-6][0]-50
    x[-6][1]=x[-6][1]-50
    x[-6][2]=x[-6][2]-50
    x[-5][0]=x[-5][0]+50
    x[-5][1]=x[-5][1]+50
    x[-5][2]=x[-5][2]+50
    for j in range(-2,0):
        x[2*j][0]=x[2*j][0]+50
        x[2*j][1]=x[2*j][1]+50
        x[2*j][2]=x[2*j][2]+50
        x[2*j+1][0]=x[2*j+1][0]-50
        x[2*j+1][1]=x[2*j+1][1]-50
        x[2*j+1][2]=x[2*j+1][2]-50
    atoms.set_positions(x)
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-8)+[1]*8
    atoms.set_angle(i2,i1,total_atoms-8,120,mask=mask)
    if state[0]==2:
        atoms.set_distance(total_atoms-8,i1,1.4040,fix=1,mask=mask)
    elif state[0]==82:
        atoms.set_distance(total_atoms-8,i1,1.4330,fix=1,mask=mask)
    else:
        atoms.set_distance(total_atoms-8,i1,1.4000,fix=1,mask=mask)
    mask[-8]=0
    atoms.set_angle(i1,i2,total_atoms-7,120,mask=mask)
    if state[0]==2:
        atoms.set_distance(total_atoms-7,i2,1.4040,fix=1,mask=mask)
    elif state[0]==82:
        atoms.set_distance(total_atoms-7,i2,1.4330,fix=1,mask=mask)
    else:
        atoms.set_distance(total_atoms-7,i2,1.4000,fix=1,mask=mask)
    if state[0]==2 or state[0]==82:
        atoms.set_dihedral(total_atoms-8,i1,i2,total_atoms-7,5,mask=mask)
    else:
        atoms.set_dihedral(total_atoms-8,i1,i2,total_atoms-7,0.5,mask=mask)
    mask[-7]=0
    atoms.set_angle(i1,total_atoms-8,total_atoms-6,120,mask=mask)
    if state[0]==82:
        atoms.set_distance(total_atoms-6,total_atoms-8,1.4330,fix=1,mask=mask)
    else:
        atoms.set_distance(total_atoms-6,total_atoms-8,1.4000,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-6,total_atoms-8,i1,i2,-0.5,mask=mask)
    mask[-6]=0
    atoms.set_angle(i2,total_atoms-7,total_atoms-5,120,mask=mask)
    if state[0]==82:
        atoms.set_distance(total_atoms-5,total_atoms-7,1.4330,fix=1,mask=mask)
    else:
        atoms.set_distance(total_atoms-5,total_atoms-7,1.4000,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-5,total_atoms-7,i2,i1,-0.5,mask=mask)
    mask[-5]=0
    atoms.set_angle(total_atoms-4,total_atoms-8,i1,120,mask=mask)
    atoms.set_distance(total_atoms-4,total_atoms-8,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(i2,i1,total_atoms-8,total_atoms-4,-179.5,mask=mask)
    mask[-4]=0
    atoms.set_angle(total_atoms-3,total_atoms-7,i2,120,mask=mask)
    atoms.set_distance(total_atoms-3,total_atoms-7,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(i1,i2,total_atoms-7,total_atoms-3,-179.5,mask=mask)
    mask[-3]=0
    atoms.set_angle(total_atoms-8,total_atoms-6,total_atoms-2,120,mask=mask)
    atoms.set_dihedral(total_atoms-2,total_atoms-6,total_atoms-8,i1,179.5,mask=mask)
    atoms.set_distance(total_atoms-2,total_atoms-6,1.0800,fix=1,mask=mask)
    mask[-2]=0
    atoms.set_angle(total_atoms-7,total_atoms-5,total_atoms-1,120,mask=mask)
    atoms.set_distance(total_atoms-1,total_atoms-5,1.0800,fix=1,mask=mask)
    atoms.set_dihedral(total_atoms-8,total_atoms-6,total_atoms-7,total_atoms-1,179.5,mask=mask)
    if place==2:
        elements=atoms.get_chemical_symbols()
        elements[-8]='N'
        atoms.set_chemical_symbols(elements)
        for k in range(10):
            atoms.set_distance(total_atoms-8,total_atoms-6,1.3390,fix=1)
            if state[0]==2:
                atoms.set_distance(total_atoms-8,i1,1.3540,fix=1)
            else:
                atoms.set_distance(total_atoms-8,i1,1.3390,fix=1)
        CnearH=near(total_atoms-8,atoms,'H')
        del atoms[CnearH]
    elif place==3:
        elements=atoms.get_chemical_symbols()
        elements[-6]='N'
        atoms.set_chemical_symbols(elements)
        for k in range(10):
            atoms.set_distance(total_atoms-6,total_atoms-8,1.3390,fix=1)
            atoms.set_distance(total_atoms-6,total_atoms-5,1.3390,fix=1)
        CnearH=near(total_atoms-6,atoms,'H')
        del atoms[CnearH]
    elif place==4:
        elements=atoms.get_chemical_symbols()
        elements[-5]='N'
        atoms.set_chemical_symbols(elements)
        for k in range(10):
            atoms.set_distance(total_atoms-5,total_atoms-6,1.3390,fix=1)
            atoms.set_distance(total_atoms-5,total_atoms-7,1.3390,fix=1)
        CnearH=near(total_atoms-5,atoms,'H')
        del atoms[CnearH]
    elif place==5:
        elements=atoms.get_chemical_symbols()
        elements[-7]='N'
        atoms.set_chemical_symbols(elements)
        for k in range(10):
            if state[0]==2:
                atoms.set_distance(total_atoms-7,i2,1.3540,fix=1)
            else:
                atoms.set_distance(total_atoms-7,i2,1.3390,fix=1)
            atoms.set_distance(total_atoms-7,total_atoms-5,1.3390,fix=1)
        CnearH=near(total_atoms-7,atoms,'H')
        del atoms[CnearH]
    return atoms
        
def numjudge(i,deleted):
    if deleted<i:
        return i-1
    else:
        return i

def bondjudge2(Cnears,atoms,iis,elements,state,mask,times=10):
    print(Cnears,iis,elements)
    if 'N' in elements:
        for k in range(times):
            if state[0]==state[1] and state[0]==2:
                atoms.set_distance(iis[0],Cnears[0],1.5095,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.3400,fix=1,mask=mask)
                atoms.set_distance(iis[0],iis[1],1.2900,fix=1,mask=mask)
        return 7
    for k in range(times):
        if state[0]==state[1] and state[0]==2:
            atoms.set_distance(iis[0],Cnears[0],1.4000,fix=1,mask=mask)
            atoms.set_distance(iis[1],Cnears[1],1.4000,fix=1,mask=mask)
            atoms.set_distance(iis[0],iis[1],1.4000,fix=1,mask=mask)
        else:
            for i in range(2):
                if state[i]==1:
                    atoms.set_distance(iis[i],Cnears[i],1.5095,fix=1,mask=mask)
                elif state[i]==2:
                    atoms.set_distance(iis[i],Cnears[i],1.4600,fix=1,mask=mask)
                elif state[i]==7:
                    atoms.set_distance(iis[i],Cnears[i],1.3400,fix=1,mask=mask)
                elif state[i]==8:
                    atoms.set_distance(iis[i],Cnears[i],1.3640,fix=1,mask=mask)
                elif state[i]==16:
                    atoms.set_distance(iis[i],Cnears[i],1.7700,fix=1,mask=mask)
                elif state[i]==82:
                    atoms.set_distance(iis[i],Cnears[i],1.4440,fix=1,mask=mask)
            atoms.set_distance(iis[0],iis[1],1.3400,fix=1,mask=mask)

def bondjudge5(Cnears,atoms,iis,state,mask,place,times=10):
    if place==2:
        for k in range(times):
            if state[0]==state[1] and state[0]==1:
                atoms.set_distance(iis[0],Cnears[0],1.5040,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.4950,fix=1,mask=mask)
                atoms.set_distance(i1,i2,1.3670,fix=1,mask=mask)
            elif state[0]==state[1] and state[0]==2:
                atoms.set_distance(iis[0],Cnears[0],1.4040,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.4000,fix=1,mask=mask)
                atoms.set_distance(i1,i2,1.4240,fix=1,mask=mask)
    elif place==3:
        for k in range(times):
            if state[0]==state[1] and state[0]==2:
                atoms.set_distance(iis[0],Cnears[0],1.4000,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.4000,fix=1,mask=mask)
            else:
                for i in range(2):
                    if state[i]==1:
                        atoms.set_distance(iis[i],Cnears[i],1.4950,fix=1,mask=mask)
                    elif state[i]==2:
                        atoms.set_distance(iis[i],Cnears[i],1.4330,fix=1,mask=mask)
                    elif state[i]==7:
                        atoms.set_distance(iis[i],Cnears[i],1.4950,fix=1,mask=mask)
                    elif state[i]==8:
                        atoms.set_distance(iis[i],Cnears[i],1.3640,fix=1,mask=mask)
                    elif state[i]==16:
                        atoms.set_distance(iis[i],Cnears[i],1.7700,fix=1,mask=mask)
                    #elif state[i]==82:
                    #    atoms.set_distance(iis[i],Cnears[i],1.4440?,fix=1,mask=mask)
            atoms.set_distance(i1,i2,1.4240,fix=1,mask=mask)
    elif place==4:
        for k in range(times):
            if state[0]==state[1] and state[0]==1:
                atoms.set_distance(iis[0],Cnears[0],1.4950,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.5040,fix=1,mask=mask)
            elif state[0]==state[1] and state[0]==2:
                atoms.set_distance(iis[0],Cnears[0],1.4000,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.4040,fix=1,mask=mask)
            atoms.set_distance(i1,i2,1.4240,fix=1,mask=mask)

def bondjudge6(Cnears,atoms,iis,state,mask,place,times=10):
    if place==2 or place==5:
        for k in range(times):
            if state[0]==state[1] and state[0]==1:
                atoms.set_distance(iis[0],Cnears[0],1.5100,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.5100,fix=1,mask=mask)
                atoms.set_distance(i1,i2,1.4000,fix=1,mask=mask)
            elif state[0]==state[1] and state[0]==2:
                atoms.set_distance(iis[0],Cnears[0],1.4040,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.4040,fix=1,mask=mask)
                atoms.set_distance(i1,i2,1.3700,fix=1,mask=mask)
    elif place==1 or place==3 or place==4:
        for k in range(times):
            if state[0]==state[1] and state[0]==2:
                atoms.set_distance(iis[0],Cnears[0],1.4040,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.4040,fix=1,mask=mask)
                atoms.set_distance(i1,i2,1.3700,fix=1,mask=mask)
            elif state[0]==state[1] and state[0]==82:
                atoms.set_distance(iis[0],Cnears[0],1.4190,fix=1,mask=mask)
                atoms.set_distance(iis[1],Cnears[1],1.4190,fix=1,mask=mask)
                atoms.set_distance(i1,i2,1.3700,fix=1,mask=mask)
            else:
                for i in range(2):
                    if state[i]==1:
                        atoms.set_distance(iis[i],Cnears[i],1.5105,fix=1,mask=mask)
                    elif state[i]==2:
                        atoms.set_distance(iis[i],Cnears[i],1.4330,fix=1,mask=mask)
                    elif state[i]==7:
                        atoms.set_distance(iis[i],Cnears[i],1.3400,fix=1,mask=mask)
                    elif state[i]==8:
                        atoms.set_distance(iis[i],Cnears[i],1.3640,fix=1,mask=mask)
                    elif state[i]==16:
                        atoms.set_distance(iis[i],Cnears[i],1.7700,fix=1,mask=mask)
                    elif state[i]==82:
                        atoms.set_distance(iis[i],Cnears[i],1.4190,fix=1,mask=mask)
                    if state[i]==2:
                        atoms.set_distance(i1,i2,1.4040,fix=1,mask=mask)
                    else:
                        atoms.set_distance(i1,i2,1.4000,fix=1,mask=mask)

def C_selection(substitution,atoms,Cnear,state=[1,1]):
    if substitution=='CO' or substitution=='C=O':
        return CO(atoms,Cnear,state)
    if substitution=='CNH' or substitution=='C=NH':
        return CNH(atoms,Cnear,state)
    if substitution=='O':
        return O(atoms,Cnear,state)
    if substitution=='S':
        return S(atoms,Cnear,state)
    if substitution=='N':
        return N(atoms,Cnear,state)
                        
def CO(atoms,Cnear,state=[1,1]):
    atoms.append('O')
    total_atoms=len(atoms)
    x=atoms.get_positions()
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    atoms.set_positions(x)
    CnearC1=near(Cnear,atoms,'C')
    CnearC2=near(Cnear,atoms,'C',ignore=[CnearC1],reverse=True)
    #print(CnearC1,CnearC2)
    mask=[0]*(total_atoms-1)+[1]
    mask[Cnear]=1
    for k in range(10):
        if state[0]==1:
            atoms.set_distance(Cnear,CnearC1,1.5220,fix=1,mask=mask)
        elif state[0]==2:
            atoms.set_distance(Cnear,CnearC1,1.4440,fix=1,mask=mask)
        elif state[0]==6:
            atoms.set_distance(Cnear,CnearC1,1.4900,fix=1,mask=mask)
        elif state[0]==7:
            atoms.set_distance(Cnear,CnearC1,1.3350,fix=1,mask=mask)
        elif state[0]==8:
            atoms.set_distance(Cnear,CnearC1,1.3640,fix=1,mask=mask)
        elif state[0]==82:
            atoms.set_distance(Cnear,CnearC1,1.5100,fix=1,mask=mask)
        if state[1]==1:
            atoms.set_distance(Cnear,CnearC2,1.5220,fix=1,mask=mask)
        elif state[1]==2:
            atoms.set_distance(Cnear,CnearC2,1.4440,fix=1,mask=mask)
        elif state[1]==6:
            atoms.set_distance(Cnear,CnearC2,1.4900,fix=1,mask=mask)
        elif state[1]==7:
            atoms.set_distance(Cnear,CnearC2,1.3350,fix=1,mask=mask)
        elif state[1]==8:
            atoms.set_distance(Cnear,CnearC2,1.3640,fix=1,mask=mask)
        elif state[1]==82:
            atoms.set_distance(Cnear,CnearC2,1.5100,fix=1,mask=mask)
    mask[Cnear]=0
    for k in range(10):
        atoms.set_angle(CnearC1,Cnear,total_atoms-1,120,mask=mask)
        atoms.set_angle(CnearC2,Cnear,total_atoms-1,120,mask=mask)
    atoms.set_distance(total_atoms-1,Cnear,1.2290,fix=1,mask=mask)
    for i in range(2):
        Hnear=near(Cnear,atoms,'H',limit=False)
        del atoms[Hnear]
    return atoms

def CNH(atoms,Cnear,state=[1,1]):
    atoms.append('N')
    atoms.append('H')
    total_atoms=len(atoms)
    x=atoms.get_positions()
    x[-2][0]=x[-2][0]+50
    x[-2][1]=x[-2][1]+50
    x[-2][2]=x[-2][2]+50
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    atoms.set_positions(x)
    CnearC1=near(Cnear,atoms,'C')
    CnearC2=near(Cnear,atoms,'C',ignore=[CnearC1],reverse=True)
    #print(CnearC1,CnearC2)
    mask=[0]*(total_atoms-2)+[2]
    mask[Cnear]=1
    for k in range(10):
        if state[0]==1:
            atoms.set_distance(Cnear,CnearC1,1.5220,fix=1,mask=mask)
        elif state[0]==2:
            atoms.set_distance(Cnear,CnearC1,1.5220,fix=1,mask=mask)
        elif state[0]==6:
            atoms.set_distance(Cnear,CnearC1,1.5220,fix=1,mask=mask)
        elif state[0]==7:
            atoms.set_distance(Cnear,CnearC1,1.3350,fix=1,mask=mask)
        elif state[0]==8:
            atoms.set_distance(Cnear,CnearC1,1.3640,fix=1,mask=mask)
        elif state[0]==82:
            atoms.set_distance(Cnear,CnearC1,1.5100,fix=1,mask=mask)
        if state[1]==1:
            atoms.set_distance(Cnear,CnearC2,1.5220,fix=1,mask=mask)
        elif state[1]==2:
            atoms.set_distance(Cnear,CnearC2,1.5220,fix=1,mask=mask)
        elif state[1]==6:
            atoms.set_distance(Cnear,CnearC2,1.5220,fix=1,mask=mask)
        elif state[1]==7:
            atoms.set_distance(Cnear,CnearC2,1.3350,fix=1,mask=mask)
        elif state[1]==8:
            atoms.set_distance(Cnear,CnearC2,1.3640,fix=1,mask=mask)
        elif state[1]==82:
            atoms.set_distance(Cnear,CnearC2,1.5100,fix=1,mask=mask)
    mask[Cnear]=0
    for k in range(10):
        atoms.set_angle(CnearC1,Cnear,total_atoms-2,120,mask=mask)
        atoms.set_angle(CnearC2,Cnear,total_atoms-2,120,mask=mask)
    mask[-2]=0
    atoms.set_angle(Cnear,total_atoms-2,total_atoms-1,120.5,mask=mask)
    atoms.set_distance(total_atoms-1,total_atoms-2,1.0100,fix=1,mask=mask)
    atoms.set_dihedral(Cnear1,Cnear,total_atoms-2,total_atoms-1,179.5,mask=mask)
    for i in range(2):
        Hnear=near(Cnear,atoms,'H',limit=False)
        del atoms[Hnear]
    return atoms

def O(atoms,Cnear,state=[1,1]):
    atoms.append('O')
    total_atoms=len(atoms)
    x=atoms.get_positions()
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    atoms.set_positions(x)
    CnearC1=near(Cnear,atoms,'C')
    CnearC2=near(Cnear,atoms,'C',ignore=[CnearC1],reverse=True)
    print(CnearC1,CnearC2)
    mask=[0]*(total_atoms-1)+[1]
    atoms.set_distance(total_atoms-1,Cnear,0,fix=1)
    for k in range(10):
        if state[0]==1:
            atoms.set_distance(total_atoms-1,CnearC1,1.410,fix=1)
        elif state[0]==2 or state[0]==6 or state[0]==82:
            atoms.set_distance(total_atoms-1,CnearC1,1.364,fix=1)
        elif state[0]==4:
            atoms.set_distance(total_atoms-1,CnearC1,1.327,fix=1)
        elif state[0]==164:
            atoms.set_distance(total_atoms-1,CnearC1,1.670,fix=1)
        if state[1]==1:
            atoms.set_distance(total_atoms-1,CnearC2,1.410,fix=1)
        elif state[1]==2 or state[1]==6 or state[1]==82:
            atoms.set_distance(total_atoms-1,CnearC2,1.364,fix=1)
        elif state[1]==4:
            atoms.set_distance(total_atoms-1,CnearC2,1.327,fix=1)
        elif state[1]==164:
            atoms.set_distance(total_atoms-1,CnearC2,1.670,fix=1)
    del atoms[Cnear]
    for i in range(2):
        total_atoms=len(atoms)
        Hnear=near(total_atoms-1,atoms,'H',limit=False)
        del atoms[Hnear]
    return atoms

def S(atoms,Cnear,state=[1,1]):
    atoms.append('S')
    total_atoms=len(atoms)
    x=atoms.get_positions()
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    atoms.set_positions(x)
    CnearC1=near(Cnear,atoms,'C')
    CnearC2=near(Cnear,atoms,'C',ignore=[CnearC1],reverse=True)
    print(CnearC1,CnearC2)
    mask=[0]*(total_atoms-1)+[1]
    atoms.set_distance(total_atoms-1,Cnear,0,fix=1)
    for k in range(10):
        if state[0]==1:
            atoms.set_distance(total_atoms-1,CnearC1,1.810,fix=1)
        elif state[0]==2:
            atoms.set_distance(total_atoms-1,CnearC1,1.770,fix=1)
        if state[1]==1:
            atoms.set_distance(total_atoms-1,CnearC2,1.810,fix=1)
        elif state[1]==2:
            atoms.set_distance(total_atoms-1,CnearC2,1.770,fix=1)
    Hnear=near(total_atoms-1,atoms,'C',limit=False)
    del atoms[Hnear]
    for i in range(2):
        total_atoms=len(atoms)
        Hnear=near(total_atoms-1,atoms,'H',limit=False)
        del atoms[Hnear]
    return atoms

def N(atoms,Cnear,state=[1,1]):
    atoms.append('N')
    total_atoms=len(atoms)
    x=atoms.get_positions()
    x[-1][0]=x[-1][0]-50
    x[-1][1]=x[-1][1]-50
    x[-1][2]=x[-1][2]-50
    atoms.set_positions(x)
    CnearC1=near(Cnear,atoms,'C')
    CnearC2=near(Cnear,atoms,'C',ignore=[CnearC1],reverse=True)
    print(CnearC1,CnearC2)
    mask=[0]*(total_atoms-1)+[1]
    atoms.set_distance(total_atoms-1,Cnear,0,fix=1)
    for k in range(10):
        if state[0]==1:
            atoms.set_distance(total_atoms-1,CnearC1,1.4480,fix=1)
        elif state[0]==2 or state[0]==6 or state[0]==82:
            atoms.set_distance(total_atoms-1,CnearC1,1.3400,fix=1)
        elif state[0]==4:
            atoms.set_distance(total_atoms-1,CnearC1,1.4900,fix=1)
        if state[1]==1:
            atoms.set_distance(total_atoms-1,CnearC2,1.4480,fix=1)
        elif state[1]==2 or state[1]==6 or state[1]==82:
            atoms.set_distance(total_atoms-1,CnearC2,1.3400,fix=1)
        elif state[1]==4:
            atoms.set_distance(total_atoms-1,CnearC2,1.4900,fix=1)
    CnearC1=numjudge(CnearC1,Cnear)
    CnearC2=numjudge(CnearC2,Cnear)
    del atoms[Cnear]
    for i in range(2):
        total_atoms=len(atoms)
        Hnear=near(total_atoms-1,atoms,'H',limit=False)
        del atoms[Hnear]
        CnearC1=numjudge(CnearC1,Hnear)
        CnearC2=numjudge(CnearC2,Hnear)
    atoms.append('H')
    total_atoms=len(atoms)
    mask=[0]*(total_atoms-1)+[1]
    CnearC3=near(CnearC1,atoms,'C',ignore=[CnearC1,CnearC2])
    for k in range(10):
        atoms.set_angle(CnearC1,total_atoms-2,total_atoms-1,109.5,mask=mask)
        atoms.set_angle(CnearC2,total_atoms-2,total_atoms-1,109.5,mask=mask)
        atoms.set_dihedral(CnearC3,CnearC1,total_atoms-2,total_atoms-1,150,mask=mask)
        atoms.set_distance(total_atoms-1,total_atoms-2,1.0100,fix=1)
    return atoms

class molecule:
    def __init__(self, file, loc='./'):
        '''
            To create a molecule platform, powered by ASE.
            molecule(file,loc='./')
            file: Your polymer name on your file.
            loc: File Location. The default is your current location.
            You can get the further information by .sub, .bind, .C, .geoview, .xyzwrite, .straight and .conformer.
        '''
        atoms=read(loc+file+'.xyz')
        self.file=file
        self.atoms=atoms
        self.loc=loc
    
    def sub(self,substitution,i,state=1,place=0,strict={1:'X'},shuffle=0,angle=-1):
        '''
            To add substitution on moiety.
            sub(substitution,i,state=1,place=0,strict={1:'X'},shuffle=0,angle=-1)
            'substitution' Include:
                Halogen (F, Cl, Br, I)
                OH, SH
                OCH3, SCH3
                NH2,NMe2
                Me, Et, Pr, i-Pr, Bu, t-Bu
                CH=CH2
                CX3 (X = F, Cl, Br, I)
                NO2
                CHO, COCH3
                COOH, COOCH3
                C=NH
                CN
                Ph
                Py (pyridine)
                pyrrole
                cyclo rings (3, 5, 6)
                SOMe
                SO2Me
                SO3H
                bicyclo rings ([2:1:2], [2:2:2])
            i: atom number
            state: the bond connecting with this moiety
            1:   with -C-C
            2:   with -C=C
            3:   with -C3C
            6:   with -Ph
            7:   with -N-R2
            8:   with -O-R
            82:  with -C(=O)-R
            16:  with -S-R
            162: with -SO-R
            164: with -SO2-R
            strict: Used on cyclo rings, to change C on ring into O, N or S. You can change two or three unconsecutive atoms at a time.
            angle: Used on -OH, -SH, -OMe, -SMe, to change the angle between the moiety and the connecting bond.
            This function is under construction, you can directly use key words 'place', 'shuffle' to change or modify your molecule structure.
        '''
        self.atoms=sub_selection(substitution,self.atoms,i,state,place,strict,shuffle,angle)
        return self.atoms
        
    
    def bind(self,ring_size,bond,state=[1,1],place=0,strict={1:'X'}, side=0, shuffle=0, less5=False, stable=False):
        '''
            To change a normal bond into a bond bound with ring.
            bind(ring_size,bond=[i1,i2],state=[state1,state2],place,strict,side,shuffle,less5,stable)
            'ring_size' Include:
                2: change C-C bond into C=C bond
                3, 5, 6: bind a 3(5,6)-membered ring on the C-C bond
                65, 66: bind a 5(6)-membered conjugated ring on the C-C bond
            bond: atom numbers in this bond.
            state: TWO bonds connecting with this bond
            1:   with -C-C
            2:   with -C=C
            6:   with -Ph
            strict: Used on cyclo rings, to change C on ring into O, N or S. You can change two or three unconsecutive atoms at a time.
            This function is under construction, you can directly use key words 'side', 'shuffle' and 'stable' to change or modify your molecule structure.
        '''
        self.atoms=bind_selection(ring_size,self.atoms,bond,state,place,strict, side, shuffle, less5, stable)
        return self.atoms
    
    def C(self,substitution,i,state=[1,1]):
        '''
            To change C into other structures.
            C_selection(substitution,atoms,i,state)
            'substitution' Include:
                O, S
                NH
                CO
                CNH
            i: C atom number
            state: TWO bonds connecting with this C atom
            1:   with -C-C
            2:   with -C=C
            3:   with -C3C
            6:   with -Ph
            7:   with -N-R2
            8:   with -O-R
            82:  with -C(=O)-R
            16:  with -S-R
            162: with -SO-R
            164: with -SO2-R
            
            This function is under construction.
        '''
        self.atoms=C_selection(substitution,self.atoms,i,state)
        return self.atoms
    
    def geoview(self):
        '''
            To see the current structures of molecule, powered by ase.
            geoview()
            Have the results as ase.visualize.view(atoms)
            Example:
                Input:
                    from MCPoly.moldraw import molecule
                    atoms=molecule('Atoms1','../')
                    atoms.sub('D5',7)
                    atoms.sub('Me',4)
                    atoms.geoview()
                Output:
                    <<3D picture powered by ase.visualize.view>>
        '''
        return view(self.atoms)
    
    def atomnum(self):
        '''
            To see the atom numbers of molecule, powered by ase.
            atomnum()
        '''
        return len(self.atoms)
    
    def xyzwrite(self, name=''):
        '''
            To save the current structures of molecule, powered by ase.
            xyzwrite(name)
            Have the results as ase.io.write(name, atoms)
            Example:
                Input:
                    from MCPoly.moldraw import molecule
                    atoms=molecule('Atoms1','../')
                    atoms.sub('D5',7)
                    atoms.sub('Me',4)
                    atoms.xyzwrite()
        '''
        if name=='':
            name=self.file
        return write(name+'.xyz',self.atoms)
    
    def straight(self,start,end):
        '''
            To turn the molecule and make sure the start and end of the molecule is parallel with x axis, powered by ase.
            straight(start,end)
            Example:
                Input:
                    from MCPoly.moldraw import molecule
                    atoms=molecule('Atoms1','../')
                    atoms.sub('D5',7)
                    atoms.sub('Me',4)
                    atoms.straight(1,12)
                    atoms.xyzwrite()
                Output in Et.xyz:
                    1 0.163 -1.262 -0.273
                    ...
                    12 4.691 -1.262 -0.273
                    ...
        '''
        pos=self.atoms.get_positions()
        Apos=pos[start]
        Bpos=pos[end]
        self.atoms.rotate(m.atan((Bpos[-1]-Apos[-1])/(Bpos[1]-Apos[1]))/m.pi*180,'x',Apos)
        pos=self.atoms.get_positions()
        Apos=pos[start]
        Bpos=pos[end]
        self.atoms.rotate(m.atan((Bpos[-1]-Apos[-1])/(Bpos[0]-Apos[0]))/m.pi*180,'y',Apos)
        pos=self.atoms.get_positions()
        Apos=pos[start]
        Bpos=pos[end]
        self.atoms.rotate(m.atan(-(Bpos[1]-Apos[1])/(Bpos[0]-Apos[0]))/m.pi*180+180,'z',Apos)
        pos=self.atoms.get_positions()
        return self.atoms
    
    def conformer(self,lowenergy=99999,highenergy=-99999,must=False):
        '''
            To create a conformer of current molecule. Powered by rdkit and py3Dmol. Not original.
            conformer(lowenergy=-99999,highenergy=99999)
            (low/high)energy: You can change it to get molecule with bigger energy or smaller.
                lowenergy:  with lower energy than this figure.
                highenergy:  with higher energy than this figure.
            TIPS: 1. Don't use 'energy' for the first time you use function 'conformer'.
                  2. The calculation method of energy is not very accurate compared with DFT methods.
            Example:
                Input:
                    from MCPoly.moldraw import molecule
                    atoms=molecule('Atoms1','../Atoms1/')
                    atoms.conformer()
                    atoms.conformer(lowenergy=atoms.energy)
                Output:
                    Energy: 291.0767968459029
                    Molecule: generated 100 conformers and 96 clusters
                    <py3Dmol strucutre>
                    Energy: 264.7508728908226
                    Molecule: generated 100 conformers and 97 clusters
                    <py3Dmol strucutre>
        '''
        opath=os.getcwd()
        os.chdir(self.loc)
        createmol(self.file)
        suppl = Chem.MolFromMolFile(self.file+'.mol')
        suppl = Chem.AddHs(suppl)
        conformerIds = gen_conformers(suppl, numConfs=100, maxAttempts=2000)
        conformerPropsDict = []
        for conformerId in conformerIds:
            props = calc_energy(suppl, conformerId, 0)[0]
            conformerPropsDict.append(props)
        try:
            j=0
            while lowenergy<conformerPropsDict[-1]:
                conformerIds = gen_conformers(suppl, numConfs=100, maxAttempts=2000)
                conformerPropsDict = []
                for conformerId in conformerIds:
                    props = calc_energy(suppl, conformerId, 0)[0]
                    conformerPropsDict.append(props)
                if j==100:
                    return "Time has run out. The programme can't find the structure with lower energy."
                j=j+1
            while highenergy>conformerPropsDict[-1]:
                conformerIds = gen_conformers(suppl, numConfs=100, maxAttempts=2000)
                conformerPropsDict = []
                for conformerId in conformerIds:
                    props = calc_energy(suppl, conformerId, 0)[0]
                    conformerPropsDict.append(props)
                if j==100:
                    return "Time has run out. The programme can't find the structure with higher energy."
                j=j+1
        except:
            pass
        print("Energy: {0}".format(conformerPropsDict[-1]))
        self.energy=conformerPropsDict[-1]
        rmsClusters = cluster_conformers(suppl)
        print("Molecule: generated", len(conformerIds), "conformers and", len(rmsClusters), "clusters")
        #print(conformerPropsDict)
        v=draw_with_spheres(suppl)
        if must==True:
            print(Chem.MolToMolBlock(suppl),file=open('{0}_Conformer.mol'.format(self.file),'w+'))
        else:
            xoy=input('Do you want to save this structure? [y/n]')
            if xoy=='y':
                print(Chem.MolToMolBlock(suppl),file=open('{0}_Conformer.mol'.format(self.file),'w+'))
        os.system('rm {0}'.format(self.file+'.mol'))
        os.chdir(opath)