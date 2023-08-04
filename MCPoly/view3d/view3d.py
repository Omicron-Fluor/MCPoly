import os
import py3Dmol
from itertools import product
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolToFile
from rdkit.Chem import rdDetermineBonds

def cleanup_qm9_xyz(fname):
    ind = open(fname).readlines()
    nAts = int(ind[0])
    # There are two smiles in the data: the one from GDB and the one assigned from the
    # 3D coordinates in the QM9 paper using OpenBabel (I think).
    gdb_smi, relax_smi = ind[-2].split()[:2]
    ind[1] = '\n'
    ind = ind[:nAts+2]
    for i in range(2, nAts+2):
        l = ind[i]
        l = l.split('\t')
        l.pop(-1)
        ind[i] = '\t'.join(l) + '\n'
    ind = ''.join(ind)
    return ind, gdb_smi, relax_smi

def draw_with_spheres(mol, width, height):
    v = py3Dmol.view(width=width, height=height)
    IPythonConsole.addMolToView(mol, v)
    v.zoomTo()
    v.setStyle({'sphere':{'radius':0.3}, 'stick':{'radius':0.1}});
    v.show()

def view3d(filename, loc = './', width=300, height=300):
    """
    A method to see the geometry structure of the ORCA optimisation, powered by py3Dmol and rdkit. Not original.
    TIPS: Make sure your .xyz file is in the document, or there will be NoFileFoundError!!!
    
    view3d(file, loc = './', width=300, height=300)
    file: File Name.
    loc: File Location. The default is your current location.
    num: The step of your convergence.
    width, height: The size of your 3D geometry molecule strcuture. Default: 300x300.
    """
    path = os.getcwd()
    os.chdir(loc)
    ind = open('{0}.xyz'.format(filename), 'r+')
    mains = ind.readlines()
    main = ''
    for i, line in enumerate(mains):
        if i == 1:
            main = main + '\n'
        else:
            main = main + mains[i]
    #print(main)
    raw_mol = Chem.MolFromXYZBlock(main)
    conn_mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineConnectivity(conn_mol)
    v = draw_with_spheres(conn_mol, width, height)
    os.chdir(path)
    return v