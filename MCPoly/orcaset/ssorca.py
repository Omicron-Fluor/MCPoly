import os
import re
import sys
from .XyzToInp import XyzToInp
from .orca import orca
from ase.io import read
import time as t

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from status import status

def checkbroken(distance, force):
    long1 = abs(distance[-2] - distance[-1])
    long2 = abs(distance[-3] - distance[-2])
    mul1 = abs(force[-2] - force[-1])
    mul2 = abs(force[-3] - force[-2])
    if long1 / long2 >= 10 * mul1 / mul2:
        return 1
    else:
        return 0

class ssorca:
    """
    The method to calculate .inp file with growing external force in one document consecutively by ORCA.
    ssorca(file, loc)
    file: Your .xyz file name.
    loc: File Location. The default is your current location.
    You can get the further information by .run.
    """
    def __init__(self, file, loc='./'):
        self.loc=loc
        self.file=file
        
    def run(self, orcaloc='./', method='B3LYP', basis_set='def2-SVP', freq=False,\
            aim=[0, 0], strain=0, forcestep={0:1}, maxiter=-1,\
            maxcore=-1, corenum=1, electron=0, state=1):
        """
    The method to calculate .inp file with growing external force by ORCA and save file at the same location with .inp files, powered by ase.
    run(self, orcaloc='./', method='B3LYP', basis_set='def2-SVP', freq=False, external_force=False, aim=[0, 0], strain=0, maxiter=-1, maxcore=-1, corenum=1, electron=0, state=1)
    orcaloc: Your location of ORCA.
             If you have set ORCA from https://www.orcasoftware.de/tutorials_orca/first_steps/trouble_install.html#path-variable, 
             you can skip this key word.
    method: Your semiempirical/ab initio/DFT calculation methods. The default is B3LYP.
    basis_set: Your basis sets. The default is def2-SVP.
    opt: To show if you want to optimise your molecule. The default is False.
    freq: To show if you want to calculate the frequency of intermolecule system and thermodynamic data. The default is False.
    aim: To express two sites of external force.
    strain: To express the beginning force energy (Unit: nN). The default is no external force.
    forcestep: The external force grow per step. You can input multiple external force grow. 
               e.g. {0:1, 5:0.5} -> from 0 nN to 5 nN, it grows 1 nN per step; when > 5 nN, it grows 0.5 nN per step.
               The default is {0:1}, growing 1 nN per step always.
    TIPS: 1. You needn't rectify if forcestep={0:0.7, 5:0.5}, thought 5.0 nN will not be calculated.
          2. The accuracy of strain and forcestep is 0.001.
    maxiter: The max iteration steps of the optimization. It must perform after key word 'optimise' is True, or it will cause 'ValueError'. The default is based on ORCA.
    maxcore: The max space you need to calculate. The default is based on ORCA.
    corenum: The core you need to calculate. The default is 1.
    TIPS: To use key word 'maxcore' and 'corenum', please make sure you have set ORCA from following two websites.
            https://www.orcasoftware.de/tutorials_orca/first_steps/trouble_install.html#path-variable
            https://www.orcasoftware.de/tutorials_orca/first_steps/parallel.html
    electron: The charge of your systems. It might be positive(>0), negative(<0) or neutral(=0). The default is neutral.
    state: The state of your systems (singlet:1, doublet:2, triplet:3 ...). The default is singlet.
    TIPS: 1. Don't put the .inp to your ORCA installment file, it may disturb your work.     
          2. If you don't input key word 'orcaloc' while you don't set ORCA, you will only get an empty .out file in the document.
             There will be a warning: sh: orca: 'command not found'.
          3. If you have set key word 'maxcore' and 'corenum', you must input 'orcaloc'.
             See in https://www.orcasoftware.de/tutorials_orca/first_steps/parallel.html
          4. You can see the error of ORCA directly from the output.
    Example 1:
        Input:
            from MCPoly.orcaset import ssorca
            polymer = ssorca('Atom1')
            orcaloc = './MCPoly/orca/'
            polymer.run(orcaloc=orcaloc, method='B3LYP', basis_set='def2-TZVP', freq=True, aim=[2, 3], forcestep={0:0.5},\
                        maxcore=4096, corenum=8)
        Output:
            File: Atom1.xyz
            Current External Force: 0.000 nN [Tue Apr 25 15:23:49 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 15:43:58 2023]
            
            Current External Force: 0.500 nN [Tue Apr 25 15:43:59 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 15:58:12 2023]
            
            Current External Force: 1.000 nN [Tue Apr 25 15:58:12 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 16:13:34 2023]
            
            Current External Force: 1.500 nN [Tue Apr 25 16:13:34 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 16:39:47 2023]
            
            The polymer is broken. Programme Completed.[Tue Apr 25 16:39:47 2023]
            [file orca_main/run.cpp, line 32327]: ORCA finished with error return - aborting the run     -> COMMAND LINES

    Example 2:
        Input:
            from MCPoly.orcaset import ssorca
            polymer = ssorca('Atom2')
            orcaloc = './MCPoly/orca/'
            polymer.run(orcaloc=orcaloc, method='B3LYP', basis_set='def2-TZVP', freq=True, aim=[2, 3], forcestep={0:1, 5:0.5},\
                        maxcore=4096, corenum=8)
        Output:
            File: Atom2.xyz
            Current External Force: 0.000 nN [Tue Apr 25 15:23:49 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 15:43:52 2023]
            
            Current External Force: 1.000 nN [Tue Apr 25 15:43:52 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 15:58:12 2023]
            
            ...
            
            Current External Force: 5.000 nN [Tue Apr 25 17:03:26 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 17:12:03 2023]
            
            Current External Force: 5.500 nN [Tue Apr 25 17:12:03 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 17:24:39 2023]
            
            Current External Force: 6.000 nN [Tue Apr 25 17:24:40 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 17:51:31 2023]
            
            The polymer is broken. Programme Completed.[Tue Apr 25 17:51:31 2023]
            [file orca_main/run.cpp, line 73453]: ORCA finished with error return - aborting the run     -> COMMAND LINES
        """
        distances = []
        force = []
        print('File: {0}.xyz'.format(self.file))
        XyzToInp(self.file, fileloc=self.loc, method=method,\
                 basis_set=basis_set, opt=True, freq=freq, \
                 external_force=True, aim=aim, strain=strain,\
                 maxiter=maxiter, maxcore=maxcore, corenum=corenum, \
                 electron=electron, state=state)
        os.system('cp {0}.inp {0}_0.000.inp'.format(self.file))
        keys = []
        for key in forcestep:
            keys.append(key)
        keys.reverse()
        key = keys[0]
        while 1:
            print('Current External Force: {0:.3f} nN [{1}]'.format(strain,\
                                                                    t.ctime(t.time())))
            orca(self.file+'_{0:.3f}'.format(strain), orcaloc, self.loc, self.loc)
            s = status(self.file+'_{0:.3f}'.format(strain),\
                       self.loc).status(statusonly=True)
            atoms = read(file+'_{0:.3f}.xyz'.format(strain))
            d = atoms.get_distance(aim[0], aim[1])
            if s == 4:
                distances.append(d)
                force.append(strain)
                xoy = checkbroken(distances, force)
                if xoy == 1:
                    print('The polymer is broken.\
                        Programme Completed. [{0}]'.format(t.ctime(t.time())))
                    break
                print('This calculation has some errors.\nThe programme\
                    will recalculate this system again. [{0}]\n'.format(t.ctime(t.time())))
                orca(self.file+'_{0:.3f}'.format(strain), orcaloc, self.loc, self.loc)
                s = status(self.file, self.loc).status(statusonly=True)
                if s == 4:
                    print('This calculation still has some errors.\n\
                        Please check\
                        your .inp and try again. [{0}]\n'.format(t.ctime(t.time())))
                    distances[-1] = None
                else:
                    distances[-1] = d
            else:
                distances.append(d)
                force.append(strain)
            for key in keys:
                if strain >= key:
                    XyzToInp(self.file+'_{0:.3f}'.format(strain),\
                             inpname=self.file+'_{0:.3f}'.format(strain+forcestep[key]), \
                             method=method, basis_set=basis_set, opt=True, freq=freq, \
                             external_force=True, aim=aim, strain = strain + forcestep[key],\
                             maxiter=maxiter, maxcore=maxcore,\
                             corenum=corenum, electron=electron, state=state)
                    strain = strain + forcestep[key]
                    break