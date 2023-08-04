import time as t
from .XyzToInp import XyzToInp
from .orca import orca

class multiorca:
    """
    The method to calculate several .inp files in one document consecutively by ORCA.
    multiorca([file1, file2, file3, ...], loc)
    file: Your .xyz file name.
    loc: File Location. The default is your current location.
    You can get the further information by .run.
    """
    def __init__(self, files=[], loc='./'):
        self.loc=loc
        self.files=files
        
    def run(self, orcaloc='./', method='B3LYP', basis_set='def2-SVP', opt=False,\
            freq=False, scan=False, external_force=False, ts=False,\
            aim=[0, 0], stretch=-1, scanstep=10, strain=-1, maxiter=-1,\
            maxcore=-1, corenum=1, electron=0, state=1):
        """
    The method to calculate several .inp files by ORCA and save file at the same location with .inp files.
    run(self, orcaloc='./', method='B3LYP', basis_set='def2-SVP', opt=False, freq=False, scan=False, external_force=False, ts=False, aim=[0, 0], stretch=-1, scanstep=10, strain=-1, maxiter=-1, maxcore=-1, corenum=1, electron=0, state=1)
    orcaloc: Your location of ORCA.
             If you have set ORCA from https://www.orcasoftware.de/tutorials_orca/first_steps/trouble_install.html#path-variable, 
             you can skip this key word.
    method: Your semiempirical/ab initio/DFT calculation methods. The default is B3LYP.
    basis_set: Your basis sets. The default is def2-SVP.
    opt: To show if you want to optimise your molecule. The default is False.
    freq: To show if you want to calculate the frequency of intermolecule system and thermodynamic data. The default is False.
    scan: To show if you want to scan the system. You must use key word 'aim', 'stretch' and 'scanstep' to express the site of external force, or it will cause 'ValueError'. 
        Herein, 'aim' is the site of external force, 'stretch' is the length you need to stretch (>0) or shrink(<0). (Unit: Å), and 'scanstep' is the number of your scan steps.
        The default is False.
    external_force: To show if you want to impose external force to the system. You must use key word 'aim' to express the site of external force AND key word 'strain' to express the force energy (Unit: nN), or it will cause 'ValueError'. The default is False.
    ts: To show if you want to calculate the transition state of the system. The default is False.
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
    Example:
        Input:
            from MCPoly.orcaset import multiorca
            files=multiorca(['Atom1', 'Atom2'])
            orcaloc='./MCPoly/orca/'
            files.run(orcaloc=orcaloc, method='B3LYP', basis_set='def2-TZVP', opt=True, freq=True,\
                        maxcore=4096, corenum=8)
        Output:
            1. Atom1.xyz [Tue Apr 25 15:26:50 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 15:47:35 2023]
            
            2. Atom2.xyz [Tue Apr 25 15:47:36 2023]
            <<COMMAND LINES>>
            Program Completed.[Tue Apr 25 16:21:55 2023]
        """
        i = 1
        for file in self.files:
            print('{0}. '.format(i)+file+'.xyz '+'[{0}]'.format(t.ctime(t.time())))
            XyzToInp(file, fileloc=self.loc, saveloc=self.loc, method=method,\
                     basis_set=basis_set, opt=opt, freq=freq, scan=scan,\
                     external_force=external_force, ts=ts, aim=aim, stretch=stretch,\
                     scanstep=scanstep, strain=strain, maxiter=maxiter, maxcore=maxcore,\
                     corenum=corenum, electron=electron, state=state)
            orca(file, orcaloc, self.loc, self.loc)
            i = i + 1