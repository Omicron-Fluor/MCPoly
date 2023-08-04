import os
import warnings
from ase.io import read

def XyzToInp(file, inpname='', fileloc='./', saveloc='./', method='B3LYP',\
             basis_set='def2-SVP', opt=False, freq=False, scan=False,\
             external_force=False, ts=False, aim=[0, 0], stretch=-1, scanstep=10,\
             strain=-1, maxiter=-1, maxcore=-1, corenum=1, electron=0, state=1,\
             MD=False, md_timestep=1.0, initvel=300, tsNHC=300, timecon=10.0,\
             dumpcon=10, runstep=100, dumpname=''):
    """
    The method to convert .xyz file into .inp file.
    XyzToInp(file, inpname='', loc='./', method='B3LYP', basis_set='def2-SVP', opt=False, freq=False, scan=False, external_force=False, ts=False, aim=[0, 0], stretch=1, scanstep=10, strain=-1, maxiter=-1, maxcore=-1, corenum=1, electron=0, state=1)
    file: Your .xyz file name.
    inpname: The name of your saved .inp file.
    fileloc: File Location. The default is your current location.
    saveloc: Input File Save Location. The default is your current location.
    method: Your semiempirical/ab initio/DFT calculation methods. The default is B3LYP.
    basis_set: Your basis set. The default is def2-SVP.
    opt: To show if you want to optimise your molecule. The default is False.
    freq: To show if you want to calculate the frequency of intermolecule system and thermodynamic data. The default is False.
    scan: To show if you want to scan the system. You must use key word 'aim', 'stretch' and 'scanstep' to express the site of external force, or it will cause 'SetError'. 
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
    electron: The charge of your system. It might be positive(>0), negative(<0) or neutral(=0). The default is neutral.
    state: The state of your system (singlet:1, doublet:2, triplet:3 ...). The default is singlet.
    Example 1:
        Input:
            from MCPoly.orcaset import XYZtoINP
            XyzToInp('Atom1')
        Output in Atom1.inp:
            #Powered by MCPoly

            ! B3LYP def2-SVP 

            *xyz 0 1
            C         -7.11449        1.33604       -0.03969
            C         -5.87285        0.46952       -0.06689
            C         -4.65267        1.33725        0.06683
            C         -3.41103        0.47072        0.03973
            F         -7.23198        2.03119        1.11508
                ...
            *
            
            
    Example 2:
        Input:
            from MCPoly.orcaset import XYZtoINP
            XyzToInp('Atom1', method='B3LYP D3BJ', basis_set='def2-TZVP', freq=True, \
                     external_force=True, aim=[0, 3], strain=1.5, maxcore=4096, corenum=8)
        Output in Atom1.inp:
            #Powered by MCPoly
            
            ! B3LYP D3BJ def2-TZVP OPT FREQ
            %maxcore 4096
            %PAL NPROCS 8 END
            %geom
            POTENTIALS
            {C 0 3 1.500}
            end
            end
            
            *xyz 0 1
            C         -7.11449        1.33604       -0.03969
            C         -5.87285        0.46952       -0.06689
            C         -4.65267        1.33725        0.06683
            C         -3.41103        0.47072        0.03973
            F         -7.23198        2.03119        1.11508
                ...
            *
    
    
    Example 3:
        Input:
            from MCPoly.orcaset import XYZtoINP
            XyzToInp('Atom1', method='B3LYP D3BJ', basis_set='def2-TZVP', scan=True, \
                     aim=[0, 3], stretch=0.8, scanstep=8, maxcore=4096, corenum=8)
        Output in Atom1.inp:
            #Powered by MCPoly
            
            ! B3LYP D3BJ def2-TZVP 
            %maxcore 4096
            %PAL NPROCS 8 END
            %geom Scan
            B 0 3 = 3.804, 4.604, 8
            end
            end
            
            *xyz 0 1
            C         -7.11449        1.33604       -0.03969
            C         -5.87285        0.46952       -0.06689
            C         -4.65267        1.33725        0.06683
            C         -3.41103        0.47072        0.03973
            F         -7.23198        2.03119        1.11508
                ...
            *
    
    Example 4:
        Input:
            from MCPoly.orcaset import XYZtoINP
            XyzToInp('Atom1', method='B3LYP D3BJ', basis_set='def2-TZVP', MD=True, \
                     MD=True, freq=True, md_timestep=1.0, initvel=350, tsNHC=350, dumpcon=25, runstep=2500, maxcore=4096, corenum=8)
        Output in Atom1.inp:
            #Powered by MCPoly
            
            ! MD B3LYP D3BJ def2-TZVP 
            %md
            Timestep 1.0_fs
            Initvel 350_K
            Thermostat NHC 350_K Timecon 10.0_fs
            Dump Position Stride 25 Filename "S4A_dump.xyz"
            Run 2500
            end
            %maxcore 4096
            %PAL NPROCS 8 END
            
            *xyz 0 1
            C         -7.11449        1.33604       -0.03969
            C         -5.87285        0.46952       -0.06689
            C         -4.65267        1.33725        0.06683
            C         -3.41103        0.47072        0.03973
            F         -7.23198        2.03119        1.11508
                ...
            *
    """
    f = open(fileloc+file+'.xyz', 'r')
    if inpname == '':
        try:
            w = open(saveloc+file+'.inp', 'x')
        except:
            w = open(saveloc+file+'.inp', 'w')
    else:
        try:
            w = open(saveloc+inpname+'.inp', 'x')
        except:
            w = open(saveloc+inpname+'.inp', 'w')
    w.write('#Powered by MCPoly\n\n')
    w.write('! ')
    if MD == True:
        w.write('MD ')
        sacn = False
        opt = False
        freq = False
        external_force = False
    
    w.write('{0} {1} '.format(method, basis_set))
        
    if MD == True:
        w.write('\n')
        if initvel < 0:
            raise ValueError("Key word 'initvel' must be positive,\
                for the unit is Kelvin.")
        if tsNHC < 0:
            raise ValueError("Key word 'tsNHC' must be positive,\
                for the unit is Kelvin.")
        w.write('%md\n')
        w.write('Timestep {0}_fs\n'.format(md_timestep))
        w.write('Initvel {0}_K\n'.format(initvel))
        w.write('Thermostat NHC {0}_K Timecon {1}_fs\n'.format(tsNHC, timecon))
        if dumpname == '':
            w.write('Dump Position Stride {0} '.format(dumpcon)
                    + 'Filename "{0}_dump.xyz"\n'.format(file))
        else:
            w.write('Dump Position Stride {0}'.format(dumpcon)
                    + 'Filename "{0}_dump.xyz"\n'.format(dumpname))
        w.write('Run {0}\n'.format(runstep))
        w.write('end\n')
        
    if opt == True:
        if scan == True:
            pass
        if ts == True:
            w.write('OPTTS ')
        else:
            w.write('OPT ')
    if external_force == True and opt == False:
        w.write('OPT ')
    if freq == True:
        if opt == False:
            warnings.warn('Your set the frequency calculation'
                          + ' without geometry optimisation.')
        w.write('FREQ')
    w.write('\n')
    
    if maxiter != -1:
        if maxiter <= 0 or type(maxiter) == float:
            raise ValueError("Key word 'maxiter' must be a positive integar.")
        w.write('%geom\n    MAXITER {0}\nend\n'.format(maxiter))
    
    if maxcore != -1:
        w.write('%maxcore {0}\n'.format(maxcore))
        if maxcore <= 2048:
            warnings.warn('Your max core space is only {0:.3f} Å.'.format(maxcore)
                + 'You might delete this key words to get the same calculation speed,'
                + 'or try bigger one.')

    if corenum != 1:
        w.write('%PAL NPROCS {0} END\n'.format(corenum))
    
    if scan == True:
        if aim == [0, 0]:
            raise ValueError("You haven't set the site of two atoms you scan.")
        elif (aim[0] < 0 or aim[1] < 0 or 
              type(aim[0]) == float or type(aim[1]) == float):
            raise ValueError("The site of two atoms must be not negative integar.")
        elif aim[0] == aim[1]:
            raise ValueError("You can't scan on the same two atoms.")
        if scanstep < 0 or type(scanstep) == float:
            raise ValueError("The scan steps must be positive integar.")
        if ts == True:
            raise AssertionError("You can't set scan and"
                                 + "transition state at the same time.")
        if freq == True:
            raise AssertionError("You can't set scan and"
                                 + "frequency calculation at the same time.")
        atoms = read(fileloc+file+'.xyz')
        distance = atoms.get_distance(aim[0], aim[1])
        if distance + stretch < 0:
            raise ValueError("You can't set the final distance"
                             + "of atoms less than 0.")
        elif distance + stretch <= 0.8:
            warnings.warn('Your final distance of two atoms'
                          + 'is only {0:.3f} Å.'.format(distance+stretch))
        w.write('%geom Scan\n')
        w.write('B {0} {1} = {2:.3f},  {3:.3f},  {4}\nend\nend\n'.format(*aim,\
                                                                       distance, distance+stretch, scanstep))

    if external_force == True:
        if aim == [0, 0]:
            raise ValueError("You haven't set the site of two atoms you scan.")
        elif (aim[0] < 0 or aim[1] < 0 or
              type(aim[0]) == float or type(aim[1]) == float):
            raise ValueError("The site of two atoms must be not negative integar.")
        elif aim[0] == aim[1]:
            raise ValueError("You can't scan on the same two atoms.")
        if scan == True:
            raise AssertionError("You can't set scan and external force at the same time.")
        if ts == True:
            raise AssertionError("You can't set external force"
                                 + "and transition state at the same time.")
        w.write('%geom\nPOTENTIALS\n')
        w.write('{C ' + '{0} {1} {2:.3f}'.format(*aim, strain) + '}\nend\nend\n')

    w.write('\n*xyz {0} {1}\n'.format(electron, state))
    i = 0
    for line in f:
        i = i + 1
        if i >= 3:
            w.write(line)
    w.write('*\n')
    f.close()
    w.close()