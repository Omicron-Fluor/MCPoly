import os
import re
import time as t

def orca(file, orcaloc='./', fileloc='./', saveloc='./'):
    """
    The method to calculate .inp file by ORCA.
    orca(file, orcaloc='./', fileloc='./', saveloc='./')
    file: Your .inp file name.
    orcaloc: Your location of ORCA.
             If you have set ORCA from https://www.orcasoftware.de/tutorials_orca/first_steps/trouble_install.html#path-variable, 
             you can skip this key word.
    fileloc: Your location of .inp file. The default is your current location.
    saveloc: Your location of .inp file. The default is your current location.
    Example:
        Input:
            from MCPoly.orcaset import XYZtoINP
            from MCPoly.orcaset import orca
            
            XYZtoINP('Atom2', method='B3LYP', basis_set='def2-TZVP', opt=True, freq=True, \
                        maxcore=4096, corenum=8)
            orca('Atom2', orcaloc='./MCPoly/orca/')
        Output:
            <<COMMAND LINES>>
            Program Completed.[Mon Apr 17 14:25:12 2023] -> See the current location to check files.
            
    TIPS: 1. Don't put the .inp to your ORCA installment file,  it may disturb your work. 
          2. If you don't input key word 'orcaloc' while you don't set ORCA,  you will only get an empty .out file in the document.
             There will be a warning: sh: orca: 'command not found'.
          3. If you have set key word 'maxcore' and 'corenum',  you must input 'orcaloc'.
             See in https://www.orcasoftware.de/tutorials_orca/first_steps/parallel.html
          4. You can see the error of ORCA directly from the output.
    """

    if orcaloc == './':
        os.system('orca {0}{1}.inp > {2}{1}.out'.format(fileloc, file, saveloc))
    else:
        path = os.getcwd() + '/'
        os.chdir(orcaloc)
        print(os.getcwd())
        if fileloc == './':
            fileloc = path
        if saveloc == './':
            saveloc = path
        os.system('./orca {0}{1}.inp > {2}{1}.out'.format(fileloc, file, saveloc))
        print('Program Completed.[{0}]\n'.format(t.ctime(t.time())))
        os.chdir(path[:-1])
    
    path = os.getcwd() + '/'
    os.chdir(fileloc)
    f = open(file+'.inp', 'r')
    for line in f:
        md = re.search('%md', line)
        if md:
            os.system('cp {0}{1}_dump.xyz {1}_dump.xyz'.format(orcaloc, file))
            os.system('rm {0}{1}_dump.xyz'.format(orcaloc, file))
    os.chdir(path)