import re
import os
import matplotlib.pyplot as plt
from ipywidgets import interact
import ipywidgets as iw
import py3Dmol
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDetermineBonds

def status_judge(s1,s2,converge):
    if s1==0 or s1==4:
        if s1==0:
            print('The Optimization is processing.')
        elif s1==4 and s2==0:
            print('The Optimization is aborted.')
        print('Last time:')
        print('Energy change: {0}'.format(converge[-5]))
        print('RMS gradient: {0}             MAX gradient: {0}'.format(converge[-4],converge[-3]))
        print('RMS step: {0}                 MAX step: {0}'.format(converge[-2],converge[-1]))
    if s1!=0 and s2<0:
        print('The Optimization was finished.')
        if s1==4:
            print('The frequency calculation is aborted.')
            return None
        if s2==-1:
            print('The CP-SCF equations are forming.')
        elif s2==-2:
            print('The CP-SCF equations are solving.')
        elif s2==-3:
            print('The Thermodynamics Calculation was processing.')
    if s2==1:
        print('The Thermodynamics Calculation was finished.')
        
def statusfig(energy,choose,num=-1):
    t=range(len(energy))
    if choose==0:
        plt.plot(t,energy,'x-')
    else:
        plt.plot(t[-choose:],energy[-choose:],'x-')
    if num>=0:
        #plt.plot(num,energy[num],'r',markersize=12)
        plt.scatter(num,energy[num],s=40,c='r')
    plt.show()

def higeo(file,num):
    energy=status(file,figureonly=True)
    statusfig(energy,0,num)
    
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

def draw_with_spheres(mol,width,height):
    v = py3Dmol.view(width=width,height=height)
    IPythonConsole.addMolToView(mol,v)
    v.zoomTo()
    v.setStyle({'sphere':{'radius':0.4},'stick':{'radius':0.1}});
    v.show()
    return v

def multicreatemol(filename,turn,width,height):
    try:
        ind = open('{0}_trj.xyz'.format(filename),'r+')
    except:
        raise FileNotFoundError('Perhaps your _trj.xyz file has some errors. Please check the location of your _trj.xyz file and then modify your file.')
    mains = ind.readlines()
    num=eval(mains[0][:-1])
    term=len(mains)/(num+2)
    main=''
    partmain=mains[turn*(num+2):(turn+1)*(num+2)]
    for i,line in enumerate(partmain):
        if i==1:
            main=main+'\n'
        else:
            main=main+partmain[i]
    #print(main)
    raw_mol = Chem.MolFromXYZBlock(main)
    conn_mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineConnectivity(conn_mol)
    v=draw_with_spheres(conn_mol,width,height)
    return v

def normalstatus(loc,file,choose=0,figureonly=False,statusonly=False):
    path = os.getcwd()
    os.chdir(loc)
    file=file+'.out'
    if figureonly==False and statusonly==False:
        print(file)
    s1=0
    s2=0
    term=0
    f=open(file,'r')
    x=0
    energy=[]
    converge=[]
    for line in f:
        a=re.search('FINAL SINGLE POINT ENERGY',line)
        if a:
            a2=re.search(r'-[0-9]+\.[0-9]+',line)
            energy.append(eval(a2.group(0)))
        x1=re.search('OPTIMIZATION RUN DONE',line)
        if x1:
            s1=1
        if s1==1:
            in1=re.search('Forming right-hand sides of CP-SCF equations',line)
            if in1:
                s2=-1
            in2=re.search('Solving the CP-SCF equations',line)
            if in2:
                s2=-2
            in3=re.search('VIBRATIONAL FREQUENCIES',line)
            if in3:
                s2=-3
        x2=re.search('ORCA TERMINATED NORMALLY',line)
        if x2:
            s2=1
        t=re.search('GEOMETRY OPTIMIZATION CYCLE',line)
        if t:
            t2=re.search(r'[0-9]+',line)
            term=t2.group(0)
        c=re.search('-|Geometry convergence|-',line)
        if c:
            x=1
        if x==2:
            yes=re.search('YES',line)
            no=re.search('NO',line)
            if yes:
                converge.append('YES')
            if no:
                converge.append('NO')
            u=re.search('-------------------------',line)
            if u:
                x=0
        if x==1:
            u=re.search('-------------------------',line)
            if u:
                x=2
        x3=re.search('aborting the run',line)
        if x3:
            s1=4
    if statusonly==True:
        if s2==1:
            return 2
        else:
            return s1
    if figureonly==False:
        print('{0} turns have been calculated.'.format(term))
        status_judge(s1,s2,converge)
        statusfig(energy,choose)
    os.chdir(path)
    return energy

class status:
    """
    A method to see the current process of the ORCA optimisation, including convergence situation and relevant energy chart.
    status(file, loc='./')
    file: File Name.
    loc: File Location. The default is your current location.
    You can get the further information by .status and .figures.
    """

    def __init__(self,file,loc='./'):
        self.file=file
        self.loc=loc
    
    def status(self,choose=0,figureonly=False,statusonly=False):
        """
    A method to see the current process of the ORCA optimisation, including convergence situation and relevant energy chart.
    status(choose=0, figureonly=False, statusonly=False,last=-1)
    choose: See the latest energy process. e.g. When choose=5, it will show the last five energy data on the chart.
    figureonly: If 'figureonly' is True, the code will show no status information. The default is False.
    statusonly: If 'statusonly' is True, the code will output the current process status only, not energy. 
                in output s, 0 means the system is optimising, 1 means the system completes the optimisation and starts frequency calculation, 2 means the process is completed, 4 means the process is aborted.
                The default is False.
                TIPS: Kill ORCA process with problems other than ORCA will get the output of 0.
    Normally, the output of this function is the energy change in optimisation, unless statusonly=True.
    Example 1:
        Input:
            from MCPoly import status
            a=status('Et','/Molecule').status()
        
        Output:
            Et.out
            8 turns have been calculated.
            The Thermodynamics Calculation was finished.
            <chart from matplotlib>
            
    Example 2:
        Input:
            from MCPoly import status
            a=status('Molecule 1','/Molecule').status(figureonly=True)
            print(a)
        
        Output:
            [-1729.508137253801, -1729.505855679262, -1729.501743481437, -1729.499301095603, -1729.499960765849, 
             -1729.499875258488, -1729.500097690162, -1729.500368245231, -1729.500825118936, -1729.501319993802, 
             -1729.503640287060, -1729.505369082598, -1729.507234351098, -1729.508676755388, ...]
    
    Example 3:
        Input:
            from MCPoly import status
            a=status('SI5_0.5','../Sulphur').status(figureonly=True,statusonly=True)
            print(a)
        
        Output:
            2
        """

        return normalstatus(self.loc,self.file,choose,figureonly,statusonly)
    
    def figure(self,num=0,width=300,height=300):
        """
    A method to see the current geometry structure of the ORCA optimisation, powered by py3Dmol and rdkit.
    TIPS: Make sure your _trj.xyz file is in the document with .out file, or there will be NoFileFoundError!!!
    
    figure(num=0, width=300, height=300)
    num: The step of your convergence.
    width, height: The size of your 3D geometry molecule strcuture. Default: 300x300.
        """
            
        figure=normalstatus(self.loc,self.file,figureonly=True)
        file=self.file
        try:
            path = os.getcwd()
            os.chdir(self.loc)
            multicreatemol(file, num, width, height)
            os.chdir(path)
        except:
            raise ValueError("The index 'num' is out of range.")
        
    def figuretraj(self,num=0,width=300,height=300):
        """
    A method to see the current geometry structure and optimization trajectory of the ORCA optimisation, powered by py3Dmol and ipywidgets package.
    TIPS: Make sure your _trj.xyz file is in the document with .out file, or there will be FileNotFoundError!!!
    
    figurestatus(num=0, width=300, height=300)
    num: The step of your convergence.
    width, height: The size of your 3D geometry molecule strcuture. Default: 300x300.
    After forming the 3D geometry molecule strcuture, you can scroll to see other structures of relevant molecules.
        """
            
        figure=normalstatus(self.loc,self.file,figureonly=True)
        file=self.file
        def turn(num):
            try:
                path = os.getcwd()
                os.chdir(self.loc)
                multicreatemol(file, num, width, height)
                os.chdir(path)
            except:
                if num==len(figure)-1:
                    print("The last step hasn't been optimised yet.")
        
        interact(turn,num=iw.IntSlider(min=0,max=len(figure)-1,step=1,value=num))