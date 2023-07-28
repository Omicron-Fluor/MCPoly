import re
import os
import sys
import py3Dmol
import matplotlib.pyplot as plt
from ipywidgets import interact
import ipywidgets as iw
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDetermineBonds

def status_judge(s1,s2,converge):
    if s1==0 or s1==4 or s1==41:
        if s1==0:
            print('The Optimization is processing.')
        elif s1==4 and s2==0:
            print('The Optimization is aborted.')
        elif s1==41:
            print('The Optimization is aborted.\n The structure includes 180 or 0 degrees angles.')
        print('Last time:')
        print('Energy change: {0}'.format(converge[-5]))
        print('RMS gradient: {0}             MAX gradient: {0}'.format(converge[-4],converge[-3]))
        print('RMS step: {0}                 MAX step: {0}'.format(converge[-2],converge[-1]))
    if s1!=0 and s2<0:
        print('The Optimization was finished.')
        if s1==4 or s1>=40:
            print('The frequency calculation is aborted.')
            if s1==41:
                print('There are some 0 or 180 degrees angles.\n Try to rebuild the system.')
            return None
        if s2==-1:
            print('The CP-SCF equations are forming.')
        elif s2==-2:
            print('The CP-SCF equations are solving.')
        elif s2==-3:
            print('The Thermodynamics Calculation was processing.')
    if s2==1:
        print('The Thermodynamics Calculation was finished.')
    if s2==9 and s1==6:
        print('The Thermodynamics Calculation was finished, but the optimization did not converge but reached the maximum number of optimization cycles.')
        
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

def multicreatemol(filename,turn,width,height,MD=False):
    if MD==True:
        try:
            ind = open('{0}_dump.xyz'.format(filename),'r+')
        except:
            raise FileNotFoundError('Perhaps your _dump.xyz file has some errors. Please check the location of your _dump.xyz file and then modify your file.')
    else:
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

def normalstatus(loc,file,choose=0,figureonly=False,statusonly=False,MD=False):
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
    if MD==False:
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
                if s1==6:
                    s2=9
                else:
                    s2=0
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
            x3=re.search('ANGLE IS APPROACHING 180',line)
            if x3:
                s1=41
                break
            x4=re.search('did not converge but reached',line)
            if x4:
                s1=6
    else:
        s2=1
        ii=0
        gap=0
        for line in f:
            x0=re.search('Dump Position Stride [0-9]+',line)
            if x0:
                y0=re.search(r'[0-9]+',x0.group(0))
                gap=eval(y0.group(0))
            x1=re.findall('\-\|',line)
            if len(x1)>3:
                if ii==0:
                    ii=1
                else:
                    break
            elif ii==1:
                alls=re.findall(r'\-*[0-9]+\.*[0-9]*',line)
                if alls[0]=='0':
                    energy.append(eval(alls[-1]))
                else:
                    if eval(alls[0])%gap==0:
                        energy.append(eval(alls[-3]))
    f.close()
    if statusonly==True:
        if s1==1:
            return 2
        else:
            return s1
    if figureonly==False and MD==False:
        print('{0} turns have been calculated.'.format(term))
        status_judge(s1,s2,converge)
        statusfig(energy,choose)
    elif figureonly==False and MD==True:
        print('{0} turns have been calculated.'.format(len(energy)))
        statusfig(energy,choose)
    os.chdir(path)
    return energy

def thermo(loc, file, keyword):
    f=open(file+'.out','r')
    for line in f:
        if keyword=='Gibbs':
            a=re.search('Final Gibbs free energy', line)
            if a:
                b=re.search(r'-[0-9]+\.[0-9]+', line)
                return '{0:.6f}'.format(eval(b.group(0)))
        if keyword=='Enthalpy':
            a=re.search('Total enthalpy', line)
            if a:
                b=re.search(r'-[0-9]+\.[0-9]+', line)
                return '{0:.6f}'.format(eval(b.group(0)))
        if keyword=='Entropy':
            a=re.search('Total entropy correction', line)
            if a:
                b=re.findall(r'-[0-9]+\.[0-9]+', line)
                return '{0:.6f}'.format(eval(b[0]))
    raise Exception("'{0}.out' may not have keyword freq, or it's aborted. Please check your input file and recalculate it.".format(file))

def charge(loc, file, num, keyword):
    f=open(file+'.out','r')
    i=0
    ac=[None]*200
    for line in f:
        if keyword=='Mulliken':
            d=re.search('Sum of atomic charges',line)
            if d:
                i=0
            a=re.search('MULLIKEN ATOMIC CHARGES',line)
            if a:
                i=1
            if i==1:
                b=re.search(r'-?[0-9]+\.[0-9]+', line)
                if b:
                    n=re.findall(r'[0-9]+', line)
                    ac[eval(n[0])]=eval(b.group(0))
        elif keyword=='Loewdin':
            a=re.search('LOEWDIN ATOMIC CHARGES',line)
            if a:
                i=1
            if i==1:
                b=re.search(r'-?[0-9]+\.[0-9]+', line)
                if b:
                    n=re.findall(r'[0-9]+', line)
                    ac[eval(n[0])]=eval(b.group(0))
            d=re.search('LOEWDIN REDUCED ORBITAL CHARGES',line)
            if d:
                i=0
    return ac[num]

def atoms_all(loc, file):
    f=open(file+'.out','r')
    i=0
    n=0
    atoms=[['NOT','99999','99999','99999']]*200
    elements=['NOT']*200
    xyzs=[['99999','99999','99999']]*200
    for line in f:
        d=re.search('CARTESIAN COORDINATES \(A\.U\.\)',line)
        if d:
            i=0
            n=0
        if i==1:
            b=re.search(r'[A-Z][a-z]*', line)
            if b:
                xyz=re.findall(r'-?[0-9]+\.[0-9]+', line)
                try:
                    elements[n]=b.group(0)
                except:
                    elements.append(b.group(0))
                try:
                    xyzs[n]=xyz
                except:
                    xyzs.append(xyz)
                try:
                    atoms[n]=[b.group(0),eval(xyz[0]),eval(xyz[1]),eval(xyz[2])]
                except:
                    atoms.append([b.group(0),eval(xyz[0]),eval(xyz[1]),eval(xyz[2])])
                n=n+1
        a=re.search('CARTESIAN COORDINATES \(ANGSTROEM\)',line)
        if a:
            i=1
    while 1:
        try:
            atoms.remove(['NOT','99999','99999','99999'])
            elements.remove('NOT')
            xyzs.remove(['99999','99999','99999'])
        except:
            break
    f.close()
    return atoms

def atoms_all_xyz(loc, file):
    f=open(file+'.xyz','r')
    atoms=[['NOT','99999','99999','99999']]*200
    elements=['NOT']*200
    xyzs=[['99999','99999','99999']]*200
    i=0
    for line in f:
        b=re.findall(r'-?[0-9]+\.[0-9]+', line)
        if len(b)==3:
            ele=re.search(r'[A-Z][a-z]?', line)
            try:
                atoms[i]=[ele.group(0),eval(b[0]),eval(b[1]),eval(b[2])]
            except:
                atoms.append([ele.group(0),eval(b[0]),eval(b[1]),eval(b[2])])
        i=i+1
    while 1:
        try:
            atoms.remove(['NOT','99999','99999','99999'])
            elements.remove('NOT')
            xyzs.remove(['99999','99999','99999'])
        except:
            break
    f.close()
    return atoms

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
    
    def energy(self):
        """
    A method to see the final(or latest) Energy of a system.
    energy()
        """
        E=normalstatus(self.loc,self.file,figureonly=True)
        return eval('{0:.6f}'.format(E[-1]))

    def converged_energy(self,unit='Eh'):
        """
    A method to see the Converged Energy of a system.
    converged_energy(unit)
    unit: Energy unit, including 'Eh' and 'kcal/mol'. The default is 'Eh' (Hartree).
        """
        E=normalstatus(self.loc,self.file,figureonly=True)
        state=normalstatus(self.loc,self.file,figureonly=True,statusonly=True)
        if state!=2:
            raise Exception("This optimisation of '{0}.out' is failed. No converged energy.".format(self.file))
        else:
            x=eval('{0:.6f}'.format(E[-1]))
            if unit=='Eh':
                return x
            elif unit=='kcal/mol':
                return x*627.509
    
    def gibbs(self,unit='Eh'):
        """
    A method to see the Gibbs Free Energy of a system.
    gibbs(unit)
    unit: Energy unit, including 'Eh' and 'kcal/mol'. The default is 'Eh' (Hartree).
    TIPS: Make sure the optimisation includes key word 'freq', which means frequency calculation.
        """
        x=eval(thermo(self.loc,self.file,keyword='Gibbs'))
        if unit=='Eh':
            return x
        elif unit=='kcal/mol':
            return x*627.509
    
    def enthalpy(self,unit='Eh'):
        """
    A method to see the total enthalpy of a system.
    enthalpy(unit)
    unit: Energy unit, including 'Eh' and 'kcal/mol'. The default is 'Eh' (Hartree).
    TIPS: Make sure the optimisation includes key word 'freq', which means frequency calculation.
        """
        x=eval(thermo(self.loc,self.file,keyword='Enthalpy'))
        if unit=='Eh':
            return x
        elif unit=='kcal/mol':
            return x*627.509
    
    def entropy_correction(self,unit='Eh'):
        """
    A method to see the total entropy correction of a system.
    enthalpy(unit)
    unit: Energy unit, including 'Eh' and 'kcal/mol'. The default is 'Eh' (Hartree).
    TIPS: Make sure the optimisation includes key word 'freq', which means frequency calculation.
        """
        x=eval(thermo(self.loc,self.file,keyword='Entropy'))
        if unit=='Eh':
            return x
        elif unit=='kcal/mol':
            return x*627.509
    
    def mulliken_charge(self,num):
        """
    A method to see the Mulliken Charge of an atom.
    mulliken_charge(num)
    num: Atom number.
    TIPS: The atom starts from 0. Don't surpass the maximum, or it will cause error.
        """
        x=charge(self.loc,self.file,num,keyword='Mulliken')
        return x
    
    def loewdin_charge(self,num):
        """
    A method to see the Loewdin Charge of an atom.
    loewdin_charge(num)
    num: Atom number.
    TIPS: The atom starts from 0. Don't surpass the maximum, or it will cause error.
        """
        x=charge(self.loc,self.file,num,keyword='Loewdin')
        return x

    def atom(self):
        """
    A method to find out the atom information of the system.
    atom()
    
    TIPS: In the following list, a piece of atom information is like that:
        ['C', 4.11199532065866, 1.93907233706568, -1.49149364116961]
        The first one is the element type, the following three numbers are location in x,y,z-axis.
        """
        try:
            return atoms_all(self.loc,self.file)
        except:
            return atoms_all_xyz(self.loc,self.file)
    
    def mass(self):
        """
    A method to find out the mass of all particles in system.
    mass()
        """
        M=0
        try:
            x=atoms_all(self.loc,self.file)
        except:
            x=atoms_all_xyz(self.loc,self.file)
        for i in x:
            if i[0]=='H':
                M=M+1.008
            elif i[0]=='C':
                M=M+12.011
            elif i[0]=='O':
                M=M+15.999
            elif i[0]=='Li':
                M=M+6.94
            elif i[0]=='Be':
                M=M+9.0122
            elif i[0]=='B':
                M=M+10.81
            elif i[0]=='N':
                M=M+14.007
            elif i[0]=='F':
                M=M+18.998
            elif i[0]=='Na':
                M=M+22.990
            elif i[0]=='Mg':
                M=M+24.305
            elif i[0]=='Al':
                M=M+26.982
            elif i[0]=='Si':
                M=M+28.085
            elif i[0]=='P':
                M=M+30.974
            elif i[0]=='S':
                M=M+32.06
            elif i[0]=='Cl':
                M=M+35.45
            elif i[0]=='Br':
                M=M+79.904
            elif i[0]=='I':
                M=M+126.90
            elif i[0]=='K':
                M=M+39.098
            elif i[0]=='Ca':
                M=M+40.078
            elif i[0]=='Ti':
                M=M+47.867
            elif i[0]=='Cr':
                M=M+51.996
            elif i[0]=='Mn':
                M=M+54.938
            elif i[0]=='Fe':
                M=M+55.845
            elif i[0]=='Co':
                M=M+58.933
            elif i[0]=='Ni':
                M=M+58.693
            elif i[0]=='Cu':
                M=M+63.546
            elif i[0]=='Zn':
                M=M+65.38
            elif i[0]=='Ru':
                M=M+101.07
            elif i[0]=='Rh':
                M=M+102.91
            elif i[0]=='Ag':
                M=M+107.87
            elif i[0]=='Ir':
                M=M+192.22
            elif i[0]=='Pt':
                M=M+195.08
            elif i[0]=='Au':
                M=M+196.97
            elif i[0]=='Hg':
                M=M+200.59
        return M
    
    def atom_num(self):
        """
    A method to find out the atom number of the system.
    atom_num()
        """
        return len(atoms_all(self.loc,self.file))
    
    def figure(self,num=0,width=300,height=300,MD=False,save=''):
        """
    A method to see the current geometry structure of the ORCA optimisation, powered by py3Dmol and rdkit.
    TIPS: Make sure your _trj.xyz file is in the document with .out file, or there will be NoFileFoundError!!!
    
    figure(num=0, width=300, height=300,MD=False,save='')
    num: The step of your convergence.
    width, height: The size of your 3D geometry molecule strcuture. Default: 300x300.
    MD: See the single status of Molecular Dynamics (MD). The default is False.
    save: The name of the single status XYZ file you want to save. The default is not saving the single status.
        """
            
        figure=normalstatus(self.loc,self.file,figureonly=True,MD=MD)
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
    num: The step of your convergence. The default is the origin structure.
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
    
    def figuremd(self,num=0,width=300,height=300):
        """
    A method to see the current geometry structure and optimization trajectory of the ORCA Molecular Dynamics, powered by py3Dmol and ipywidgets package.
    TIPS: Make sure your _dump.xyz file is in the document with .out file, or there will be FileNotFoundError!!!
    
    figuremd(num=0, width=300, height=300)
    num: The step of your convergence. The default is the origin structure.
    width, height: The size of your 3D geometry molecule strcuture. Default: 300x300.
    After forming the 3D geometry molecule strcuture, you can scroll to see other structures of relevant molecules.
        """
            
        figure=normalstatus(self.loc,self.file,figureonly=True,MD=True)
        file=self.file
        def turn(num):
            try:
                path = os.getcwd()
                os.chdir(self.loc)
                multicreatemol(file, num, width, height, MD=True)
                os.chdir(path)
            except:
                if num==len(figure)-1:
                    print("The last step hasn't been optimised yet.")
        interact(turn,num=iw.IntSlider(min=0,max=len(figure),step=1,value=num))
    
    def figurecharge(self,width=300,height=300):
        """
    A method to see the Mulliken Charge of each atom base on 3D geometry structure of the system.
    TIPS: Make sure your .xyz file is in the document with .out file, or there will be FileNotFoundError!!!
    
    width, height: The size of your 3D geometry molecule strcuture. Default: 300x300.
    After forming the 3D geometry molecule strcuture, you can scroll to see the charge of the relevant atom.
        """
        output=iw.Output()
        mydir = os.path.dirname( __file__ )
        orcadir = os.path.join(mydir, '..', 'orcaset')
        sys.path.append(orcadir)
        from view3dchoose import view3dchoose
        
        chargetype=iw.Dropdown(description='Atomic Charge Type',options=['Mulliken','Loewdin'])
        num=iw.IntSlider(value=0,min=0,max=len(atoms_all(self.loc,self.file))-1)
        
        def GUI(chargetype, num):
            output.clear_output()
            with output:
                alls=atoms_all(self.loc,self.file)
                view3dchoose(self.file,self.loc,choose=[num],width=width,height=height)
                print('Atom: {0}, Element: {1}, Charge: {2}'.format(num,alls[num][0],charge(self.loc,self.file,num,chargetype))) 
            
        geoout=iw.interactive_output(GUI, {'chargetype': chargetype, 'num': num})
        display(iw.VBox([chargetype,num,output]))
