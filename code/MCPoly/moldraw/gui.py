# import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.visualize import view
from MCPoly import view3d
from MCPoly.moldraw import molecule
from MCPoly.orcaset import view3dchoose
import ipywidgets as widgets
from IPython.display import display
from rdkit import Chem
from rdkit.Chem import AllChem
from ipywidgets import interact
import os
import re
import py3Dmol

def smi2conf(smiles):
    '''Convert SMILES to rdkit.Mol with 3D coordinates'''
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        return mol
    else:
        return None

def MolTo3DView(mol,size=(400,550),surface=False,opacity=0.5):
    """Draw molecule in 3D (Not Origin)
    
    Args:
    ----
        mol: rdMol, molecule to show
        size: tuple(int, int), canvas size
        style: str, type of drawing molecule
               style can be 'line', 'stick', 'sphere', 'cartoon'
        surface, bool, display SAS
        opacity, float, opacity of surface, range 0.0-1.0
    Return:
    ----
        viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
    """
    mblock = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=size[0],height=size[1])
    viewer.addModel(mblock,'mol')
    viewer.setStyle({'sphere':{'radius':0.4},'stick':{'radius':0.1}})
    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    Chem.MolToXYZFile(mol,'Mole_Untitled.xyz')
    return viewer

def gui():
    """
        The method to draw molecule and save it in the folder.
        You will know how to do as long as you get into the GUI platform.
    
        TIPS: 1. This drawing function is under construction. There might be some errors in the programme.
              2. 'shuffle', 'less than 5' are used to modify some geometry structure.
              3. The visualisation will also risk the flexibility of the task, to make you programme more flexible, please directly use the code.
    """
    
    def outputmol1(button):
        output.clear_output()
        output2.clear_output()
        placeoutput.clear_output()
        strictoutput1.clear_output()
        strictoutput2.clear_output()
        strictoutput3.clear_output()
        angleoutput.clear_output()
        shuffleoutput.clear_output()
        strictadd2.value=False
        strictadd3.value=False
        with output:
            i.max=atoms.atomnum()-1
            if button.description=='pyrrole':
                place.max=5
                with placeoutput:
                    display(place)
            if button.description=='Py':
                place.max=6
                with placoutput:
                    display(place)
            if button.description=='cyclo3':
                strictnum1.max=3
            elif button.description=='cyclo5':
                strictnum1.max=5
                strictnum2.max=5
            elif button.description=='cyclo6':
                strictnum1.max=6
                strictnum2.max=6
                strictnum3.max=6
            if button.description in ['C2H3','C=NH','Ph','Py','pyrrole','cyclo5','cyclo6','Bpin','bicyclo[2.2.1]','bicyclo[2.2.2]']:
                with shuffleoutput:
                    display(shuffle)
            if button.description!='cyclo3' and button.description!='cyclo5' and button.description!='cyclo6':
                pass
            elif button.description=='cyclo3':
                with strictoutput1:
                    display(widgets.VBox([strictnum1,widgets.HBox([strictele1])]))
            else:
                with strictoutput1:
                    display(widgets.VBox([strictnum1,widgets.HBox([strictele1,strictadd2])]))
            if button.description=='OH' or button.description=='OMe' or button.description=='SH' or button.description=='SMe':
                with angleoutput:
                    display(widgets.HBox([angle,angleavail]))
            buttonname.value=button.description
            display(widgets.VBox([i,state,placeoutput,strictoutput1,strictoutput2,strictoutput3,shuffleoutput,angleoutput]))
            display(widgets.HBox([preview,settle]))
    
    def outputmol2(button):
        output.clear_output()
        output2.clear_output()
        placeoutput.clear_output()
        strictoutput1.clear_output()
        strictoutput2.clear_output()
        strictoutput3.clear_output()
        angleoutput.clear_output()
        shuffleoutput.clear_output()
        strictadd2.value=False
        strictadd3.value=False
        with output:
            i1.max=atoms.atomnum()-1
            i2.max=atoms.atomnum()-1
            state1.options=[('1',1),('2',2),('6',6)]
            state2.options=[('1',1),('2',2),('6',6)]
            if button.description=='bind3':
                strictnum1.max=1
            elif button.description=='bind5' or button.description=='bindconj5':
                strictnum1.max=3
                strictnum2.max=3
            elif button.description=='bind6' or button.description=='bindconj6':
                strictnum1.max=4
                strictnum2.max=4
            if button.description not in ['bind2','bindconj5','bindconj6']:
                with shuffleoutput:
                    display(shuffle)
            if button.description=='bind2':
                pass
            elif button.description=='bind3':
                with strictoutput1:
                    display(widgets.VBox([strictnum1,widgets.HBox([strictele1])]))
            else:
                with strictoutput1:
                    display(widgets.VBox([strictnum1,widgets.HBox([strictele1,strictadd2])]))
            buttonname.value=button.description
            display(widgets.VBox([widgets.HBox([i1,state1]),widgets.HBox([i2,state2]),strictoutput1,strictoutput2,side,shuffle,less5,stable]))
            display(widgets.HBox([preview,settle]))
            
    def outputmol3(button):
        output.clear_output()
        with output:
            i.max=atoms.atomnum()-1
            state1.options=[('1',1),('2',2),('6',6),('7',7),('8',8),('82',82)]
            state2.options=[('1',1),('2',2),('6',6),('7',7),('8',8),('82',82)]
            buttonname.value=button.description
            display(widgets.VBox([i,state1,state2]))
            display(widgets.HBox([preview,settle]))
        
    def clear(button):
        output.clear_output()
        output2.clear_output()

    def openf(button):
        os.system('cp {0}.xyz Mole_Untitled.xyz'.format(file.value))
        file.value='Mole_Untitled'
        with output3:
            view3d(file.value,loc.value,width=400,height=500)
            display(widgets.HBox([savename,savefile]))
        output0.clear_output()
        global atoms
        atoms=molecule(file.value,loc.value)
        with output0:
            display(widgets.VBox([allchoose,stack]))

    def smartshow(smart):
        os.system('rm Mole_Untitled.xyz')
        output3.clear_output()
        with output3:
            try:
                conf = smi2conf(smart)
                return MolTo3DView(conf).show()
            except:
                return None
    
    def smartset(button):
        output3.clear_output()
        with output3:
            view3d('Mole_Untitled',loc.value,width=400,height=500)
            display(widgets.HBox([savename,savefile]))
        output0.clear_output()
        file.value='Mole_Untitled'
        global atoms
        atoms=molecule('Mole_Untitled',loc.value)
        with output0:
            display(widgets.VBox([allchoose,stack]))
    
    def saveset(button):
        output0.clear_output()
        output.clear_output()
        output3.clear_output()
        atoms.xyzwrite(savename.value)
        os.system('rm Mole_Untitled.xyz')
        file.value=''
        with output0:
            display(widgets.VBox([loc,widgets.HBox([file,openfile]),widgets.Label('OR',style=dict(font_weight='bold',text_decoration='underline')),\
                          widgets.HBox([smart,smartbutton])]))
    
    def viewset(button):
        global atoms
        output3.clear_output()
        atoms_b=atoms
        if choose.value=='Substrates':
            if buttonname.value=='cyclo3':
                if strictele1.value[0]=='X':
                    atoms_b.sub(buttonname.value,i.value,state.value,place.value,shuffle=shuffle.value,angle=angle.value)
                else:
                    atoms_b.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1},shuffle=shuffle.value,angle=angle.value)
            elif buttonname.value=='cyclo5':
                if strictele1.value[0]=='X':
                    atoms_b.sub(buttonname.value,i.value,state.value,place.value,shuffle=shuffle.value,angle=angle.value)
                elif strictadd2.value==False:
                    atoms_b.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1},shuffle=shuffle.value,angle=angle.value)
                else:
                    atoms_b.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1,strictnum2:strictele2},shuffle=shuffle.value,angle=angle.value)
            elif buttonname.value=='cyclo6':
                if strictele1.value[0]=='X':
                    atoms_b.sub(buttonname.value,i.value,state.value,place.value,shuffle=shuffle.value,angle=angle.value)
                elif strictadd2.value==False:
                    atoms_b.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1},shuffle=shuffle.value,angle=angle.value)
                elif strictadd3.value==False:
                    atoms_b.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1,strictnum2:strictele2},shuffle=shuffle.value,angle=angle.value)
                else:
                    atoms_b.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1,strictnum2:strictele2,strictnum3:strictele3},shuffle=shuffle.value,angle=angle.value)
            else:
                atoms_b.sub(buttonname.value,i.value,state.value,place.value,shuffle=shuffle.value,angle=angle.value)
        elif choose.value=='Ringbinders':
            if buttonname.value=='bind3':
                atoms_b.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,strict={strictnum1:strictele1},side=side.value,shuffle=shuffle.value, less5=less5.value, stable=stable.value)
            elif buttonname.value in ['SmallChanges','bindconj5','bindconj6']:
                if strictele1.value[0]=='X':
                    atoms_b.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,side=side.value,shuffle=shuffle.value,less5=less5.value,stable=stable.value)
                elif strictadd2.value==False:
                    atoms_b.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,strict={strictnum1:strictele1},side=side.value,shuffle=shuffle.value,less5=less5.value,stable=stable.value)
                else:
                    atoms_b.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,strict={strictnum1:strictele1,strictnum2:strictele2},side=side.value,shuffle=shuffle.value,less5=less5.value,stable=stable.value)
            else:
                atoms_b.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,side=side.value,shuffle=shuffle.value,less5=less5.value,stable=stable.value)
        elif choose.value=='SmallChanges':
            atoms_b.C(buttonname.value,i.value,state=[state1.value,state2.value])
        atoms_b.xyzwrite('Mole_Untitled_beta')
        with output3:
            view3d('Mole_Untitled_beta',loc.value,width=400,height=500)
            display(widgets.HBox([savename,savefile]))
        file.value='Mole_Untitled'
        os.system('rm Mole_Untitled_beta.xyz')
        atoms=molecule('Mole_Untitled',loc.value)

    def settlechange(button):
        global atoms
        output3.clear_output()
        if choose.value=='Substrates':
            if buttonname.value=='cyclo3':
                if strictele1.value[0]=='X':
                    atoms.sub(buttonname.value,i.value,state.value,place.value,shuffle=shuffle.value,angle=angle.value)
                else:
                    atoms.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1},shuffle=shuffle.value,angle=angle.value)
            elif buttonname.value=='cyclo5':
                if strictele1.value[0]=='X':
                    atoms.sub(buttonname.value,i.value,state.value,place.value,shuffle=shuffle.value,angle=angle.value)
                elif strictadd2.value==False:
                    atoms.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1},shuffle=shuffle.value,angle=angle.value)
                else:
                    atoms.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1,strictnum2:strictele2},shuffle=shuffle.value,angle=angle.value)
            elif buttonname.value=='cyclo6':
                if strictele1.value[0]=='X':
                    atoms.sub(buttonname.value,i.value,state.value,place.value,shuffle=shuffle.value,angle=angle.value)
                elif strictadd2.value==False:
                    atoms.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1},shuffle=shuffle.value,angle=angle.value)
                elif strictadd3.value==False:
                    atoms.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1,strictnum2:strictele2},shuffle=shuffle.value,angle=angle.value)
                else:
                    atoms.sub(buttonname.value,i.value,state.value,place.value,strict={strictnum1:strictele1,strictnum2:strictele2,strictnum3:strictele3},shuffle=shuffle.value,angle=angle.value)
            else:
                atoms.sub(buttonname.value,i.value,state.value,place.value,shuffle=shuffle.value,angle=angle.value)
        elif choose.value=='Ringbinders':
            if buttonname.value=='bind3':
                atoms.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,strict={strictnum1:strictele1},side=side.value,shuffle=shuffle.value, less5=less5.value, stable=stable.value)
            elif buttonname.value in ['SmallChanges','bindconj5','bindconj6']:
                if strictele1.value[0]=='X':
                    atoms.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,side=side.value,shuffle=shuffle.value,less5=less5.value,stable=stable.value)
                elif strictadd2.value==False:
                    atoms.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,strict={strictnum1:strictele1},side=side.value,shuffle=shuffle.value,less5=less5.value,stable=stable.value)
                else:
                    atoms.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,strict={strictnum1:strictele1,strictnum2:strictele2},side=side.value,shuffle=shuffle.value,less5=less5.value,stable=stable.value)
            else:
                atoms.bind(buttonname.value,[i1.value,i2.value],state=[state1.value,state2.value],place=place.value,side=side.value,shuffle=shuffle.value,less5=less5.value,stable=stable.value)
        elif choose.value=='SmallChanges':
            atoms.C(buttonname.value,i.value,state=[state1.value,state2.value])
        atoms.xyzwrite('Mole_Untitled')
        atoms=molecule('Mole_Untitled',loc.value)
        file.value='Mole_Untitled'
        with output3:
            view3d('Mole_Untitled',loc.value,width=400,height=500)
            display(widgets.HBox([savename,savefile]))
    
    def strictout2(strictadd2):
        strictoutput2.clear_output()
        if strictadd2==True:
            with strictoutput2:
                if buttonname.value=='cyclo5' or buttonname.value=='bind5' or buttonname.value=='bind6' or buttonname.value=='bindconj5' or buttonname.value=='bindconj6':
                    display(widgets.VBox([strictnum2,strictele2]))
                else:
                    display(widgets.VBox([strictnum2,widgets.HBox([strictele2,strictadd3])]))
    
    def strictout3(strictadd3):
        strictoutput3.clear_output()
        if strictadd3==True:
            with strictoutput3:
                display(widgets.VBox([strictnum3,strictele3]))
    
    def angleavailout(angleavail):
        if angleavail==True:
            angle.disabled=True
        else:
            angle.disabled=False

    def chooseshow(i,i1,i2,choose):
        output3.clear_output()
        with output3:
            try:
                if choose=='Substrates':
                    view3dchoose(file.value,loc.value,choose=[i],width=400,height=500)
                elif choose=='Ringbinders':
                    view3dchoose(file.value,loc.value,choose=[i1,i2],width=400,height=500)
                elif choose=='SmallChanges':
                    view3dchoose(file.value,loc.value,choose=[i],width=400,height=500)
                display(widgets.HBox([savename,savefile]))
            except:
                pass
    
    def default(strictele1):
        if strictele1=='X(Default)':
            strictnum1.value=1
            strictadd2.disabled=True
        else:
            strictadd2.disabled=False
    
    output0 = widgets.Output()
    output = widgets.Output()
    output2 = widgets.Output()
    output3 = widgets.Output()
    strictoutput1 = widgets.Output()
    strictoutput2 = widgets.Output()
    strictoutput3 = widgets.Output()
    angleoutput = widgets.Output()
    placeoutput = widgets.Output()
    shuffleoutput = widgets.Output()
    
    options=[]
    loc=widgets.Text(description='Location:',value='./')
    for path in os.listdir(loc.value):
        if os.path.isfile(os.path.join(loc.value, path)):
            a=re.search('.xyz', path)
            if a:
                options.append(path[:-4])
    file=widgets.Combobox(placeholder='Input the file name',options=options,
             description='xyz File:')
    
    openfile=widgets.Button(description='Set')
    openfile.on_click(openf)
    
    smart=widgets.Text(description='SMART:',placeholder='Input SMART structure')
    smartbutton=widgets.Button(description='Set')
    smartbutton.on_click(smartset)
    
    i=widgets.BoundedIntText(description='Atom Number:',min=0,style={'description_width': 'initial'})
    state=widgets.Dropdown(description='Environment:',options=[('1',1),('2',2),('4',4),('6',6),('7',7),('8',8),('16',16),('82',82),('162',162),('164',164)],value=1,style={'description_width': 'initial'})
    place=widgets.BoundedIntText(description='Connect Place:',min=0,value=0,style={'description_width': 'initial'})
    strictnum1=widgets.BoundedIntText(description='Hetero Place 1:',min=1,value=1,style={'description_width': 'initial'})
    strictele1=widgets.Dropdown(description='Hetero Elements 1:',options=['N','O','S','X(Default)'],value='X(Default)',style={'description_width': 'initial'})
    #strictdefault=widgets.ToggleButton(description='default',value=True)
    strictnum2=widgets.BoundedIntText(description='Hetero Place 2:',min=3,value=3,style={'description_width': 'initial'})
    strictele2=widgets.Dropdown(description='Hetero Elements 2:',options=['N','O','S'],style={'description_width': 'initial'})
    strictadd2=widgets.ToggleButton(description='Add',value=False)
    strictnum3=widgets.BoundedIntText(description='Hetero Place 3:',min=5,value=5,style={'description_width': 'initial'})
    strictele3=widgets.Dropdown(description='Hetero Elements 3:',options=['N','O','S'],style={'description_width': 'initial'})
    strictadd3=widgets.ToggleButton(description='Add',value=False)
    i1=widgets.BoundedIntText(description='Atom Number 1:',min=0,style={'description_width': 'initial'})
    state1=widgets.Dropdown(description='Environment:',options=[('1',1),('2',2),('6',6)],value=1,style={'description_width': 'initial'})
    i2=widgets.BoundedIntText(description='Atom Number 2:',min=0,style={'description_width': 'initial'})
    state2=widgets.Dropdown(description='Environment:',options=[('1',1),('2',2),('6',6)],value=1,style={'description_width': 'initial'})
    shuffle=widgets.BoundedIntText(description='Shuffle:',min=0,max=3,value=0,style={'description_width': 'initial'})
    angle=widgets.BoundedFloatText(min=0,max=360,step=0.1,description='Angle:',style={'description_width': 'initial'})
    angleavail=widgets.ToggleButton(value=False,description='Default',style={'description_width': 'initial'})
    side=widgets.BoundedIntText(description='Conformer:',min=0,max=3,value=0,style={'description_width': 'initial'})
    shuffle=widgets.BoundedIntText(description='Shuffle:',min=0,value=0,style={'description_width': 'initial'})
    less5=widgets.Checkbox(description='Connect with small rings (member<=5)?',style={'description_width': 'initial'})
    stable=widgets.Checkbox(description='Most Stable structure?',style={'description_width': 'initial'})
    preview=widgets.Button(description='View')
    preview.on_click(viewset)
    settle=widgets.Button(description='Set')
    settle.on_click(settlechange)
    savename=widgets.Text(description='Save Name:',value='Mole_Untitled')
    savefile=widgets.Button(description='Save')
    savefile.on_click(saveset)
    buttonname=widgets.Text()
    
    #file=widgets.Button(description='Set')
    
    strict2=widgets.interactive_output(strictout2, {'strictadd2': strictadd2})
    strict3=widgets.interactive_output(strictout3, {'strictadd3': strictadd3})
    angleout=widgets.interactive_output(angleavailout, {'angleavail': angleavail})
    smartout=widgets.interactive_output(smartshow, {'smart': smart})
    
    substrates=['F','Cl','Br','I','OH','OMe','SH','SMe','NH2','NMe2','Me','Et','n-Pr','i-Pr','n-Bu','t-Bu','CF3','CCl3','CBr3','CI3','NO2','C2H3','C=NH',
          'CHO','COOH','COMe','CN','Ph','pyrrole','Py','cyclo3','cyclo5','cyclo6','SOMe','SO2Me','SO3H','Bpin','bicyclo[2.2.1]','bicyclo[2.2.2]']
    buttons1=[]
    for x in substrates:
        button=widgets.Button(layout=widgets.Layout(height='30px', min_width='120px'),description=x)
        buttons1.append(button)
        button.on_click(outputmol1)
    Substrates=widgets.HBox(children=[buttons1[x] for x in range(len(substrates))])
    
    binds=['bind2','bind3','bind5','bind6','bindconj5','bindconj6']
    buttons2=[]
    for x in binds:
        button=widgets.Button(layout=widgets.Layout(height='30px', min_width='120px'),description=x)
        buttons2.append(button)
        button.on_click(outputmol2)
    Ringbinders=widgets.HBox(children=[buttons2[x] for x in range(len(binds))])
    
    Cs=['C=O','C=NH','O','S','N']
    buttons3=[]
    for x in Cs:
        button=widgets.Button(layout=widgets.Layout(height='30px', min_width='120px'),description=x)
        buttons3.append(button)
        button.on_click(outputmol3)
    SmallChanges=widgets.HBox(children=[buttons3[x] for x in range(len(Cs))])
        
    stack = widgets.Stack([Substrates,Ringbinders,SmallChanges])
    
    choose=widgets.Dropdown(options=['Substrates','Ringbinders','SmallChanges'])
    widgets.jslink((choose,'index'),(stack,'selected_index'))
    
    chooseout=widgets.interactive_output(chooseshow, {'i': i, 'i1': i1, 'i2': i2, 'choose': choose})
    defaultout=widgets.interactive_output(default, {'strictele1': strictele1})
    
    box_layout = widgets.Layout(overflow='scroll hidden', border='2px solid black',
                    width='65%',height='550px',display='flex')
    box_layout2 = widgets.Layout(overflow='scroll hidden', border='2px solid black',
                    width='35%',height='550px',display='flex')
    output3.layout=box_layout2
    clean=widgets.Button(description='Clean Data')
    clean.on_click(clear)
    
    with output0:
        display(widgets.VBox([loc,widgets.HBox([file,openfile]),widgets.Label('OR',style=dict(font_weight='bold',text_decoration='underline')),\
                          widgets.HBox([smart,smartbutton])]))
    
    allchoose=widgets.HBox([widgets.Label('Choose your step: '),choose,clean])
    display(widgets.HBox([widgets.VBox([output0,output],layout=box_layout),output3]))