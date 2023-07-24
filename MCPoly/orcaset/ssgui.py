import ipywidgets as widgets
import sys
import os
import re
from ase.io import read
from ase import Atoms
from IPython import display

mydir = os.path.dirname( __file__ )
orcadir = os.path.join(mydir, '..', 'orcaset')
sys.path.append(orcadir)

from view3dchoose import view3dchoose
from XYZtoINP import XYZtoINP
from orca import orca
import time as t

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from status import status

mydir = os.path.dirname( __file__ )
viewdir = os.path.join(mydir, '..', 'view3d')
sys.path.append(viewdir)
from view3d import view3d

mydir = os.path.dirname( __file__ )
ssdir = os.path.join(mydir, '..', 'sscurve')
sys.path.append(ssdir)
from single import single

def ssgui():
    """
        The method to use GUI to set ORCA input files, calculate the file with external force by ORCA and draw stress-strain curve in handy!
        You will know how to do as long as you get into the GUI platform.

        TIPS: 1. This function can calculate ORCA input files, and draw single stress-strain curve and save data texts as well.
              2. Take caution when you set External Force site. Don't set number out of range.
              3. To see the status of your optimisation, you can set another window and use status.figuretraj .
              4. The visualisation will also risk the flexibility of the task, to make you programme more flexible, please directly use the code.
    """
    output=widgets.Output()
    output2=widgets.Output()
    addoutput1=widgets.Output()
    addoutput2=widgets.Output()
    statusoutput=widgets.Output()
    geooutput=widgets.Output()
    curveoutput=widgets.Output()
    fromoutput=widgets.Output()
    case=widgets.Output()

    def checkbroken(file,strain,forcestep,aim):
        if strain<=forcestep or strain<=1:
            return 0
        atoms=read(file+'_{0:.3f}.xyz'.format(strain))
        atoms0=read(file+'_{0:.3f}.xyz'.format(strain-forcestep))
        d=atoms.get_distance(aim[0],aim[1])
        d0=atoms0.get_distance(aim[0],aim[1])
        if d/d0>=1.4 or d/d0<=0.7:
            return 1
        else:
            return 0

    def showfrom(fromshow,file,loc,aim1,aim2):
        fromoutput.clear_output()
        opath=os.getcwd()
        polymers0=[]
        for path in os.listdir(loc):
            if os.path.isfile(os.path.join(loc, path)):
                a=re.search('{0}_'.format(file), path)
                if a:
                    b=re.search('_trj.xyz', path)
                    if b:
                        c=re.match('{0}_trj.xyz'.format(file), path)
                        if not c:
                            polymers0.append(path[:-8])
        polynums=[]
        polymers=[]
        polymers0.sort()
        for polymer in polymers0:
            a=re.split('_',polymer)
            polynums.append(eval(a[1]))
            polyname=a[0]
        polynums.sort()
        
        for num in polynums:
            polymers.append('{0}_{1:.3f}'.format(polyname,num))
        
        with fromoutput:
            for i,polymer in enumerate(polymers):
                print(polymer+': {0:.6f} Hartree'.format(status(polymer,loc).status(figureonly=True)[-1]))
                if fromshow==False:
                    view3d(polymer,loc,width=450,height=400)
                else:
                    view3dchoose(polymer,loc,choose=[aim1,aim2],width=450,height=400)
        os.chdir(opath)

    def curveshow(curve,file,loc):
        output2.clear_output()
        curveoutput.clear_output()
        with curveoutput:
            if curve==True:
                single(file,loc).autocurve(savefig=False,savedata=False)
                display(widgets.HBox([savechart,savedata]))

    def addshow1(addox1):
        addoutput1.clear_output()
        if addox1==True:
            addox2.disabled=False
            with addoutput1:
                display(strain2,addforce2)
        elif addox1==False:
            addox2.value=False
            addox2.disabled=True

    def addshow2(addox2):
        addoutput2.clear_output()
        if addox2==True:
            with addoutput2:
                display(strain3,addforce3)

    def geoshow(file,loc,aim1,aim2,method,bs):
        atomsx=read(loc+file+'.xyz')
        num=len(atomsx)
        output.clear_output()
        with output:
            try:
                if aim1==0 and aim2==0:
                    view3d(file,loc,width=450,height=400)
                else:
                    if aim1>=num or aim2>=num:
                        print('Setting Error. The maximum number of this system is {0}.\n We will automatically change the number into the maximum number.'.format(num-1))
                    if aim1>=num:
                        aim1=num-1
                    if aim2>=num:
                        aim2=num-1
                    aim=[aim1,aim2]
                    print('Atom1: {0}      Atom2: {1}'.format(aim1,aim2))
                    view3dchoose(file,loc,choose=aim,width=450,height=400)
            except:
                if file=='':
                    print('Please input the file to show the geometry structure on this screen.')
                else:
                    print('No such file.')
            if method!='' and bs!='':
                display(widgets.HBox([widgets.Label('ORCA Location:'),orcaloc]),startorca)
    
    def savechartx(button):
        single(file.value,loc.value).autocurve(savedata=False)
        
    def savedatax(button):
        single(file.value,loc.value).autocurve(savefig=False)
    
    def intro_status3d(button):
        atomsx=read(loc.value+file.value+'.xyz')
        num=len(atomsx)
        if aim1.value>=num:
            aim1.value=num-1
        if aim2.value>=num:
            aim2.value=num-1
        output.clear_output()
        with output2:
            display(widgets.Label('File: {0}.xyz'.format(file.value),style=dict(font_weight='bold')))
            if roomox.value==False and coreox.value==False:
                XYZtoINP(file.value,loc=loc.value,method=method.value,basis_set=bs.value,opt=True,freq=freq.value,\
                 external_force=True,aim=[aim1.value,aim2.value],strain=strain.value,maxcore=room.value*1024,corenum=corenum.value)
            elif roomox.value==True and coreox.value==True:
                XYZtoINP(file.value,loc=loc.value,method=method.value,basis_set=bs.value,opt=True,freq=freq.value,\
                 external_force=True,aim=[aim1.value,aim2.value],strain=strain.value)
            elif roomox.value==False and coreox.value==True:
                XYZtoINP(file.value,loc=loc.value,method=method.value,basis_set=bs.value,opt=True,freq=freq.value,\
                 external_force=True,aim=[aim1.value,aim2.value],strain=strain.value,maxcore=room.value*1024)
            elif roomox.value==True and coreox.value==False:
                XYZtoINP(file.value,loc=loc.value,method=method.value,basis_set=bs.value,opt=True,freq=freq.value,\
                 external_force=True,aim=[aim1.value,aim2.value],strain=strain.value,corenum=corenum.value)
        os.system('cp {0}.inp {0}_0.000.inp'.format(file.value))
        keys=[]
        if addforce.value==0:
            raise ValueError("You can't set the force add per step zero!")
        forcestep={strain.value:addforce.value}
        if addox1.value==True:
            if addforce2.value==0:
                raise ValueError("You can't set the force add per step zero!")
            forcestep={strain.value: addforce.value, strain2.value: addforce2.value}
        elif addox1.value==True and addox2.value==True:
            if addforce3.value==0:
                raise ValueError("You can't set the force add per step zero!")
            forcestep={strain.value: addforce.value, strain2.value: addforce2.value, strain3.value: addforce3.value}
        for key in forcestep:
            keys.append(key)
        keys.reverse()
        key=keys[0]
        ii=-1
        while 1:
            ii=ii+1
            with case:
                display(widgets.HBox([widgets.Label(file.value+'_{0:.3f}'.format(strain.value)), widgets.Label('  Processing...'), widgets.Label('[{0}]'.format(t.ctime(t.time())))]))
            orca(file.value+'_{0:.3f}'.format(strain.value),orcaloc.value,loc.value,loc.value)
            s=status(file.value+'_{0:.3f}'.format(strain.value),loc.value).status(statusonly=True)
            
            case.clear_output()
            if s==4:
                xoy=checkbroken(file.value,strain.value,forcestep[key],[aim1.value,aim2.value])
                if xoy==1:
                    with output2:
                        display(widgets.HBox([widgets.Valid(description=file.value+'_{0:.3f}'.format(strain.value),value=True),\
                                      widgets.Label('Broken.'), widgets.Label('[{0}]'.format(t.ctime(t.time())))]))
                    fromoutput.clear_output()
                    with fromoutput:
                        for i,poly in enumerate(polys):
                            print(poly+': {0:.6f} Hartree'.format(status(poly,loc.value).status(figureonly=True)[-1]))
                            if fromshow==False:
                                view3d(poly,loc.value,width=450,height=400)
                            else:
                                view3dchoose(poly,loc.value,choose=[aim1.value,aim2.value],width=450,height=400)
                    break
                else:
                    case.clear_output()
                    with case:
                        display(widgets.Valid(description=file.value+'_{0:.3f}'.format(strain.value),value=False,\
                                readout='Retrying...'), widgets.Label('[{0}]'.format(t.ctime(t.time()))))
                    orca(file.value+'_{0:.3f}'.format(strain.value),orcaloc.value,loc.value,loc.value)
                    s=status(file.value,loc.value).status(statusonly=True)
                    if s==4:
                        case.clear_output()
                        with case:
                            display(widgets.Valid(description=file.value+'_{0:.3f}'.format(strain.value),value=False,\
                                readout='Error!'), widgets.Label('[{0}]'.format(t.ctime(t.time()))))
            elif s==2:
                with output2:
                    display(widgets.HBox([widgets.Valid(description=file.value+'_{0:.3f}'.format(strain.value),value=True),\
                                  widgets.Label('Completed.'), widgets.Label('[{0}]'.format(t.ctime(t.time())))]))
            if ii==0:
                output.clear_output()
                with output:
                    display(fromoutput)
            
            opath=os.getcwd()
            polys=[]
            try:
                for path in os.listdir(loc.value):
                    if os.path.isfile(os.path.join(loc.value, path)):
                        a=re.search('{0}_'.format(file.value), path)
                        if a:
                            b=re.search('_trj.xyz', path)
                            if b:
                                c=re.search('{0}_trj.xyz'.format(file.value), path)
                                if not c:
                                    polys.append(path[:-8])
                polys.sort()
                fromoutput.clear_output()
                with fromoutput:
                    for i,poly in enumerate(polys):
                        print(poly+': {0:.6f} Hartree'.format(status(poly,loc.value).status(figureonly=True)[-1]))
                        if fromshow==False:
                            view3d(poly,loc.value,width=450,height=400)
                        else:
                            view3dchoose(poly,loc.value,choose=[aim1.value,aim2.value],width=450,height=400)
            except:
                pass
            os.chdir(opath)
            
            for key in keys:
                if strain.value>=key:
                    XYZtoINP(file.value+'_{0:.3f}'.format(strain.value),inpname=file.value+'_{0:.3f}'.format(strain.value+forcestep[key]),\
                         method=method.value,basis_set=bs.value,opt=True,freq=freq.value,\
                         external_force=True,aim=[aim1.value,aim2.value],strain=strain.value+forcestep[key],maxcore=room.value*1024,\
                         corenum=corenum.value)
                    strain.value=strain.value+forcestep[key]
                    break
        with output:
            display(fromshow)
        with output2:
            display(curve)

    def statusresult(button):
        statusoutput.clear_output()
        with statusoutput:
            while 1:
                if s==1:
                    break
                t.sleep(30)
            
    options=[]
    loc=widgets.Text(description='Location:',value='./')
    for path in os.listdir(loc.value):
        if os.path.isfile(os.path.join(loc.value, path)):
            a=re.search('.xyz', path)
            if a:
                options.append(path[:-4])
    file=widgets.Combobox(placeholder='Input the file name',options=options,
         description='xyz File:')
    ssfiles=[]
    for path in os.listdir(loc.value):
        if os.path.isfile(os.path.join(loc.value, path)):
            a=re.search('.xyz', path)
            if a:
                b=re.search(file.value,path)
                if b:
                    ssfiles.append(path[:-4])
    method=widgets.Text(description='Method:')
    bs=widgets.Text(description='Basis Set:')
    freq=widgets.ToggleButton(value=False,description='Frequency',layout=widgets.Layout(width='33.3%'))
    room=widgets.IntText(description='Space(GB):',value=4,min=1,layout=widgets.Layout(width='50%'))
    roomox=widgets.ToggleButton(value=False,description='Default Space')
    corenum=widgets.IntText(description='Core:',value=4,min=1,layout=widgets.Layout(width='50%'))
    coreox=widgets.ToggleButton(value=False,description='Default CoreNumber')
    aim1=widgets.BoundedIntText(description='Atom 1:',value=0,min=0,layout=widgets.Layout(width='50%'))
    aim2=widgets.BoundedIntText(description='Atom 2:',value=0,min=0,layout=widgets.Layout(width='47.5%'))
    strain=widgets.BoundedFloatText(description='Start Force(nN):',value=0.000,min=0.000,step=0.100,style = {'description_width': 'initial'})
    addforce=widgets.BoundedFloatText(description='Force Add per Step(nN):',value=0.100,min=0.000,step=0.100,style = {'description_width': 'initial'})
    addox1=widgets.ToggleButton(value=False,description='The Second Stage',style=dict(font_weight='bold'))
    strain2=widgets.BoundedFloatText(description='Till Force(nN):',value=0.000,min=0.000,step=0.100,style = {'description_width': 'initial'})
    addforce2=widgets.BoundedFloatText(description='Force Add per Step(nN):',value=0.100,min=0.000,step=0.100,style = {'description_width': 'initial'})
    addox2=widgets.ToggleButton(value=False,description='The Third Stage',style=dict(font_weight='bold'))
    strain3=widgets.BoundedFloatText(description='Till Force(nN):',value=0.000,min=0.000,step=0.100,style = {'description_width': 'initial'})
    addforce3=widgets.BoundedFloatText(description='Force Add per Step(nN):',value=0.100,min=0.000,step=0.100,style = {'description_width': 'initial'})
    startorca=widgets.Button(description='Start ORCA Calculation',layout=widgets.Layout(width='99%'),style=dict(font_weight='bold',text_decoration='underline'))
    orcaloc=widgets.Text(value='./')
    curve=widgets.ToggleButton(description='Stress-Strain Curve',layout=widgets.Layout(width='99%'),style=dict(font_weight='bold',text_decoration='underline'))
    fromshow=widgets.ToggleButton(description='Show the site of Force',value=True,layout=widgets.Layout(width='99%'),style=dict(font_weight='bold',text_decoration='underline'))
    savechart=widgets.Button(description='Save chart ( .png)',layout=widgets.Layout(width='50%'))
    savedata=widgets.Button(description='Save data (.txt)',layout=widgets.Layout(width='47.5%'))
    
    startorca.on_click(intro_status3d)
    savechart.on_click(savechartx)
    savedata.on_click(savedatax)
    
    title1=widgets.Label('File Details',style=dict(font_weight='bold'))
    title2=widgets.Label('Calculation Methods',style=dict(font_weight='bold'))
    title3=widgets.Label('Constant External Force Settings',style=dict(font_weight='bold'))
    
    addout1=widgets.interactive_output(addshow1, {'addox1': addox1})
    addout2=widgets.interactive_output(addshow2, {'addox2': addox2})
    geoout=widgets.interactive_output(geoshow, {'file': file, 'loc': loc, 'aim1': aim1, 'aim2': aim2, 'method': method, 'bs': bs})
    ssout=widgets.interactive_output(curveshow, {'curve': curve, 'file': file, 'loc': loc})
    fromout=widgets.interactive_output(showfrom, {'fromshow': fromshow, 'file': file, 'loc': loc, 'aim1': aim1, 'aim2': aim2})
    
    widgets.jslink((roomox,'value'),(room,'disabled'))
    widgets.jslink((coreox,'value'),(corenum,'disabled'))
    
    box_layout1 = widgets.Layout(display='flex', flew_flow='column', overflow='scroll hidden',\
                    align_items='stretch', border='solid', height='550px', width='31%')
    box_layout2 = widgets.Layout(display='flex', flew_flow='column', overflow='scroll hidden',\
                    align_items='stretch', border='solid', height='550px', width='33%')
    
    frame1=widgets.VBox([title1,loc,file,widgets.HBox([title2,freq]),method,bs,widgets.HBox([corenum,coreox]),\
                     widgets.HBox([room,roomox]),title3,widgets.HBox([aim1,aim2]),strain,addforce,addox1,addoutput1,addox2,addoutput2],layout=box_layout1)
    frame3=widgets.VBox([output2,case,curveoutput])
    output.layout=box_layout2
    frame3.layout=box_layout2
    curveoutput.layout=widgets.Layout(display='flex', flew_flow='column', overflow='scroll hidden',\
                    align_items='stretch', height='550px', width='100%')
    output2.layout=widgets.Layout(display='flex', flew_flow='column', overflow='scroll hidden',\
                    align_items='stretch', width='100%')
    display(widgets.HBox([frame1,output,frame3]))
