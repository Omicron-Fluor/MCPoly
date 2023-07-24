import ipywidgets as widgets
from .XYZtoINP import XYZtoINP
from .orca import orca
from .multiorca import multiorca
from .ssorca import ssorca
from ase.io import read
from ase import Atoms
import sys
import os
import re
from .view3dchoose import view3dchoose
from IPython import display

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from status import status

mydir = os.path.dirname( __file__ )
viewdir = os.path.join(mydir, '..', 'view3d')
sys.path.append(viewdir)
from view3d import view3d

def gui():
    """
        The method to use GUI to set ORCA input files and calculate the file by ORCA get in handy!
        You will know how to do as long as you get into the GUI platform.

        TIPS: 1. Try to directly input the figure of external force instead of scrolling it, or you will see a lot of input files.
              2. If you specialized in calcuating mechanical properties, orcaset.ssgui is a better option.
              3. Take caution when you set Scanning or External Force site. Don't set number out of range.
              4. To see the status of your optimisation, you can set another window and use status.figuretraj .
              5. The visualisation will also risk the flexibility of the task, to make you programme more flexible, please directly use the code.
    """
    output=widgets.Output()
    output2=widgets.Output()
    scanoutput=widgets.Output()
    forceoutput=widgets.Output()
    eoutput=widgets.Output()
    geooutput=widgets.Output()
        
    def scanshow(scanox):
        aim1.value=0
        aim2.value=0
        scanoutput.clear_output()
        opt.disabled=False
        ts.disabled=False
        freq.disabled=False
        forceox.disabled=False
        if scanox==True:
            opt.disabled=True
            ts.disabled=True
            freq.disabled=True
            forceox.value=False
            forceox.disabled=True
            with scanoutput:
                display(widgets.VBox([widgets.HBox([aim1,aim2]),stretch,scanstep]))

    def forceshow(forceox):
        aim1.value=0
        aim2.value=0
        forceoutput.clear_output()
        ts.disabled=False
        scanox.disabled=False
        if forceox==True:
            ts.disabled=True
            scanox.value=False
            scanox.disabled=True
            opt.value=True
            with forceoutput:
                display(widgets.VBox([widgets.HBox([aim1,aim2]),strain]))

    def eshow(eox):
        eoutput.clear_output()
        if eox==True:
            with eoutput:
                if scanox.value==True:
                    status(file.value+'_scan',loc.value).status()
                elif forceox.value==True:
                    status(file.value+'_{0:.3f}'.format(strain.value),loc.value).status()
                else:
                    status(file.value,loc.value).status()

    def geoshow(geoox):
        geooutput.clear_output()
        if geoox==True:
            with geooutput:
                if scanox.value==True:
                    status(file.value+'_scan',loc.value).figuretraj()
                elif forceox.value==True:
                    status(file.value+'_{0:.3f}'.format(strain.value),loc.value).figuretraj()
                else:
                    status(file.value,loc.value).figuretraj()
    
    def guiconvert(file,loc,method,bs,opt,freq,ts,scanox,forceox,aim1,aim2,scanstep,stretch,strain,room,roomox,corenum,coreox,e,state):
        if scanox==True or forceox==True:
            atomsx=read(loc+file+'.xyz')
            num=len(atomsx)
            if aim1>=num:
                aim1=num-1
            if aim2>=num:
                aim2=num-1
        try:
            if roomox==True and coreox==True:
                if scanox==True:
                    XYZtoINP(file,inpname=file+'_scan',loc=loc,method=method,basis_set=bs,opt=opt,scan=True,aim=[aim1,aim2],stretch=stretch,scanstep=scanstep,electron=e,state=state)
                elif forceox==True:
                    XYZtoINP(file,inpname='{0}_{1:.3f}'.format(file,strain),loc=loc,method=method,basis_set=bs,opt=opt,freq=freq,external_force=True,aim=[aim1,aim2],strain=strain,electron=e,state=state)
                else:
                    XYZtoINP(file,loc=loc,method=method,basis_set=bs,opt=opt,freq=freq,ts=ts,electron=e,state=state)
            elif roomox==False and coreox==False:
                if scanox==True:
                    XYZtoINP(file,inpname=file+'_scan',loc=loc,method=method,basis_set=bs,opt=opt,scan=True,aim=[aim1,aim2],stretch=stretch,scanstep=scanstep,maxcore=room*1024,corenum=corenum,electron=e,state=state)
                elif forceox==True:
                    XYZtoINP(file,inpname='{0}_{1:.3f}'.format(file,strain),loc=loc,method=method,basis_set=bs,opt=opt,freq=freq,external_force=True,aim=[aim1,aim2],strain=strain,maxcore=room*1024,corenum=corenum,electron=e,state=state)
                else:
                    XYZtoINP(file,loc=loc,method=method,basis_set=bs,opt=opt,freq=freq,ts=ts,maxcore=room*1024,corenum=corenum,electron=e,state=state)
            elif roomox==True and coreox==False:
                if scanox==True:
                    XYZtoINP(file,inpname=file+'_scan',loc=loc,method=method,basis_set=bs,opt=opt,scan=True,aim=[aim1,aim2],stretch=stretch,scanstep=scanstep,corenum=corenum,electron=e,state=state)
                elif forceox==True:
                    XYZtoINP(file,inpname='{0}_{1:.3f}'.format(file,strain),loc=loc,method=method,basis_set=bs,opt=opt,freq=freq,external_force=True,aim=[aim1,aim2],strain=strain,corenum=corenum,electron=e,state=state)
                else:
                    XYZtoINP(file,inpname=file+'_scan',loc=loc,method=method,basis_set=bs,opt=opt,freq=freq,ts=ts,corenum=corenum,electron=e,state=state)
            elif roomox==False and coreox==True:
                if scanox==True:
                    XYZtoINP(file,inpname=file+'_scan',loc=loc,method=method,basis_set=bs,opt=opt,scan=True,aim=[aim1,aim2],stretch=stretch,scanstep=scanstep,maxcore=room*1024,electron=e,state=state)
                elif forceox==True:
                    XYZtoINP(file,inpname='{0}_{1:.3f}'.format(file,strain),loc=loc,method=method,basis_set=bs,opt=opt,freq=freq,external_force=True,aim=[aim1,aim2],strain=strain,maxcore=room*1024,electron=e,state=state)
                else:
                    XYZtoINP(file,loc=loc,method=method,basis_set=bs,opt=opt,freq=freq,ts=ts,maxcore=room*1024,electron=e,state=state)
        
        except:
            print('Setting Error.')
            return None
        print()
        if scanox==True:
            print(file+'_scan'+'.inp\n')
            f=open(loc+file+'_scan'+'.inp','r')
        elif forceox==True:
            print('{0}_{1:.3f}.inp\n'.format(file,strain))
            f=open(loc+'{0}_{1:.3f}.inp'.format(file,strain),'r')
        else:
            print(file+'.inp\n')
            f=open(loc+file+'.inp','r')
        for line in f:
            print(line[:-1])
        f.close()
        if forceox==True:
            os.system('rm {2}/{0}_{1:.3f}.inp'.format(file,strain,loc))

    def screen3d(file,loc,aim1,aim2,scanox,forceox):
        aim=[]
        try:
            if scanox==True or forceox==True:
                atomsx=read(loc+file+'.xyz')
                num=len(atomsx)
                if aim1>=num or aim2>=num:
                    print('Setting Error. The maximum number of this system is {0}.\n We will automatically change the number into the maximum number.'.format(num-1))
                if aim1>=num:
                    aim1=num-1
                if aim2>=num:
                    aim2=num-1
                aim=[aim1,aim2]
                print('Atom1: {0}      Atom2: {1}'.format(aim1,aim2))
                view3dchoose(file,loc,choose=aim,width=400,height=400)
            else:
                view3d(file,loc,width=400,height=500)
        except:
            print('Setting Error.')

    def createinp(button):
        output.clear_output()
        output2.clear_output()
        with output:
            out=widgets.interactive_output(guiconvert, {'file': file, 'loc': loc, 'method': method, 'bs': bs, 'opt': opt, 'freq': freq,\
                                                    'ts': ts, 'scanox': scanox, 'forceox': forceox, 'aim1': aim1, 'aim2': aim2,\
                                                    'stretch': stretch, 'scanstep': scanstep, 'strain': strain, 'room': room,\
                                                    'roomox': roomox, 'corenum': corenum, 'coreox': coreox, 'e': e, 'state': state})
            startorca=widgets.Button(description='Start ORCA Calculation',layout=widgets.Layout(width='99%'),style=dict(font_weight='bold',text_decoration='underline'))
            startorca.on_click(intro_status3d)
            display(widgets.VBox([out,widgets.HBox([widgets.Label('ORCA Location:'),orcaloc]),startorca]))
        with output2:
            out2=widgets.interactive_output(screen3d, {'file': file, 'loc': loc, 'aim1': aim1, 'aim2': aim2, 'scanox': scanox, 'forceox': forceox}) 
            display(out2)

    def intro_status3d(button):
        if scanox.value==True or forceox.value==True:
            atomsx=read(loc.value+file.value+'.xyz')
            num=len(atomsx)
            if aim1.value>=num:
                aim1.value=num-1
            if aim2.value>=num:
                aim2.value=num-1
        with output:
            button.description='Processing...'
        if scanox.value==True:
            orca(file.value+'_scan',orcaloc=orcaloc.value)
        elif forceox.value==True:
            if roomox.value==True and coreox.value==True:
                XYZtoINP(file,inpname='{0}_{1:.3f}'.format(file.value,strain.value),loc=loc.value,method=method.value,basis_set=bs.value,opt=opt.value,freq=freq.value,external_force=True,aim=[aim1.value,aim2.value],strain=strain.value,electron=e.value,state=state.value)
            elif roomox.value==False and coreox.value==False:
                XYZtoINP(file.value,inpname='{0}_{1:.3f}'.format(file.value,strain.value),loc=loc.value,method=method.value,basis_set=bs.value,opt=opt.value,freq=freq.value,external_force=True,aim=[aim1.value,aim2.value],strain=strain.value,maxcore=room.value*1024,corenum=corenum.value,electron=e.value,state=state.value)
            else:
                XYZtoINP(file.value,inpname='{0}_{1:.3f}'.format(file.value,strain.value),loc=loc.value,method=method.value,basis_set=bs.value,opt=opt.value,freq=freq.value,external_force=True,aim=[aim1.value,aim2.value],strain=strain.value,maxcore=room.value*1024,electron=e.value,state=state.value)
            orca('{0}_{1:.3f}'.format(file.value,strain.value),orcaloc=orcaloc.value)
        else:
            orca(file.value,orcaloc=orcaloc.value)
        with output2:
            if scanox.value==True:
                if status(file.value+'_scan',loc.value).status(figureonly=True,statusonly=True)==2:
                    display(widgets.HBox([widgets.Valid(description=file.value+'_scan',value=True),widgets.Label('  Calculation completed!')]))
                elif status(file.value+'_scan',loc.value).status(figureonly=True,statusonly=True)==4:
                    display(widgets.Valid(description=file.value+'_scan',value=False,readout='Calculation aborted!'))
            elif forceox.value==True:
                if status('{0}_{1:.3f}'.format(file.value,strain.value),loc.value).status(figureonly=True,statusonly=True)==2:
                    display(widgets.HBox([widgets.Valid(description='{0}_{1:.3f}'.format(file.value,strain.value),value=True),widgets.Label('  Calculation completed!')]))
                elif status('{0}_{1:.3f}'.format(file.value,strain.value),loc.value).status(figureonly=True,statusonly=True)==4:
                    display(widgets.Valid(description='{0}_{1:.3f}'.format(file.value,strain.value),value=False,readout='Calculation aborted!'))
            else:
                if status(file.value,loc.value).status(figureonly=True,statusonly=True)==2:
                    display(widgets.HBox([widgets.Valid(description=file.value,value=True),widgets.Label('  Calculation completed!')]))
                elif status(file.value,loc.value).status(figureonly=True,statusonly=True)==4:
                    display(widgets.Valid(description=file.value,value=False,readout='Calculation aborted!'))
            display(introstatus)
        output.clear_output()
        
    def status3d(button):
        output2.clear_output()
        with output2:
            display(widgets.VBox([eox,eoutput,geoox,geooutput]))

    options=[]
    loc=widgets.Text(description='Location:',value='./')
    for path in os.listdir(loc.value):
        if os.path.isfile(os.path.join(loc.value, path)):
            a=re.search('.xyz', path)
            if a:
                options.append(path[:-4])
    file=widgets.Combobox(placeholder='Input the file name',options=options,
         description='xyz File:')
    method=widgets.Text(description='Method:')
    bs=widgets.Text(description='Basis Set:')
    opt=widgets.ToggleButton(value=False,description='Optimisation',layout=widgets.Layout(width='33.3%'))
    ts=widgets.ToggleButton(value=False,description='TS',layout=widgets.Layout(width='33.3%'))
    freq=widgets.ToggleButton(value=False,description='Frequency',layout=widgets.Layout(width='33.3%'))
    room=widgets.IntText(description='Space(GB):',value=4,min=1,layout=widgets.Layout(width='50%'))
    roomox=widgets.ToggleButton(value=False,description='Default Space')
    corenum=widgets.IntText(description='Core:',value=4,min=1,layout=widgets.Layout(width='50%'))
    coreox=widgets.ToggleButton(value=False,description='Default CoreNumber')
    e=widgets.IntText(description='Charge:',value=0,layout=widgets.Layout(width='50%'))
    state=widgets.IntText(description='State:',value=1,layout=widgets.Layout(width='47.5%'))
    scanox=widgets.ToggleButton(value=False,description='Scan Settings',style=dict(font_weight='bold'),layout=widgets.Layout(width='98.5%'))
    forceox=widgets.ToggleButton(value=False,description='Constant External Force Settings',style=dict(font_weight='bold'),layout=widgets.Layout(width='98.5%'))
    aim1=widgets.BoundedIntText(description='Atom 1:',value=0,min=0,layout=widgets.Layout(width='50%'))
    aim2=widgets.BoundedIntText(description='Atom 2:',value=0,min=0,layout=widgets.Layout(width='47.5%'))
    strain=widgets.BoundedFloatText(description='Force(nN):',value=0.000,min=0.000,step=0.100)
    stretch=widgets.FloatText(description='Distance(Å):',value=1.000,step=0.100)
    scanstep=widgets.BoundedIntText(description='Steps:',value=10,min=0)
    scanout=widgets.interactive_output(scanshow, {'scanox': scanox})
    forceout=widgets.interactive_output(forceshow, {'forceox': forceox})
    eox=widgets.ToggleButton(value=False,description='Energy Status',style=dict(font_weight='bold'),layout=widgets.Layout(width='98.5%'))
    geoox=widgets.ToggleButton(value=False,description='Structure Status',style=dict(font_weight='bold'),layout=widgets.Layout(width='98.5%'))
    
    eout=widgets.interactive_output(eshow, {'eox': eox})
    geoout=widgets.interactive_output(geoshow, {'geoox': geoox})
    orcaloc=widgets.Text(value='./')
    
    title1=widgets.Label('File Details',style=dict(font_weight='bold'))
    title2=widgets.Label('Calculation Methods',style=dict(font_weight='bold'))
    
    box_layout1 = widgets.Layout(display='flex', flew_flow='column', overflow='scroll hidden',\
                    align_items='stretch', border='solid', height='550px', width='31%')
    box_layout2 = widgets.Layout(display='flex', flew_flow='column', overflow='scroll hidden',\
                    align_items='stretch', border='solid', height='550px', width='33%')
    widgets.jslink((roomox,'value'),(room,'disabled'))
    widgets.jslink((coreox,'value'),(corenum,'disabled'))
    create=widgets.Button(description='Preview ORCA Input File',style=dict(font_weight='bold',text_decoration='underline'),
                      layout=widgets.Layout(width='98.5%'))
    introstatus=widgets.Button(description='See Current Calculation Status',style=dict(font_weight='bold',text_decoration='underline'),
                      layout=widgets.Layout(width='98.5%'))

    frame1=widgets.VBox([title1,loc,file,title2,widgets.HBox([opt,ts,freq]),method,bs,widgets.HBox([corenum,coreox]),\
                     widgets.HBox([room,roomox]),widgets.HBox([e,state]),scanox,scanoutput,forceox,forceoutput,create],layout=box_layout1)
    output.layout=box_layout2
    output2.layout=box_layout2
    create.on_click(createinp)
    introstatus.on_click(status3d)

    display(widgets.HBox([frame1,output,output2]))
