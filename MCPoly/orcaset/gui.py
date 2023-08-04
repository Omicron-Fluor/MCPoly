import sys
import os
import re
import ipywidgets as widgets
from .XyzToInp import XyzToInp
from .orca import orca
from .multiorca import multiorca
from .ssorca import ssorca
from ase.io import read
from ase import Atoms
from .view3dchoose import view3dchoose
from IPython.display import display

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
    output = widgets.Output()
    output2 = widgets.Output()
    scanoutput = widgets.Output()
    forceoutput = widgets.Output()
    eoutput = widgets.Output()
    geooutput = widgets.Output()
    mdoutput = widgets.Output()
    global geore
        
    def scanshow(scanox):
        aim1.value = 0
        aim2.value = 0
        scanoutput.clear_output()
        opt.disabled = False
        ts.disabled = False
        freq.disabled = False
        forceox.disabled = False
        if scanox == True:
            opt.disabled = True
            ts.disabled = True
            freq.disabled = True
            forceox.value = False
            forceox.disabled = True
            with scanoutput:
                display(widgets.VBox([widgets.HBox([aim1, aim2]), stretch, scanstep]))

    def forceshow(forceox):
        aim1.value = 0
        aim2.value = 0
        forceoutput.clear_output()
        ts.disabled = False
        scanox.disabled = False
        if forceox == True:
            ts.disabled = True
            scanox.value = False
            scanox.disabled = True
            opt.value = True
            with forceoutput:
                display(widgets.VBox([widgets.HBox([aim1, aim2]), strain]))

    def eshow(eox):
        eoutput.clear_output()
        if eox == True:
            with eoutput:
                if scanox.value == True:
                    status(file.value+'_scan', loc.value).status()
                elif forceox.value == True:
                    status(file.value+'_{0:.3f}'.format(strain.value), loc.value).status()
                else:
                    status(file.value, loc.value).status()

    def geoshow(geoox):
        geooutput.clear_output()
        if geoox == True:
            with geooutput:
                if scanox.value == True:
                    geore = status(file.value+'_scan', loc.value)
                    savetotal = geore.steps()
                    savenumchoose.max = savetotal
                    def save3Dshow(savenumchoose):
                        savenum = savenumchoose
                        geore.figure(num=savenumchoose)
                    savenumchooseout = widgets.interactive_output(save3Dshow,\
                                                                  {'savenumchoose': savenumchoose})
                    display(savenumchoose)
                    display(savenumchooseout)
                elif forceox.value == True:
                    geore = status(file.value+'_{0:.3f}'.format(strain.value),\
                                   loc.value)
                    savetotal = geore.steps()
                    savenumchoose.max = savetotal
                    def save3Dshow(savenumchoose):
                        savenum = savenumchoose
                        geore.figure(num=savenumchoose)
                    savenumchooseout = widgets.interactive_output(save3Dshow,\
                                                                  {'savenumchoose': savenumchoose})
                    display(savenumchoose)
                    display(savenumchooseout)
                elif mdox.value == True:
                    geore = status(file.value+'_{0:.3f}'.format(strain.value),\
                                   loc.value)
                    savetotal = geore.steps(MD=True)
                    savenumchoose.max = savetotal
                    def save3Dshow(savenumchoose):
                        savenum = savenumchoose
                        geore.figure(num=savenumchoose)
                    savenumchooseout = widgets.interactive_output(save3Dshow,\
                                                                  {'savenumchoose': savenumchoose})
                    display(savenumchoose)
                    display(savenumchooseout)
                else:
                    geore = status(file.value, loc.value)
                    savetotal = geore.steps()
                    savenumchoose.max = savetotal
                    def save3Dshow(savenumchoose):
                        savenum = savenumchoose
                        geore.figure(num=savenumchoose)
                    savenumchooseout = widgets.interactive_output(save3Dshow,\
                                                                  {'savenumchoose': savenumchoose})
                    display(savenumchoose)
                    display(savenumchooseout)
    
    def mdshow(mdox):
        mdoutput.clear_output()
        if mdox == True:
            ts.disabled = True
            scanox.disabled = True
            opt.disabled = True
            freq.disabled = True
            forceox.disabled = True
            with mdoutput:
                display(widgets.VBox([widgets.HBox([tstep, initvel]),\
                                      widgets.HBox([tsNHC, timecon]), dumpcon, runstep]))
        else:
            ts.disabled = False
            scanox.disabled = False
            opt.disabled = False
            freq.disabled = False
            forceox.disabled = False
    
    def guiconvert(file, loc, method, bs, opt, freq, ts, scanox,\
                   forceox, aim1, aim2, scanstep, stretch, strain,\
                   room, roomox, corenum, coreox, e, state, maxiter,\
                   mdox, tstep, initvel, tsNHC, timecon, dumpcon, runstep):
        if scanox == True or forceox == True:
            atomsx = read(loc+file+'.xyz')
            num = len(atomsx)
            if aim1 >= num:
                aim1 = num - 1
            if aim2 >= num:
                aim2 = num - 1
        try:
            if roomox == True and coreox == True:
                if scanox == True:
                    XyzToInp(file, inpname=file+'_scan', fileloc=loc, method=method,\
                             basis_set=bs, opt=opt, scan=True, aim = [aim1, aim2],\
                             stretch=stretch, scanstep=scanstep, maxiter=maxiter,\
                             electron=e, state=state, MD=mdox, md_timestep=tstep,\
                             initvel=initvel, tsNHC=tsNHC, timecon=timecon,\
                             dumpcon=dumpcon, runstep=runstep)
                elif forceox == True:
                    XyzToInp(file, inpname = '{0}_{1:.3f}'.format(file, strain),\
                             fileloc=loc, method=method, basis_set=bs, opt=opt, freq=freq,\
                             external_force=True, aim = [aim1, aim2], strain=strain,\
                             maxiter=maxiter, electron=e, state=state, MD=mdox, md_timestep=tstep,\
                             initvel=initvel, tsNHC=tsNHC, timecon=timecon,\
                             dumpcon=dumpcon, runstep=runstep)
                else:
                    XyzToInp(file, fileloc=loc, method=method, basis_set=bs, opt=opt,\
                             freq=freq, ts = ts, electron=e, state=state, MD=mdox,\
                             md_timestep=tstep, initvel=initvel, tsNHC=tsNHC, timecon=timecon,\
                             dumpcon=dumpcon, runstep=runstep)
            elif roomox == False and coreox == False:
                if scanox == True:
                    XyzToInp(file, inpname = file+'_scan', fileloc=loc, method=method,\
                             basis_set=bs, opt=opt, scan=True, aim = [aim1, aim2],\
                             stretch=stretch, scanstep=scanstep, maxiter=maxiter,\
                             maxcore = room*1024, corenum=corenum, electron=e, state=state,\
                             MD=mdox, md_timestep=tstep, initvel=initvel, tsNHC=tsNHC,\
                             timecon=timecon, dumpcon=dumpcon, runstep=runstep)
                elif forceox == True:
                    XyzToInp(file, inpname = '{0}_{1:.3f}'.format(file, strain), fileloc=loc,\
                             method=method, basis_set=bs, opt=opt, freq=freq, external_force=True,\
                             aim = [aim1, aim2], strain=strain, maxiter=maxiter,\
                             maxcore = room*1024, corenum=corenum, electron=e, state=state,\
                             MD=mdox, md_timestep=tstep, initvel=initvel, tsNHC=tsNHC,\
                             timecon=timecon, dumpcon=dumpcon, runstep=runstep)
                else:
                    XyzToInp(file, fileloc=loc, method=method, basis_set=bs, opt=opt, freq=freq,\
                             ts=ts, maxiter=maxiter, maxcore = room*1024, corenum=corenum,\
                             electron=e, state=state, MD=mdox, md_timestep=tstep, initvel=initvel,\
                             tsNHC=tsNHC, timecon=timecon, dumpcon=dumpcon, runstep=runstep)
            elif roomox == True and coreox == False:
                if scanox == True:
                    XyzToInp(file, inpname = file+'_scan', fileloc=loc, method=method,\
                             basis_set=bs, opt=opt, scan=True, aim = [aim1, aim2],\
                             stretch=stretch, scanstep=scanstep, maxiter=maxiter, corenum=corenum,\
                             electron=e, state=state, MD=mdox, md_timestep=tstep, initvel=initvel,\
                             tsNHC=tsNHC, timecon=timecon, dumpcon=dumpcon, runstep=runstep)
                elif forceox == True:
                    XyzToInp(file, inpname = '{0}_{1:.3f}'.format(file, strain), fileloc=loc,\
                             method=method, basis_set=bs, opt=opt, freq=freq, external_force=True,\
                             aim = [aim1, aim2], strain=strain, maxiter=maxiter, corenum=corenum,\
                             electron=e, state=state, MD=mdox, md_timestep=tstep, initvel=initvel,\
                             tsNHC=tsNHC, timecon=timecon, dumpcon=dumpcon, runstep=runstep)
                else:
                    XyzToInp(file, inpname = file+'_scan', fileloc=loc, method=method,\
                             basis_set=bs, opt=opt, freq=freq, ts=ts, maxiter=maxiter,\
                             corenum=corenum, electron=e, state=state, MD=mdox, md_timestep=tstep,\
                             initvel=initvel, tsNHC=tsNHC, timecon=timecon,\
                             dumpcon=dumpcon, runstep=runstep)
            elif roomox == False and coreox == True:
                if scanox == True:
                    XyzToInp(file, inpname = file+'_scan', fileloc=loc, method=method,\
                             basis_set=bs, opt=opt, scan=True, aim = [aim1, aim2],\
                             stretch=stretch, scanstep=scanstep, maxiter=maxiter,\
                             maxcore = room*1024, electron=e, state=state, MD=mdox,\
                             md_timestep=tstep, initvel=initvel, tsNHC=tsNHC, timecon=timecon,\
                             dumpcon=dumpcon, runstep=runstep)
                elif forceox == True:
                    XyzToInp(file, inpname = '{0}_{1:.3f}'.format(file, strain), fileloc=loc,\
                             method=method, basis_set=bs, opt=opt, freq=freq, external_force=True,\
                             aim = [aim1, aim2], strain=strain, maxiter=maxiter,\
                             maxcore = room*1024, electron=e, state=state, MD=mdox,\
                             md_timestep=tstep, initvel=initvel, tsNHC=tsNHC,\
                             timecon=timecon, dumpcon=dumpcon, runstep=runstep)
                else:
                    XyzToInp(file, fileloc=loc, method=method, basis_set=bs, opt=opt,\
                             freq=freq, ts=ts, maxiter=maxiter, maxcore = room*1024,\
                             electron=e, state=state, MD=mdox, md_timestep=tstep, initvel=initvel,\
                             tsNHC=tsNHC, timecon=timecon, dumpcon=dumpcon, runstep=runstep)

        except:
            print('Setting Error.')
            return None
        if scanox == True:
            print(file+'_scan'+'.inp\n')
            f = open(loc+file+'_scan'+'.inp', 'r')
        elif forceox == True:
            print('{0}_{1:.3f}.inp\n'.format(file, strain))
            f = open(loc+'{0}_{1:.3f}.inp'.format(file, strain), 'r')
        else:
            print(file+'.inp\n')
            f = open(loc+file+'.inp', 'r')
        for line in f:
            print(line[:-1])
        f.close()
        if forceox == True:
            os.system('rm {2}/{0}_{1:.3f}.inp'.format(file, strain, loc))

    def screen3d(file, loc, aim1, aim2, scanox, forceox):
        aim = []
        try:
            if scanox == True or forceox == True:
                atomsx = read(loc+file+'.xyz')
                num = len(atomsx)
                if aim1 >= num or aim2 >= num:
                    print('Setting Error. The maximum number of this system is {0}.\n\
                        We will automatically change\
                        the number into the maximum number.'.format(num-1))
                if aim1 >= num:
                    aim1 = num - 1
                if aim2 >= num:
                    aim2 = num - 1
                aim = [aim1, aim2]
                print('Atom1: {0}      Atom2: {1}'.format(aim1, aim2))
                view3dchoose(file, loc, choose=aim, width=400, height=400)
            else:
                view3d(file, loc, width=400, height=500)
        except:
            print('Setting Error.')

    def createinp(button):
        output.clear_output()
        output2.clear_output()
        with output:
            out = widgets.interactive_output(guiconvert, {'file': file, 'loc': loc,\
                                                          'method': method, 'bs': bs, 'opt': opt, 'freq': freq, \
                                                          'ts': ts, 'scanox': scanox, 'forceox': forceox, 'aim1': aim1,\
                                                          'aim2': aim2,'stretch': stretch, 'scanstep': scanstep,\
                                                          'strain': strain, 'maxiter': maxiter, 'room': room, \
                                                          'roomox': roomox, 'corenum': corenum, 'coreox': coreox,\
                                                          'e': e, 'state': state, 'mdox': mdox, 'tstep': tstep,\
                                                          'initvel': initvel, 'tsNHC': tsNHC, 'timecon': timecon,\
                                                          'dumpcon': dumpcon, 'runstep': runstep})
            startorca = widgets.Button(description = 'Start ORCA Calculation',\
                                       layout = widgets.Layout(width='99%'),\
                                       style = dict(font_weight='bold', text_decoration='underline'))
            startorca.on_click(intro_status3d)
            display(widgets.VBox([out, widgets.HBox([widgets.Label('ORCA Location:'),\
                                                     orcaloc]), startorca]))
        with output2:
            out2 = widgets.interactive_output(screen3d, {'file': file, 'loc': loc,\
                                                         'aim1': aim1, 'aim2': aim2, 'scanox': scanox, 'forceox': forceox}) 
            display(out2)

    def intro_status3d(button):
        if scanox.value == True or forceox.value == True:
            atomsx = read(loc.value+file.value+'.xyz')
            num = len(atomsx)
            if aim1.value >= num:
                aim1.value = num - 1
            if aim2.value >= num:
                aim2.value = num - 1
        with output:
            button.description = 'Processing...'
        if scanox.value == True:
            orca(file.value+'_scan', orcaloc=orcaloc.value)
        elif forceox.value == True:
            if roomox.value == True and coreox.value == True:
                XyzToInp(file.value,\
                         inpname = '{0}_{1:.3f}'.format(file.value, strain.value),\
                         fileloc=loc.value, method=method.value, basis_set=bs.value,\
                         opt=opt.value, freq=freq.value, external_force=True,\
                         aim = [aim1.value, aim2.value], strain=strain.value,\
                         maxiter=maxiter.value, electron=e.value, state=state.value)
            elif roomox.value == False and coreox.value == False:
                XyzToInp(file.value,\
                         inpname = '{0}_{1:.3f}'.format(file.value, strain.value),\
                         fileloc=loc.value, method=method.value, basis_set=bs.value,\
                         opt=opt.value, freq=freq.value, external_force=True,\
                         aim = [aim1.value, aim2.value], strain=strain.value,\
                         maxiter=maxiter.value, maxcore=room.value*1024,\
                         corenum=corenum.value, electron=e.value, state=state.value)
            else:
                XyzToInp(file.value,\
                         inpname = '{0}_{1:.3f}'.format(file.value, strain.value),\
                         fileloc=loc.value, method=method.value, basis_set=bs.value,\
                         opt=opt.value, freq=freq.value, external_force=True,\
                         aim = [aim1.value, aim2.value], strain=strain.value,\
                         maxiter=maxiter.value, maxcore=room.value*1024,\
                         electron=e.value, state=state.value)
            orca('{0}_{1:.3f}'.format(file.value, strain.value), orcaloc=orcaloc.value)
        else:
            orca(file.value, orcaloc=orcaloc.value)
        with output2:
            if scanox.value == True:
                if status(file.value+'_scan', loc.value).status(figureonly=True,\
                                                                statusonly=True) == 2:
                    display(widgets.HBox([widgets.Valid(description = file.value+'_scan',\
                                                        value=True), widgets.Label('  Calculation completed!')]))
                elif status(file.value+'_scan', loc.value).status(figureonly=True,\
                                                                  statusonly=True) == 4:
                    display(widgets.Valid(description = file.value+'_scan', value=False,\
                                          readout='Calculation aborted!'))
            elif forceox.value == True:
                if status('{0}_{1:.3f}'.format(file.value, strain.value),\
                          loc.value).status(figureonly=True, statusonly=True) == 2:
                    display(widgets.HBox([
                        widgets.Valid(description = '{0}_{1:.3f}'.format(file.value, strain.value),
                                      value=True), widgets.Label('  Calculation completed!')
                    ]))
                elif status('{0}_{1:.3f}'.format(file.value, strain.value),\
                            loc.value).status(figureonly=True, statusonly=True) == 4:
                    display(widgets.Valid(description = '{0}_{1:.3f}'.format(file.value,\
                                                                             strain.value), value=False, readout='Calculation aborted!'))
            else:
                if status(file.value, loc.value).status(figureonly=True,\
                                                        statusonly=True) == 2:
                    display(widgets.HBox([widgets.Valid(description = file.value,\
                                                        value=True), widgets.Label('  Calculation completed!')]))
                elif status(file.value, loc.value).status(figureonly=True,\
                                                          statusonly=True) == 4:
                    display(widgets.Valid(description = file.value,\
                                          value=False, readout='Calculation aborted!'))
            display(introstatus)
        output.clear_output()
        
    def status3d(button):
        output2.clear_output()
        with output2:
            display(widgets.VBox([eox, eoutput, geoox, geooutput,\
                                  widgets.HBox([widgets.Label('Save Name:'), savename, save])]))
    
    def saved(button):
        status(file.value, loc.value).figure(num=savenumchoose.value,\
                                             save=savename.value)

    options = []
    loc = widgets.Text(description = 'Location:', value='./')
    for path in os.listdir(loc.value):
        if os.path.isfile(os.path.join(loc.value, path)):
            a = re.search('.xyz', path)
            if a:
                options.append(path[:-4])
    file = widgets.Combobox(placeholder='Input the file name', options=options,\
                            description = 'xyz File:')
    method = widgets.Text(description = 'Method:')
    bs = widgets.Text(description = 'Basis Set:')
    opt = widgets.ToggleButton(value=False, description = 'Optimisation',\
                               layout = widgets.Layout(width='33.3%'))
    ts = widgets.ToggleButton(value=False, description = 'TS',\
                              layout = widgets.Layout(width='33.3%'))
    freq = widgets.ToggleButton(value=False, description = 'Frequency',\
                                layout = widgets.Layout(width='33.3%'))
    room = widgets.IntText(description = 'Space(GB):', value= 4, min=1,\
                           layout = widgets.Layout(width='50%'))
    roomox = widgets.ToggleButton(value=False, description = 'Default Space')
    corenum=widgets.IntText(description = 'Core:', value=4, min=1,\
                            layout = widgets.Layout(width='50%'))
    coreox = widgets.ToggleButton(value=False, description = 'Default CoreNumber')
    e = widgets.IntText(description = 'Charge:', value=0,\
                        layout = widgets.Layout(width='50%'))
    state = widgets.IntText(description = 'State:', value=1,\
                            layout = widgets.Layout(width='47.5%'))
    scanox = widgets.ToggleButton(value=False, description = 'Scan Settings',\
                                  style = dict(font_weight = 'bold'), layout = widgets.Layout(width='98.5%'))
    forceox = widgets.ToggleButton(value=False,\
                                   description = 'Constant External Force Settings',\
                                   style = dict(font_weight = 'bold'), layout = widgets.Layout(width='98.5%'))
    mdox = widgets.ToggleButton(value=False, description = 'Molecular Dynamics',\
                                style = dict(font_weight = 'bold'), layout = widgets.Layout(width='98.5%'))
    aim1 = widgets.BoundedIntText(description = 'Atom 1:', value=0,\
                                  min=0, max=300, layout = widgets.Layout(width='50%'))
    aim2 = widgets.BoundedIntText(description = 'Atom 2:', value=0,\
                                  min=0, max=300, layout = widgets.Layout(width='47.5%'))
    strain=widgets.BoundedFloatText(description = 'Force(nN):',\
                                    value=0.000, min=0.000, step=0.100)
    stretch=widgets.FloatText(description = 'Distance(Å):',\
                              value=1.000, step=0.100)
    scanstep=widgets.BoundedIntText(description = 'Steps:', value=10, min=0)
    scanout = widgets.interactive_output(scanshow, {'scanox': scanox})
    forceout = widgets.interactive_output(forceshow, {'forceox': forceox})
    eox = widgets.ToggleButton(value=False, description = 'Energy Status',\
                               style = dict(font_weight = 'bold'), layout = widgets.Layout(width='98.5%'))
    geoox = widgets.ToggleButton(value=False, description = 'Structure Status',\
                                 style = dict(font_weight = 'bold'), layout = widgets.Layout(width='98.5%'))
    maxiter=widgets.IntText(description = 'MaxIteration:', value=200,\
                            layout = widgets.Layout(width='98.5%'))
    tstep=widgets.BoundedFloatText(description = 'Timestep:', min=0, value=1.0,\
                                   step=0.100, layout = widgets.Layout(width='47.5%'))
    initvel = widgets.BoundedFloatText(description = 'Initial T (K):', min=0,\
                                       max=2000, value=300, step=0.100, layout = widgets.Layout(width='47.5%'))
    tsNHC = widgets.BoundedFloatText(description = 'NHC T (K):', min=0, max=2000,\
                                     value=300, step=1.000, layout = widgets.Layout(width='47.5%'))
    timecon = widgets.BoundedFloatText(description = 'Time Gap:', min=0,\
                                       value=10.0, step=0.100, layout = widgets.Layout(width='47.5%'))
    dumpcon = widgets.BoundedIntText(description = 'Dump Gap:', min=1,\
                                     value=10, step=1, layout = widgets.Layout(width='98.5%'))
    runstep=widgets.BoundedIntText(description = 'Run Steps:', min=1,\
                                   value=100, max=1000000, step=1, layout = widgets.Layout(width='98.5%'))
    savename = widgets.Text()
    save = widgets.Button(description = 'Save')
    savenumchoose = widgets.IntSlider(min=0, step=1)
    
    eout = widgets.interactive_output(eshow, {'eox': eox})
    geoout = widgets.interactive_output(geoshow, {'geoox': geoox})
    mdout = widgets.interactive_output(mdshow, {'mdox': mdox})
    orcaloc = widgets.Text(value='./')
    
    title1 = widgets.Label('File Details', style = dict(font_weight = 'bold'))
    title2 = widgets.Label('Calculation Methods',\
                           style = dict(font_weight = 'bold'))
    
    box_layout1 = widgets.Layout(display = 'flex', flew_flow = 'column',\
                                 overflow = 'scroll hidden', align_items = 'stretch',\
                                 border = 'solid', height = '650px', width='31%')
    box_layout2 = widgets.Layout(display = 'flex', flew_flow = 'column',\
                                 overflow = 'scroll hidden', align_items = 'stretch',\
                                 border = 'solid', height = '650px', width='33%')
    widgets.jslink((roomox, 'value'), (room, 'disabled'))
    widgets.jslink((coreox, 'value'), (corenum, 'disabled'))
    create = widgets.Button(description = 'Preview ORCA Input File',\
                            style = dict(font_weight = 'bold', text_decoration = 'underline'),\
                            layout = widgets.Layout(width='98.5%'))
    introstatus = widgets.Button(description = 'See Current Calculation Status',\
                                 style = dict(font_weight = 'bold', text_decoration = 'underline'),\
                                 layout = widgets.Layout(width='98.5%'))

    frame1 = widgets.VBox([title1, loc, file, title2,\
                           widgets.HBox([opt, ts, freq]), method, bs,\
                           widgets.HBox([corenum, coreox]), widgets.HBox([room, roomox]),\
                           widgets.HBox([e, state]), maxiter, scanox, scanoutput, forceox,\
                           forceoutput, mdox, mdoutput, create], layout = box_layout1)
    output.layout = box_layout2
    output2.layout = box_layout2
    create.on_click(createinp)
    save.on_click(saved)
    introstatus.on_click(status3d)

    display(widgets.HBox([frame1, output, output2]))