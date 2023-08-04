import os
import re
import sys
import time as t
import ipywidgets as widgets
from IPython.display import display

mydir = os.path.dirname( __file__ )
orcadir = os.path.join(mydir, '..', 'orcaset')
sys.path.append(orcadir)
from XyzToInp import XyzToInp
from orca import orca

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from status import status

mydir = os.path.dirname( __file__ )
viewdir = os.path.join(mydir, '..', 'view3d')
sys.path.append(viewdir)
from view3d import view3d

def mgui():
    """
        The method to use GUI to set ORCA input files and calculate multiple files consecutively by ORCA.
        You will know how to do as long as you get into the GUI platform.

        TIPS: 1. To see the status of your optimisation,  you can set another window and use status.figuretraj .
              2. This visualisation will greatly restrict the flexibility of the task,  to make you programme more flexible,  please directly use the code.
    """

    print('TIPS: Multiple values can be selected with Shift and/or Ctrl (or Command) pressed and mouse clicks or arrow keys.')
    
    output = widgets.Output()
    output2 = widgets.Output()
    scanoutput = widgets.Output()
    forceoutput = widgets.Output()
    eoutput = widgets.Output()
    geooutput = widgets.Output()
    case = widgets.Output()
    
    def geoshow(files, loc, method, bs):
        output.clear_output()
        with output:
            if files == []:
                print('Please input the file to show the geometry structure on this screen.')
            for file in files:
                try:
                    print('File: {0}.xyz'.format(file))
                    view3d(file, loc, width=450, height=400)
                except:
                    print("'{0}.xyz' is not found.".format(file))
            if method != '' and bs != '':
                display(widgets.HBox([widgets.Label('ORCA Location:'),\
                                      orcaloc]), startorca)
    
    def intro_status3d(button):
        ii = 0
        for file in files.value:
            ii = ii + 1
            if roomox.value == False and coreox.value == False:
                XyzToInp(file, fileloc=loc.value, method=method.value,\
                         basis_set=bs.value, opt=opt.value, freq=freq.value,\
                         maxiter=maxiter.value, maxcore=room.value*1024,\
                         corenum=corenum.value, electron=e.value, state=state.value)
            elif roomox.value == False and coreox.value == True:
                XyzToInp(file, fileloc=loc.value, method=method.value,\
                         basis_set=bs.value, opt=opt.value, freq=freq.value,\
                         maxiter=maxiter.value, maxcore=room.value*1024,\
                         electron=e.value, state=state.value)
            elif roomox.value == True and coreox.value == False:
                XyzToInp(file, fileloc=loc.value, method=method.value,\
                         basis_set=bs.value, opt=opt.value, freq=freq.value,\
                         maxiter=maxiter.value, corenum=corenum.value,\
                         electron=e.value, state=state.value)
            elif roomox.value == True and coreox.value == True:
                XyzToInp(file, fileloc=loc.value, method=method.value,\
                         basis_set=bs.value, opt=opt.value, freq=freq.value,\
                         maxiter=maxiter.value, electron=e.value, state=state.value)
            with case:
                display(widgets.HBox([widgets.Label(file), widgets.Label('  Processing...'),\
                                      widgets.Label('[{0}]'.format(t.ctime(t.time())))]))
            orca(file, orcaloc.value, loc.value, loc.value)
            s = status(file, loc.value).status(statusonly=True)
            if s == 4:
                case.clear_output()
                with case:
                    display(widgets.Valid(description = '     ' + file, value=False, \
                            readout = 'Retrying...'),\
                            widgets.Label('[{0}]'.format(t.ctime(t.time()))))
                orca(file, orcaloc.value, loc.value, loc.value)
                s = status(file, loc.value).status(statusonly=True)
                if s == 4:
                    case.clear_output()
                    with output2:
                        display(widgets.Valid(description = '     ' + file, value=False, \
                                readout = 'Error!'),\
                                widgets.Label('[{0}]'.format(t.ctime(t.time()))))
            elif s == 2:
                case.clear_output()
                with output2:
                    display(widgets.HBox([widgets.Valid(description = '     ' + file,\
                                                        value=True),\
                                          widgets.Label('Completed.'),\
                                          widgets.Label('[{0}]'.format(t.ctime(t.time())))]))
            if ii == 1:
                output.clear_output()
            with output:
                print('File: {0}.xyz    Energy: {1:.6f} Hartree'.format(file,\
                                                                        status(file, loc.value).status(figureonly = True)[-1]))
                view3d(file, loc.value, width=450, height=400)
    
    options = []
    loc = widgets.Text(description = 'Location:', value='./')
    for path in os.listdir(loc.value):
        if os.path.isfile(os.path.join(loc.value,  path)):
            a = re.search('.xyz', path)
            if a:
                b = re.search('_trj', path)
                if b:
                    pass
                else:
                    options.append(path[:-4])
    options.sort()
    files = widgets.SelectMultiple(options=options, description = 'xyz Files:')
    method = widgets.Text(description = 'Method:')
    bs = widgets.Text(description = 'Basis Set:')
    opt = widgets.ToggleButton(value=False, description = 'Optimisation',\
                               layout = widgets.Layout(width='50%'))
    freq = widgets.ToggleButton(value=False, description = 'Frequency',\
                                layout = widgets.Layout(width='50%'))
    room = widgets.IntText(description = 'Space(GB):', value=4, min=1,\
                           layout = widgets.Layout(width='50%'))
    roomox = widgets.ToggleButton(value=False, description = 'Default Space')
    corenum = widgets.IntText(description = 'Core:', value=4, min=1,\
                              layout = widgets.Layout(width='50%'))
    coreox = widgets.ToggleButton(value=False, description = 'Default CoreNumber')
    e = widgets.IntText(description = 'Charge:', value=0,\
                        layout = widgets.Layout(width='50%'))
    state = widgets.IntText(description = 'State:', value=1,\
                            layout = widgets.Layout(width='47.5%'))
    geoox = widgets.ToggleButton(value=False, description = 'Structure Status',\
                                 style = dict(font_weight='bold'), layout = widgets.Layout(width='98.5%'))
    startorca = widgets.Button(description = 'Start ORCA Calculation',\
                               layout = widgets.Layout(width='99%'),\
                               style = dict(font_weight = 'bold', text_decoration = 'underline'))
    geoout = widgets.interactive_output(geoshow, {'files': files, 'loc': loc,\
                                                  'method': method, 'bs': bs})
    orcaloc = widgets.Text(value = './')
    maxiter = widgets.IntText(description = 'MaxIteration:', value=200,\
                              layout = widgets.Layout(width='98.5%'))
    
    title1 = widgets.Label('File Details', style = dict(font_weight='bold'))
    title2 = widgets.Label('Calculation Methods',\
                           style = dict(font_weight='bold'))
    
    box_layout1  =  widgets.Layout(display='flex', flew_flow='column',\
                                   overflow = 'scroll hidden', align_items = 'stretch',\
                                   border='solid',  height='650px',  width='31%')
    box_layout2  =  widgets.Layout(display='flex',  flew_flow='column',\
                                   overflow = 'scroll hidden', align_items = 'stretch',\
                                   border='solid',  height='650px',  width='33%')
    widgets.jslink((roomox, 'value'), (room, 'disabled'))
    widgets.jslink((coreox, 'value'), (corenum, 'disabled'))
    #introstatus = widgets.Button(description = 'See Current Calculation Status',\
    #                             style = dict(font_weight = 'bold', text_decoration = 'underline'), 
    #                  layout = widgets.Layout(width='98.5%'))

    frame1 = widgets.VBox([title1, loc, files, title2, widgets.HBox([opt, freq]),\
                           method, bs, widgets.HBox([corenum, coreox]), \
                           widgets.HBox([room, roomox]), maxiter, widgets.HBox([e, state])],\
                          layout = box_layout1)
    frame3 = widgets.VBox([output2, case])
    output.layout = box_layout2
    frame3.layout = box_layout2
    startorca.on_click(intro_status3d)

    display(widgets.HBox([frame1, output, frame3]))