import sys
import os
import ipywidgets as widgets

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from status import status
from gui import gui as gui3

mydir = os.path.dirname( __file__ )
orcadir = os.path.join(mydir, '..')#, 'orcaset')
sys.path.append(orcadir)
from orcaset.gui import gui as gui1
from sscurve.gui import gui as gui2
from moldraw.gui import gui as gui4
from sscurve.gui2 import gui2 as gui5

mydir = os.path.dirname( __file__ )
orcadir = os.path.join(mydir, '..', 'orcaset')
sys.path.append(orcadir)
from ssgui import ssgui

mydir = os.path.dirname( __file__ )
orcadir = os.path.join(mydir, '..', 'orcaset')
sys.path.append(orcadir)
from mgui import mgui

#mydir = os.path.dirname( __file__ )
#ssdir = os.path.join(mydir, '..', 'sscurve')
#sys.path.append(ssdir)
#from gui import gui as gui2

def vis():
    """
        The method to summon all four GUI functions.
        You can see how to use by input 'MCPoly.orcaset.gui?'(Normal ORCA), 'MCPoly.orcaset.ssgui?' (Mechanical ORCA), 'MCPoly.orcaset.mgui?' (Multiple ORCA) and 'MCPoly.sscurve.gui?' (Stress-Strain Curve)
    """
    screen1 = widgets.Output()
    screen2 = widgets.Output()
    screen3 = widgets.Output()
    screen4 = widgets.Output()
    screen5 = widgets.Output()
    screen6 = widgets.Output()
    screen7 = widgets.Output()
    
    with screen1:
        gui1()
    with screen2:
        ssgui()
    with screen3:
        mgui()
    with screen4:
        gui2()
    with screen5:
        gui5()
    with screen6:
        gui3()
    with screen7:
        gui4()
    
    tab_contents = ['Normal ORCA', 'Mechanical ORCA', 'Multiple ORCA', 'Stress-Strain Curve', 'Distance Curve', 'Energy Diagram', 'Molecule Designer']
    tab = widgets.Tab(style=dict(font_weight='bold'))
    #children = [widgets.Label(name,\
    #            style=dict(font_weight='bold')) for name in tab_contents]
    children = [screen1, screen2, screen3, screen4, screen5, screen6, screen7]
    tab.children = children
    tab.titles = tab_contents
    tab.layout = widgets.Layout(layout=widgets.Layout(width='50%'))
    display(tab)