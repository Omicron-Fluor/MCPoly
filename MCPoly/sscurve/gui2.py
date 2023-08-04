import os
import re
import warnings
import sys
import time as t
import pandas as pd
import seaborn as sns
import ipywidgets as widgets
import matplotlib.pyplot as plt
from IPython.display import display
from .YModulus import YModulus
from .single import single

mydir = os.path.dirname( __file__ )
viewdir = os.path.join(mydir, '..', 'orcaset')
sys.path.append(viewdir)
from view3dchoose import view3dchoose

mydir = os.path.dirname( __file__ )
viewdir = os.path.join(mydir, '..', 'orcaset')
sys.path.append(viewdir)
from XyzToInp import XyzToInp

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from status import status

def gui2():
    """
        The method to use GUI to assemble mechanical datas and get the multiple stress-strain curve on a chart in handy!
        You will know how to do as long as you get into the GUI platform.

        TIPS: 1. Make sure all your created .xyz files from orcaset.ssorca or orcaset.ssgui is on the direct folder, or the programme will get error.
              2. You'd better have relevant .inp or .out file on the direct folder to automatically set the site of external force. Though GUI has manual auxiliary, it will cost you much more time.
              3. The visualisation will also risk the flexibility of the task, to make you programme more flexible, please directly use the code.
    """
    def multiple_gui(polymers, loc = './', savedata = False, savefig=False, pngname = '', csvname = ''):
        opath = os.getcwd()
        os.chdir(loc)
        try:
            f = open('ot.txt', 'x')
        except:
            f = open('ot.txt', 'w')
        
        name = []
        i = 0
        
        site = []
        
        for polymer in polymers:
            polymer = polymer + '_Result.txt'
            try:
                datum = pd.read_csv(polymer)
            except:
                warnings.warn("'{0}' is not found.".format(polymer))
                continue
            datum = datum[:].sort_values(by=['Distance(A)'])
            for num in datum['Force(nN)'][:-1]:
                if num not in site:
                    site.append(num)
        site.sort()
        f.write('Force(nN)\n')
        for num in site:
            f.write('{0:.3f}\n'.format(num))
        f.close()
        ot = pd.read_csv('ot.txt')
        for polymer in polymers:
            polymer = polymer + '_Result.txt'
            try:
                datum = pd.read_csv(polymer)
            except:
                warnings.warn('{0} is not the file aboue stress-strain curve.'.format(polymer))
                continue
            datum = datum[:].sort_values(by = ['Distance(A)'])
            ax = plt.plot(datum['Force(nN)'][:-1], datum['Distance(A)'][:-1], '*-')
            datum = datum.rename(
                columns = {'Distance(A)':'Distance(A)({0})'.format(polymer[:-11])})
            ot = pd.merge_ordered(ot, datum[:-1], fill_method = "ffill",\
                                  left_by = "Force(nN)")
            name.append(polymer[:-11])
            i = i + 1
    
        ot = ot.set_index('Force(nN)')
        if savedata == True or savefig == True:
            if savedata == True:
                ot.to_csv('{0}.csv'.format(csvname))
        else:
            output3.clear_output()
            with output3:
                display(ot)
                display(widgets.HBox([widgets.Label('CSV File Name:'), name2]))
                display(save2)
        
        plt.xlabel('Force(nN)')
        plt.ylabel('Distance(A)')
        plt.legend(name, bbox_to_anchor = [1.1, 0.9])
        _ = plt.title('Distance Relationship')
        if savedata == True or savefig == True:
            if savefig == True:
                plt.savefig('{0}.png'.format(pngname), bbox_inches = 'tight')
        else:
            with output2:
                plt.show()
                display(widgets.HBox([widgets.Label('PNG File Name:'), name1]))
                display(save1)
    
        os.system('rm ot.txt')
        os.chdir(opath)
    
    print('TIPS: Multiple values can be selected with\
        Shift and/or Ctrl (or Command) pressed\
        and mouse clicks or arrow keys.')
    
    output2 = widgets.Output()
    output3 = widgets.Output()
    judgeoutput = widgets.Output()
    
    def geoshow(aim1, aim2, geoi, files):
        output2.clear_output()
        with output2:
            try:
                print("'{0}.xyz' Geometry Structure:".format(files[geoi]))
            except:
                return None
            aim = [aim1, aim2]
            print('Atom1: {0}      Atom2: {1}'.format(aim1, aim2))
            try:
                view3dchoose(files[geoi]+'_0.000', loc.value,\
                             choose=aim, width=450, height=400)
            except:
                view3dchoose(files[geoi]+'_0', loc.value, choose=aim, width=450, height=400)
    
    def filejudge(files, loc):
        judgeoutput.clear_output()
        output2.clear_output()
        output3.clear_output()
        opath = os.getcwd()
        os.chdir(loc)
        #geoi = widgets.IntText(value=0)
        t0 = 0
        try:
            geoi.value = 0
            for file in files:
                t0 = 0
                for path in os.listdir(loc):
                    if os.path.isfile(os.path.join(loc, path)):
                        a = re.match('{0}_Result.txt'.format(file), path)
                        if a:
                            output2.clear_output()
                            geoi.value = geoi.value + 1
                            t0 = 1
                            break
                if t0 == 0:
                    output2.clear_output()
                    single(file, loc = loc).autocurve2(savefig=False)
                    geoi.value = geoi.value + 1
                #single(file, loc = loc).autocurve(savefig=False)
                #geoi.value = geoi.value + 1
        except:
            xtitle = widgets.Label(
                'Extra Information about {0}:'.format(files[geoi.value]),\
                style = dict(font_weight = 'bold'))
            #aim1 = widgets.BoundedIntText(description = 'Atom 1:',\
            #                              value=0, min=0, layout = widgets.Layout(width='48%'))
            #aim2 = widgets.BoundedIntText(description = 'Atom 2:',\
            #                              value=0, min=0, layout = widgets.Layout(width='48%'))
            aim1.value = 0
            aim2.value = 0
            OK = widgets.Button(description = 'OK')
            with judgeoutput:
                display(widgets.VBox([xtitle, widgets.HBox([aim1, aim2]), OK]))
            with output2:
                #print("'{0}' Geometry Structure:".format(file))
                #geoout = widgets.interactive_output(geoshow, {'aim1': aim1, 'aim2': aim2, 'geoi': geoi})
                OK.on_click(OKss)
        os.chdir(opath)
    
    def OKss(button):
        opath = os.getcwd()
        os.chdir(loc.value)
        #with output3:
        #    print(geoi.value, aim1.value, aim2.value)
        single(files.value[geoi.value],\
               loc = loc.value).curve(aim1.value, aim2.value, savefig=False)
        judgeoutput.clear_output()
        output2.clear_output()
        geoi.value = geoi.value + 1
        try:
            if geoi.value >= len(files.value):
                return None
            for i in range(geoi.value, len(files.value)):
                for path in os.listdir(loc.value):
                    if os.path.isfile(os.path.join(loc.value, path)):
                        a = re.match(files.value[geoi.value]+'_Result.rxt', path)
                        if a:
                            geoi.value = geoi.value + 1
                            return None
                single(files.value[i], loc=loc.value).autocurve2(savefig=False)
                geoi.value = geoi.value + 1
        except:
            xtitle = widgets.Label(
                'Extra Information about {0}:'.format(files.value[geoi.value]),\
                style = dict(font_weight = 'bold'))
            #aim1 = widgets.BoundedIntText(description = 'Atom 1:',\
            #                              value=0, min=0, layout = widgets.Layout(width='48%'))
            #aim2 = widgets.BoundedIntText(description = 'Atom 2:',\
            #                              value=0, min=0, layout = widgets.Layout(width='48%'))
            OK = widgets.Button(description = 'OK')
            with judgeoutput:
                display(widgets.VBox([xtitle, widgets.HBox([aim1, aim2]), OK]))
            with output2:
                #print("'{0}.xyz' Geometry Structure:".format(files.value[geoi.value]))
                #geoout = widgets.interactive_output(geoshow,\
                #                                    {'aim1': aim1, 'aim2': aim2, 'geoi': geoi})
                #display(geoout)
                OK.on_click(OKss)
        os.chdir(opath)
    
    def intro_result(button):
        for file in files.value:
            t0 = 0
            try:
                t0 = 0
                for path in os.listdir(loc.value):
                    if os.path.isfile(os.path.join(loc.value, path)):
                        a = re.match(file+'_Result.rxt', path)
                        if a:
                            t0 = 1
                            break
                if t0 == 0:
                    single(file, loc=loc.value).autocurve2(savefig=False)
            except:
                pass
    #            xtitle = widgets.Label('Extra Information about {0}:'.format(file),\
    #                                   style = dict(font_weight = 'bold'))
    #            aim1 = widgets.BoundedIntText(description = 'Atom 1:', value=0,\
    #                                          min=0, layout = widgets.Layout(width='48%'))
    #            aim2 = widgets.BoundedIntText(description = 'Atom 2:', value=0,\
    #                                          min=0, layout = widgets.Layout(width='48%'))
    #            OK = widgets.Button(description = 'OK')
    #            with output:
    #                display(widgets.VBox[xtitle, widgets.HBox([aim1, aim2]), OK])
    #            geoout = widgets.interactive_output(geoshow, {'aim1': aim1, 'aim2': aim2})
    #            with outout2:
    #                display(geoout)
    #            OK.on_click(OKss)
    #            t.sleep(15)
        multiple_gui(files.value, loc=loc.value)
        
    def pngsave(button):
        multiple_gui(files.value, loc=loc.value, savefig=True, pngname=name1.value)
    
    def csvsave(button):
        multiple_gui(files.value, loc=loc.value, savedata=True, csvname=name1.value)
    
    options = []
    chooses = []
    title1 = widgets.Label('File Settings:', style = dict(font_weight = 'bold'))
    loc = widgets.Text(description = 'Location:', value='./')
    for path in os.listdir(loc.value):
        if os.path.isfile(os.path.join(loc.value, path)):
            a = re.search('.xyz', path)
            if a:
                b = re.search('_', path)
                if b:
                    c = re.search('.', path[:-4])
                    if c:
                        fragment = re.split('_', path)
                        if fragment[0] not in options:
                            options.append(fragment[0])
                        if fragment[1][:-4] not in chooses:
                            chooses.append(fragment[1][:-4])
    options.sort()
    files = widgets.SelectMultiple(options=options, description = 'Files Series:')
    #ym = widgets.Button(description = "Young's Modulus",\
    #                    layout = widgets.Layout(width='100%'))
    resultshow = widgets.Button(description = 'Show Chart and table',\
                                layout = widgets.Layout(width='100%'),\
                                style = dict(font_weight = 'bold', text_decoration = 'underline'))
    save1 = widgets.Button(description = 'Save',\
                           layout = widgets.Layout(width='100%'))
    save2 = widgets.Button(description = 'Save',\
                           layout = widgets.Layout(width='100%'))
    name1 = widgets.Text(value = 'Result', layout = widgets.Layout(width='60%'))
    name2 = widgets.Text(value = 'Result', layout = widgets.Layout(width='60%'))
    aim1 = widgets.BoundedIntText(description = 'Atom 1:', value=0, min=0,\
                                  layout = widgets.Layout(width='48%'))
    aim2 = widgets.BoundedIntText(description = 'Atom 2:', value=0, min=0,\
                                  layout = widgets.Layout(width='48%'))
    OK = widgets.Button(description = 'OK')
    geoi = widgets.IntText(value = 0)
    
    output = widgets.VBox([title1, loc, files, resultshow])
    
    resultshow.on_click(intro_result)
    save1.on_click(pngsave)
    save2.on_click(csvsave)
    
    box_layout1 = widgets.Layout(display='flex', flew_flow='column',\
                                 overflow='scroll hidden', align_items='stretch',\
                                 border='solid', height='550px', width='26%')
    box_layout2 = widgets.Layout(display='flex', flew_flow='column',\
                                 overflow='scroll hidden', align_items='stretch',\
                                 border='solid', height='550px', width='37%')
    judgeoutput.layout = widgets.Layout(display='flex', flew_flow='column',\
                                        overflow='scroll hidden', align_items='stretch',\
                                        height='325px', width='100%')
    
    filejudgeout = widgets.interactive_output(filejudge,\
                                              {'files': files, 'loc': loc})
    geoout = widgets.interactive_output(geoshow, {'aim1': aim1, 'aim2': aim2,\
                                                  'geoi': geoi, 'files': files})
    
    output = widgets.VBox([title1, loc, files, resultshow])
    output2.layout = box_layout2
    output3.layout = box_layout2
    frame1 = widgets.VBox([output, judgeoutput])
    frame1.layout = box_layout1
    
    resultshow.on_click(intro_result)
    save1.on_click(pngsave)
    save2.on_click(csvsave)
    
    display(widgets.HBox([frame1, output2, output3]))