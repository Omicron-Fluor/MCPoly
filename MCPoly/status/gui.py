import os
import re
import sys
import warnings
import time as t
import numpy as np
import pandas as pd
import seaborn as sns
import ipywidgets as widgets
import matplotlib.pyplot as plt

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from status import status

def gui():
    """
        The method to draw the reaction energy diagram.
        You will know how to do as long as you get into the GUI platform.
    
        TIPS: 1. You can combine the table with the function with 'Result Combination'.
              2. You can skip the a file name if you need.
              3. For some limitation, the chart created on this function might not very good.
              4. The 'Result Combination' can only create line chart.
              5. The visualisation will also risk the flexibility of the task, to make you programme more flexible, please directly use the code.
    """
    
    def fig(button):
        opath = os.getcwd()
        os.chdir(loc.value)
        os.system('cp {0}.png {1}.png'.format(figoldname.value, figname.value))
        os.system('rm {0}.png'.format(figoldname.value))
        os.system('rm ot.png'.format(csv6.value))
        if delete.value == True:
            os.system('rm {0}.png'.format(csv1.value))
            os.system('rm {0}.png'.format(csv2.value))
            os.system('rm {0}.png'.format(csv3.value))
            os.system('rm {0}.png'.format(csv4.value))
            os.system('rm {0}.png'.format(csv5.value))
            os.system('rm {0}.png'.format(csv6.value))
        os.chdir(opath)
        output2.clear_output()
        controloutput1.clear_output()
        file1.value = ''
        file2.value = ''
        file3.value = ''
        file4.value = ''
        file5.value = ''
        file6.value = ''
        file7.value = ''
        file8.value = ''
        file9.value = ''
        file10.value = ''
        csv1.value = ''
        csv2.value = ''
        csv3.value = ''
        csv4.value = ''
        csv5.value = ''
        csv6.value = ''
        figname.value = 'Untitled'
    
    def data(button):
        opath = os.getcwd()
        os.chdir(loc.value)
        os.system('cp {0}.csv {1}.csv'.format(dataoldname.value, dataname.value))
        os.system('rm {0}.csv'.format(dataoldname.value))
        os.system('rm ot.csv'.format(dataoldname.value))
        if delete.value == True:
            os.system('rm {0}.csv'.format(csv1.value))
            os.system('rm {0}.csv'.format(csv2.value))
            os.system('rm {0}.csv'.format(csv3.value))
            os.system('rm {0}.csv'.format(csv4.value))
            os.system('rm {0}.csv'.format(csv5.value))
            os.system('rm {0}.csv'.format(csv6.value))
        os.chdir(opath)
        output3.clear_output()
        file1.value = ''
        file2.value = ''
        file3.value = ''
        file4.value = ''
        file5.value = ''
        file6.value = ''
        file7.value = ''
        file8.value = ''
        file9.value = ''
        file10.value = ''
        csv1.value = ''
        csv2.value = ''
        csv3.value = ''
        csv4.value = ''
        csv5.value = ''
        csv6.value = ''
        dataname.value = 'Untitled'
    
    def pattern1(fig_pattern, datatype):
        output2.clear_output()
        output3.clear_output()
        opath = os.getcwd()
        os.chdir(loc.value)
        sonbox = []
        files = [file1.value, file2.value, file3.value, file4.value, file5.value,\
                 file6.value, file7.value, file8.value, file9.value, file10.value]
        if files == ['', '', '', '', '', '', '', '', '', '']:
            return None
        for file in files:
            try:
                if datatype == 'Energy':
                    data2 = status(file).converged_energy()
                elif datatype == 'Gibbs':
                    data2 = status(file).gibbs()
                sonbox.append(data2)
            except:
                sonbox.append(49999)
        num = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        Ebar = np.array(sonbox)
        Ebar = Ebar - sonbox[0]
        Ebar = Ebar * 627.509
        Ebar.tolist()
        for i in range(len(Ebar)):
            if Ebar[i] >= 40000:
                Ebar[i] = None
        Emax = -99999
        for E in Ebar:
            try:
                if E > Emax:
                    Emax = E
            except:
                pass
        Emin = min(Ebar)
        gap = Emax - Emin
        if fig_pattern == 'line':
            ax = plt.plot(num, Ebar, '*-')
        elif fig_pattern == 'bar':
            ax = plt.bar(num, Ebar)
        plt.xlabel('Step')
        if datatype == 'Energy':
            plt.ylabel('∆E(kcal/mol)')
        elif datatype == 'Gibbs':
            plt.ylabel('∆G(kcal/mol)')
        for i in range(len(Ebar)):
            if files[i] == '':
                continue
            if Emax - Ebar[i] <= 5.0:
                plt.text(num[i]-0.5, Ebar[i] - 0.05*gap,\
                         '{0}\n{1:.1f}'.format(files[i], Ebar[i]))
            elif Emax - Ebar[i] > 5.0:
                plt.text(num[i]-0.2, Ebar[i] + 0.02*gap,\
                         '{0}\n{1:.1f}'.format(files[i], Ebar[i]))
        #plt.ylim(min(Ebar)-4, max(EBar)+4)
        plt.savefig('{0}.png'.format(figname.value), bbox_inches = 'tight')
        with output2:
            plt.show()
        #with controloutput1:
            #display(widgets.HBox([widgets.Label('Figure Pattern: '), fig_pattern]))
            #display(widgets.VBox([widgets.HBox(
            #    [widgets.Label('Curve Name: '), figname]), savefig]))
        try:
            f = open('{0}.csv'.format(dataname.value), 'x')
        except:
            f = open('{0}.csv'.format(dataname.value), 'w')
        f.write('Step, ')
        f.write('File, ')
        if datatype == 'Energy':
            f.write('∆E')
        elif datatype == 'Gibbs':
            f.write('∆G')
        #if hartree == True:
        #    f.write('(Hartree)')
        #else:
        f.write('(kcal/mol)\n')
        for i in range(len(num)):
            if Ebar[i] == None:
                f.write('{0}, , \n'.format(num[i]))
            else:
                f.write('{0}, {1}, {2}\n'.format(num[i], files[i], Ebar[i]))
        f.close()
        with output3:
            csv = pd.read_csv('{0}.csv'.format(dataname.value))
            display(csv)
        with output3:
            display(widgets.VBox([widgets.HBox([widgets.Label('CSV Table Name'),\
                                                dataname]), savedata]))
        dataoldname.value = dataname.value
        figoldname.value = figname.value
        os.chdir(opath)
        
    def shows(button):
        controloutput1.clear_output()
        output2.clear_output()
        output3.clear_output()
        opath = os.getcwd()
        os.chdir(loc.value)
        with controloutput1:
            display(widgets.HBox([widgets.Label('Figure Pattern: '), fig_pattern]))
            display(widgets.VBox([widgets.HBox([widgets.Label('Curve Name'),\
                                                figname]), savefig]))
        dataoldname.value = dataname.value
        figoldname.value = figname.value
        dataout = widgets.interactive_output(pattern1,\
                                             {'fig_pattern': fig_pattern, 'datatype': datatype})
        os.chdir(opath)
    
    def shows2(button):
        controloutput1.clear_output()
        output2.clear_output()
        output3.clear_output()
        opath = os.getcwd()
        os.chdir(loc.value)
        sonbox2 = []
        csvs = [csv1.value, csv2.value, csv3.value,\
                csv4.value, csv5.value, csv6.value]
        allcsv = []
        site = []
        name = []
        try:
            f = open('ot.csv', 'x')
        except:
            f = open('ot.csv', 'w')
        
        for csv in csvs:
            csv = csv + '.csv'
            try:
                datum = pd.read_csv(csv)
            except:
                warnings.warn("'{0}' is not found.".format(csv))
                continue
            try:
                datum = datum[:].sort_values(by = ['Step'])
            except:
                warnings.warn("'{0}' is not the file aboue stress-strain curve.".format(csv))
                continue
            for num in datum['Step']:
                if num not in site:
                    site.append(num)
        site.sort()
        f.write('Step\n')
        for num in site:
            f.write('{0:.3f}\n'.format(num))
        f.close()
        
        figoldname.value = 'ot'
        dataoldname.value = 'ot'
        
        ot = pd.read_csv('ot.csv')
        for csv in csvs:
            csv = csv + '.csv'
            xo = 0
            try:
                r = open(csv, 'r')
            except:
                continue
            for line in r:
                print(line)
                z = re.search('G', line)
                if z:
                    xo = 1
                    break
            r.close()
            try:
                datum = pd.read_csv(csv)
            except:
                continue
            datum = datum[:].sort_values(by = ['Step'])
            if xo == 0:
                ax = plt.plot(datum['Step'], datum['∆E(kcal/mol)'], '*-')
                #elif fig_pattern == 'bar':
                #    
                datum = datum.rename(
                    columns = {'∆E(kcal/mol)':'∆E_{0}(kcal/mol)'.format(csv[:-4])})
            elif xo == 1:
                ax = plt.plot(datum['Step'], datum['∆G(kcal/mol)'], '*-')
                datum = datum.rename(
                    columns = {'∆G(kcal/mol)':'∆G_{0}(kcal/mol)'.format(csv[:-4])})
            datum = datum.rename(columns = {'File':'File_{0}'.format(csv[:-4])})
            #if hartree == True:
            ot = pd.merge_ordered(ot, datum, fill_method = "ffill" , left_by = "Step")
            name.append(csv[:-4])
        
        with output3:
            ot = ot.set_index('Step')
            display(ot)
            ot.to_csv('ot.csv')
            display(widgets.VBox([widgets.HBox([widgets.Label('CSV Table Name'),\
                                                dataname]), savedata,\
                                  widgets.HBox([widgets.Label('Delete all combined files'), delete])]))
        
        plt.xlabel('Steps')
        plt.ylabel('Energy(kcal/mol)')
        plt.legend(name, bbox_to_anchor = [1.1, 0.9])
        with output2:
            plt.savefig('ot.png', bbox_inches = 'tight')
            plt.show()
        with controloutput1:
            #display(widgets.HBox([widgets.Label('Figure Pattern: '), fig_pattern]))
            display(widgets.VBox([widgets.HBox([widgets.Label('Curve Name'), figname]),\
                                  savefig, widgets.HBox([widgets.Label('Delete all combined files'),\
                                                         delete])]))
            
    output = widgets.Output()
    output2 = widgets.Output()
    output3 = widgets.Output()
    output4 = widgets.Output()
    controloutput1 = widgets.Output()
    
    loc = widgets.Text(description = 'Location:', value='./')
    options = []
    title1 = widgets.Label('File Settings:', style = dict(font_weight='bold'))
    loc = widgets.Text(description = 'Location:', value='./')
    for path in os.listdir(loc.value):
        if os.path.isfile(os.path.join(loc.value, path)):
            a = re.search('.out', path)
            if a:
                options.append(path[:-4])
    options.sort()
    
    file1 = widgets.Combobox(placeholder = 'Input the file name', options=options,\
                             description = 'ORCA File 1', layout = widgets.Layout(width='98.5%'))
    file2 = widgets.Combobox(placeholder = 'Input the file name', options=options,\
                             description = 'ORCA File 2', layout = widgets.Layout(width='98.5%'))
    file3 = widgets.Combobox(placeholder = 'Input the file name', options=options,\
                             description = 'ORCA File 3', layout = widgets.Layout(width='98.5%'))
    file4 = widgets.Combobox(placeholder = 'Input the file name', options=options,\
                             description = 'ORCA File 4', layout = widgets.Layout(width='98.5%'))
    file5 = widgets.Combobox(placeholder = 'Input the file name', options=options,\
                             description = 'ORCA File 5', layout = widgets.Layout(width='98.5%'))
    file6 = widgets.Combobox(placeholder = 'Input the file name', options=options,\
                             description = 'ORCA File 6', layout = widgets.Layout(width='98.5%'))
    file7 = widgets.Combobox(placeholder = 'Input the file name', options=options,\
                             description = 'ORCA File 7', layout = widgets.Layout(width='98.5%'))
    file8 = widgets.Combobox(placeholder = 'Input the file name', options=options,\
                             description = 'ORCA File 8', layout = widgets.Layout(width='98.5%'))
    file9 = widgets.Combobox(placeholder = 'Input the file name', options=options,\
                             description = 'ORCA File 9', layout = widgets.Layout(width='98.5%'))
    file10 = widgets.Combobox(placeholder = 'Input the file name',\
                              options=options,\
             description = 'ORCA File 10', layout = widgets.Layout(width='98.5%'))
    
    options2 = []
    title2 = widgets.Label('Combination Settings:',\
                           style = dict(font_weight='bold'))
    loc = widgets.Text(description='Location:', value='./')
    for path in os.listdir(loc.value):
        if os.path.isfile(os.path.join(loc.value, path)):
            a = re.search('.csv', path)
            if a:
                options2.append(path[:-4])
    options2.sort()
    csv1 = widgets.Combobox(placeholder = 'Input the file name', options=options2,\
             description = 'CSV File 1', layout = widgets.Layout(width='98.5%'))
    csv2 = widgets.Combobox(placeholder = 'Input the file name', options=options2,\
             description = 'CSV File 2', layout = widgets.Layout(width='98.5%'))
    csv3 = widgets.Combobox(placeholder = 'Input the file name', options=options2,\
             description = 'CSV File 3', layout = widgets.Layout(width='98.5%'))
    csv4 = widgets.Combobox(placeholder = 'Input the file name', options=options2,\
             description = 'CSV File 4', layout = widgets.Layout(width='98.5%'))
    csv5 = widgets.Combobox(placeholder = 'Input the file name', options=options2,\
             description = 'CSV File 5', layout = widgets.Layout(width='98.5%'))
    csv6 = widgets.Combobox(placeholder = 'Input the file name', options=options2,\
             description = 'CSV File 6', layout = widgets.Layout(width='98.5%'))
    
    showcurve = widgets.Button(description = 'Create Curve and Table',\
                               layout = widgets.Layout(width='98.5%'),\
                               style = dict(font_weight='bold', text_decoration='underline'))
    showcurve2 = widgets.Button(description = 'Create Curve and Table',\
                                layout = widgets.Layout(width='98.5%'),\
                                style = dict(font_weight='bold', text_decoration='underline'))
    datatype = widgets.Dropdown(description = 'Data Need:',\
                                options = ['Energy', 'Gibbs'], value='Energy',\
                                layout = widgets.Layout(width='98.5%'))
    fig_pattern = widgets.Dropdown(options = ['line', 'bar'], value='line',\
                                   layout = widgets.Layout(width='80%'))
    dataname = widgets.Text(value='Untitled')
    dataoldname = widgets.Text()
    savedata = widgets.Button(description = 'Save CSV Table',\
                              layout = widgets.Layout(width='98.5%'))
    figname = widgets.Text(value='Untitled')
    figoldname = widgets.Text()
    savefig = widgets.Button(description = 'Save Curve',\
                             layout = widgets.Layout(width='98.5%'))
    delete = widgets.Checkbox()
    
    showcurve.on_click(shows)
    showcurve2.on_click(shows2)
    savedata.on_click(data)
    savefig.on_click(fig)
    
    with output:
        display(widgets.VBox([loc, title1, file1, file2, file3, file4,\
                              file5, file6, file7, file8, file9, file10, datatype, showcurve]))
    
    with output4:
        display(widgets.VBox([loc, title2, csv1, csv2, csv3, csv4, csv5,\
                              csv6, showcurve2]))
    
    box_layout1 = widgets.Layout(display='flex', flew_flow='column',\
                                 overflow = 'scroll hidden', align_items='stretch',\
                                 border='solid', height='550px', width='26%')
    box_layout2 = widgets.Layout(display='flex', flew_flow='column',\
                                 overflow = 'scroll hidden', align_items='stretch',\
                                 border='solid', height='550px', width='37%')
    frame1 = widgets.Output()
    #frame1 = widgets.VBox([output1, output4])
    frame2 = widgets.VBox([output2, controloutput1])
    frame1.layout = box_layout1
    frame2.layout = box_layout2
    output3.layout = box_layout2

    l = widgets.jslink((dataname, 'value'), (figname, 'value'))
    
    accordion = widgets.Accordion(children = [output, output4],\
                                  titles = ('Create Result', 'Result Combination'),\
                                  layout = widgets.Layout(width='98.5%'))
    with frame1:
        display(accordion)

    display(widgets.HBox([frame1, frame2, output3]))