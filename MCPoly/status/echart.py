import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from status import status

def echart(files=[], loc='./', energy_pattern='Energy', hartree=False,\
           absolute=False, fig_pattern='line', figdata=True, savefig=True,\
           figname='Result', savedata=True, dataname='Result', xx=False):
    """
        A method to make a chart with calculated Electronic Energy datas or Gibbs Free Energy data.
        echart(files=[], loc='./', energy_pattern='Energy', hartree=False, absolute=False, figname='Result', fig_pattern='line', savefig=False, dataname='Result', savedata=True)
        files: The name of ORCA output files already calculated.
        loc: File Location. The default is your current location.
        energy_pattern: Energy you need. There are 2 options, 'Energy' for normal Electronic Energy and 'Gibbs' for Gibbs Free Energy. If you want to use 'Gibbs', please make sure your output file including frequency and thermodynamics calculation.
        hartree: Energy unit is normally kcal/mol, if you need unit 'Hartree', please input 'hartree=True'
                 TIPS: 1 Hartree ≈ 627.509 kcal/mol
        absolute: Normally we'll automatical use Relative Energy to show energy change, if you want Absolute Energy, please input 'absolute=True'
        fig_pattern: Pattern of the figure. There are 2 options, 'line' means Line Chart and 'bar' means Bar Chart.
        figdata: Show the File Name and Energy on the chart. The defaulr is True.
        savefig / savedata: To show if you want to save the chart and .txt result. The default is True.
        figname / savename: To set the name of the chart and .txt result. The default is 'Result'. 
        TIPS: If you use it to collect Gibbs Free Energy, please make sure you have frequency calculation. If the programme can't find out thermodynamics calculation, it will return None(NaN).
        
        Example:
            Input:
                from MCPoly.status import echart
                a = echart(files['Molecule1', 'Molecule2', 'Molecule3', 'Molecule4', 'Molecule5', 'Molecule6'], fig_pattern='line', \
                    figdata=True)
                print(a)
            Output:
                <<matplotlib.plot>>
                [ 0.          0.79379889  3.3396029   8.29064891 16.99670877 33.95325697]
            In Result.csv:
                Step, File, ∆E(kcal/mol)
                0, Molecule1, 0.0
                1, Molecule2, 0.7937988850646964
                2, Molecule3, 3.3396028980577155
                3, Molecule4, 8.290648908041941
                4, Molecule5, 16.996708774033824
                5, Molecule6, 33.95325697202637
        
        TIPS: For some limitation, the chart created on this function might not very good.
    """
    datas = []
    opath = os.getcwd()
    if files == []:
        for path in os.listdir(loc):
            if os.path.isfile(os.path.join(loc, path)):
                a = re.search('_Result.txt', path)
                if a:
                    files.append(path[:-11])
        files.sort()
    #print(files)
    if energy_pattern == 'Energy':
        for file in files:
            #try:
            #print(os.getcwd())
            data = status(file, loc).converged_energy()
            datas.append(data)
            #except:
            #    datas.append(None)
    elif energy_pattern == 'Gibbs':
        for file in files:
            #try:
            data = status(file, loc).gibbs()
            datas.append(data)
            #except:
            #    datas.append(None)
    os.chdir(opath)
    os.chdir(loc)
    Ebar = np.array(datas)
    if absolute == False:
        Ebar = Ebar - datas[0]
    if hartree == False:
        Ebar = Ebar * 627.509
    Ebar.tolist()
    Emax = -99999
    for E in Ebar:
        try:
            if E > Emax:
                Emax = E
        except:
            pass
    Emin = min(Ebar)
    gap = Emax - Emin
    if xx == False:
        if figdata == True:
            if fig_pattern == 'line':
                for i in range(len(Ebar)):
                    if files[i] == '':
                        continue
                    if Emax - Ebar[i] <= 5.0:
                        plt.text(i-0.5, Ebar[i] - 0.05*gap, '{0:.1f}'.format(Ebar[i]))
                    elif Emax - Ebar[i] > 5.0:
                        plt.text(i-0.2, Ebar[i] + 0.02*gap, '{0:.1f}'.format(Ebar[i]))
            elif fig_pattern == 'bar':
                for i in range(len(Ebar)):
                    plt.text(i, Ebar[i], '{0:.1f}'.format(Ebar[i]))
        if fig_pattern == 'line':
            ax = plt.plot(files, Ebar, '*-')
            if energy_pattern == 'Energy' and hartree == True:
                plt.ylabel('∆E(Eh)')
            elif energy_pattern == 'Energy' and hartree == False:
                plt.ylabel('∆E(kcal/mol)')
            elif energy_pattern == 'Gibbs' and hartree == True:
                plt.ylabel('∆G(Eh)')
            elif energy_pattern == 'Gibbs' and hartree == False:
                plt.ylabel('∆G(kcal/mol)')
        elif fig_pattern == 'bar':
            ax=plt.bar(files, Ebar)
            if energy_pattern == 'Energy' and hartree == True:
                plt.ylabel('∆E(Eh)')
            elif energy_pattern == 'Energy' and hartree == False:
                plt.ylabel('∆E(kcal/mol)')
            elif energy_pattern == 'Gibbs' and hartree == True:
                plt.ylabel('∆G(Eh)')
            elif energy_pattern == 'Gibbs' and hartree == False:
                plt.ylabel('∆G(kcal/mol)')
    if savefig == True:
        plt.savefig('{0}.png'.format(figname))
    if xx == False:
        plt.show()
    if savedata == True:
        try:
            f = open('{0}.csv'.format(dataname), 'x')
        except:
            f = open('{0}.csv'.format(dataname), 'w')
        f.write('Step,File,')
        if energy_pattern == 'Energy':
            f.write('∆E')
        elif energy_pattern == 'Gibbs':
            f.write('∆G')
        if hartree == True:
            f.write('(Eh)')
        else:
            f.write('(kcal/mol)')
        f.write('\n')
        for i in range(len(files)):
            f.write('{0},{1},{2}\n'.format(i, files[i], Ebar[i]))
        f.close()
    os.chdir(opath)
    return Ebar.tolist()