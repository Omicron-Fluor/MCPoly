import re
import os
import sys
import warnings
import math as m
from ase.io import read, write

mydir = os.path.dirname( __file__ )
lmpdir = os.path.join(mydir, '..', 'lmpset')
sys.path.append(lmpdir)
from chain import chain as normalchain
from infchain import infchain as inf
from DataToMolTxt import DataToMolTxt as ctxt
from DataToXyz import DataToXyz as cxyz
from rebuild import rebuild as rb

def replication(filename, x, xgap, y, ygap, z, zgap, absolute = False,\
                startloc = [0, 0, 0], loc = './', chain = False):
    if absolute == True:
        warnings.warn('Make sure you avoid all the overlaps of molecules!')
    opath = os.getcwd()
    os.chdir(loc)
    xloc = startloc[0]
    yloc = startloc[1]
    zloc = startloc[2]
    all_x = []
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    l = []
    real = [0, 0, 0]
    try:
        f  =  open('{0}_{1}{2}{3}.data'.format(filename, x, y, z), 'x')
    except:
        f  =  open('{0}_{1}{2}{3}.data'.format(filename, x, y, z), 'w')
    inp  =  open('{0}.data'.format(filename), 'r')
    for line in inp:
        a1 = re.search(r' atoms', line)
        a2 = re.search(r' bonds', line)
        a3 = re.search(r' angles', line)
        a4 = re.search(r' dihedrals', line)
        a5 = re.search(r' impropers', line)
        a6 = re.search(r' xlo xhi', line)
        a7 = re.search(r' ylo yhi', line)
        a8 = re.search(r' zlo zhi', line)
        c1 = re.search('Atoms', line)
        c2 = re.search('Bonds', line)
        c3 = re.search('Angles', line)
        c4 = re.search('Dihedrals', line)
        c5 = re.search('Impropers', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            f.write('{0} atoms\n'.format(eval(b1.group(0)) * x * y * z))
            num = b1.group(0)
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            f.write('{0} bonds\n'.format(eval(b2.group(0)) * x * y * z))
            num2 = b2.group(0)
        elif a3:
            b3 = re.search(r'[0-9]+', line)
            f.write('{0} angles\n'.format(eval(b3.group(0)) * x * y * z))
            num3 = b3.group(0)
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            f.write('{0} dihedrals\n'.format(eval(b4.group(0)) * x * y * z))
            num4 = b4.group(0)
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            f.write('{0} impropers\n'.format(eval(b5.group(0)) * x * y * z))
            num5 = b5.group(0)
        elif a6:
            x1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            xd = eval(x1[1])-eval(x1[0])
            if absolute == True:
                real[0] = xd-1/x*xgap
            f.write('{0:.5f} {1:.5f} xlo xhi\n'.format(xloc,\
                                                       xloc + (xd-real[0])*x + xgap*(x-1)))
            if absolute == True:
                real[0] = xd
        elif a7:
            y1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            yd = eval(y1[1])-eval(y1[0])
            if absolute == True:
                real[1] = yd-1/y*ygap
            f.write('{0:.5f} {1:.5f} ylo yhi\n'.format(yloc,\
                                                       yloc + (yd-real[1])*y + ygap*(y-1)))
            if absolute == True:
                real[1] = yd
        elif a8:
            z1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            zd = eval(z1[1])-eval(z1[0])
            if absolute == True:
                real[2] = zd-1/z*zgap
            f.write('{0:.5f} {1:.5f} zlo zhi\n'.format(zloc,\
                                                       zloc + (zd-real[2])*z + zgap*(z-1)))
            if absolute == True:
                real[2] = zd
        elif c1:
            w1 = 1
            f.write(line)
            f.write('\n')
        elif c2:
            w2 = 1
            f.write(line)
            f.write('\n')
        elif c3:
            w3 = 1
            f.write(line)
            f.write('\n')
        elif c4:
            w4 = 1
            f.write(line)
            f.write('\n')
        elif c5:
            w5 = 1
            f.write(line)
            f.write('\n')
        elif w1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(l1) >= 6:
                l.append(l1)
                if l1[0] == num:
                    w1 = 2
                    j = 1
        elif w1 == 2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,\
                                                                                           eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),\
                                                                                           eval(l[i4][-3])+(xd-real[0]+xgap)*i1-eval(x1[0]),\
                                                                                           eval(l[i4][-2])+(yd-real[1]+ygap)*i2-eval(y1[0]),\
                                                                                           eval(l[i4][-1])+(zd-real[2]+zgap)*i3-eval(z1[0])))
                            if i1 + i2 + i3 == 0:
                                all_x.append(eval(l[i4][-3]))
                            j = j + 1
            l = []
            j = 1
            w1 = 0
            f.write('\n')
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if l1[0] == num2:
                    w2 = 2
                    j = 1
        elif w2 == 2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            if chain == True:
                                if i4 == len(l) - 1:
                                    if all_x[eval(l[i4][2])-1] <= all_x[eval(l[i4][3])-1]:
                                        if i1 == 0:
                                            if eval(l[i4][2]) > eval(l[i4][3])+eval(num)*k:
                                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                                   eval(l[i4][2]) + eval(num)*k + eval(num)*(x-1)*y*z),\
                                                        eval(l[i4][3]) + eval(num)*k)
                                            else:
                                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                                   eval(l[i4][2]) + eval(num)*k,\
                                                                                   eval(l[i4][3]) + eval(num)*k + eval(num)*(x-1)*y*z))
                                        else:
                                            if eval(l[i4][2]) > eval(l[i4][3]) + eval(num)*k:
                                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                                   eval(l[i4][2]) + eval(num)*k - eval(num)*y*z,\
                                                                                   eval(l[i4][3]) + eval(num)*k))
                                            else:
                                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                                   eval(l[i4][2]) + eval(num)*k,
                                                                                   eval(l[i4][3]) + eval(num)*k - eval(num)*y*z))
                                    elif all_x[eval(l[i4][2])-1] > all_x[eval(l[i4][3])-1]:
                                        if i1 == 0:
                                            if eval(l[i4][2]) > eval(l[i4][3]) + eval(num)*k:
                                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                                   eval(l[i4][2]) + eval(num)*k,\
                                                                                   eval(l[i4][3]) + eval(num)*k + eval(num)*(x-1)*y*z))
                                            else:
                                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                                   eval(l[i4][2]) + eval(num)*k + eval(num)*(x-1)*y*z,\
                                                                                   eval(l[i4][3]) + eval(num)*k))
                                        else:
                                            if eval(l[i4][2]) > eval(l[i4][3]) + eval(num)*k:
                                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                                   eval(l[i4][2]) + eval(num)*k,\
                                                                                   eval(l[i4][3]) + eval(num)*k - eval(num)*y*z))
                                            else:
                                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                                   eval(l[i4][2]) + eval(num)*k - eval(num)*y*z,\
                                                                                   eval(l[i4][3]) + eval(num)*k))
                                else:
                                    f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                       eval(l[i4][2]) + eval(num)*k,\
                                                                       eval(l[i4][3]) + eval(num)*k))
                            else:
                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                   eval(l[i4][2]) + eval(num)*k,\
                                                                   eval(l[i4][3]) + eval(num)*k))
                            j = j + 1
                        k = k + 1
            l = []
            j = 1
            w2 = 0
            k = 0
            f.write('\n')
        elif w3 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 5:
                l.append(l1)
                if l1[0] == num3:
                    w3 = 2
                    j = 1
        elif w3 == 2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            f.write('{0} {1} {2} {3} {4}\n'.format(j, eval(l[i4][1]),\
                                                                   eval(l[i4][2]) + eval(num)*k,\
                                                                   eval(l[i4][3]) + eval(num)*k,\
                                                                   eval(l[i4][4]) + eval(num)*k))
                            j = j + 1
                        k = k + 1
            l = []
            j = 1
            k = 0
            w3 = 0
            f.write('\n')
        elif w4 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num4:
                    w4 = 2
                    j = 1
        elif w4 == 2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                                       eval(l[i4][2]) + eval(num)*k,\
                                                                       eval(l[i4][3]) + eval(num)*k,\
                                                                       eval(l[i4][4]) + eval(num)*k,\
                                                                       eval(l[i4][5]) + eval(num)*k))
                            j = j + 1
                        k = k + 1
            l = []
            j = 1
            k = 0
            w4 = 0
            f.write('\n')
        elif w5 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num5:
                    w5 = 2
                    j = 1
        elif w5 == 2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                                       eval(l[i4][2]) + eval(num)*k,\
                                                                       eval(l[i4][3]) + eval(num)*k,\
                                                                       eval(l[i4][4]) + eval(num)*k,\
                                                                       eval(l[i4][5])+eval(num)*k))
                            j = j + 1
                        k = k + 1
            l = []
            j = 1
            k = 0
            w4 = 0
            f.write('\n')
        else:
            f.write(line)
    f.close()
    os.chdir(opath)

def replication_brick(filename, x, xgap, y, ygap, z, zgap,\
                      xpattern = '0', ypattern = '0', zpattern = '0',\
                      shuffle = 0, absolute = False, startloc = [0, 0, 0], loc = './'):
    if xpattern == '0' and ypattern == '0' and zpattern == '0':
        warnings.warn("You didn't set any patterns!")
    xloc = startloc[0]
    yloc = startloc[1]
    zloc = startloc[2]
    opath = os.getcwd()
    os.chdir(loc)
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    half = 0.0
    l = []
    real = [0, 0, 0]
    pattern = ''
    if xpattern == 'y':
        pattern = 'xy'
    if xpattern == 'z':
        pattern = 'xz'
    if ypattern == 'x':
        pattern = 'yx'
    if ypattern == 'z':
        pattern = 'yz'
    if zpattern == 'x':
        pattern = 'zx'
    if zpattern == 'y':
        pattern = 'zy'
    try:
        f  =  open('{0}_{1}{2}{3}_brick{4}.data'.format(filename,\
                                                        x, y, z, pattern), 'x')
    except:
        f  =  open('{0}_{1}{2}{3}_brick{4}.data'.format(filename,\
                                                        x, y, z, pattern), 'w')
    inp  =  open('{0}.data'.format(filename), 'r')
    for line in inp:
        i1, i2, i3 = None, None, None
        a1 = re.search(r' atoms', line)
        a2 = re.search(r' bonds', line)
        a3 = re.search(r' angles', line)
        a4 = re.search(r' dihedrals', line)
        a5 = re.search(r' impropers', line)
        a6 = re.search(r' xlo xhi', line)
        a7 = re.search(r' ylo yhi', line)
        a8 = re.search(r' zlo zhi', line)
        c1 = re.search('Atoms', line)
        c2 = re.search('Bonds', line)
        c3 = re.search('Angles', line)
        c4 = re.search('Dihedrals', line)
        c5 = re.search('Impropers', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            if xpattern != '0':
                if xpattern == 'y':
                    axpattern = 'z'
                else:
                    axpattern = 'y'
                f.write('{0:.0f} atoms\n'.format(eval(b1.group(0))*(x*y*z
                                                                    - int(x/2) * eval(axpattern)
                                                 - int(eval(axpattern)/2) * (eval(axpattern)%2))))
            if ypattern != '0':
                if ypattern == 'z':
                    aypattern = 'x'
                else:
                    aypattern = 'z'
                f.write('{0:.0f} atoms\n'.format(eval(b1.group(0)) * (x*y*z
                                                                      - int(y/2) * eval(aypattern)
                                                 - int(eval(aypattern)/2) * (eval(aypattern)%2))))
            if zpattern != '0':
                if zpattern == 'y':
                    azpattern = 'x'
                else:
                    azpattern = 'y'
                f.write('{0:.0f} atoms\n'.format(eval(b1.group(0)) * (x*y*z
                                                                      - int(z/2) * eval(azpattern)
                                                 - int(eval(azpattern)/2) * (eval(azpattern)%2))))
            num = b1.group(0)
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            if xpattern != '0':
                f.write('{0:.0f} bonds\n'.format(eval(b2.group(0)) * (x*y*z
                                                                      - int(x/2) * eval(axpattern)
                                                                      - int(eval(axpattern)/2) * (eval(axpattern)%2))))
            if ypattern != '0':
                f.write('{0:.0f} bonds\n'.format(eval(b2.group(0)) * (x*y*z
                                                                      - int(y/2) * eval(aypattern)
                                                                      - int(eval(aypattern)/2) * (eval(aypattern)%2))))
            if zpattern != '0':
                f.write('{0:.0f} bonds\n'.format(eval(b2.group(0)) * (x*y*z
                                                                      - int(z/2) * eval(azpattern)
                                                                      - int(eval(azpattern)/2) * (eval(azpattern)%2))))
            num2 = b2.group(0)
        elif a3:
            b3 = re.search(r'[0-9]+', line)
            if xpattern != '0':
                f.write('{0:.0f} angles\n'.format(eval(b3.group(0)) * (x*y*z
                                                                       - int(x/2) * eval(axpattern)
                                                                       - int(eval(axpattern)/2) * (eval(axpattern)%2))))
            if ypattern != '0':
                f.write('{0:.0f} angles\n'.format(eval(b3.group(0)) * (x*y*z
                                                                       - int(y/2) * eval(aypattern)
                                                                       - int(eval(aypattern)/2) * (eval(aypattern)%2))))
            if zpattern != '0':
                f.write('{0:.0f} angles\n'.format(eval(b3.group(0)) * (x*y*z
                                                                       - int(z/2) * eval(azpattern)
                                                                       - int(eval(azpattern)/2) * (eval(azpattern)%2))))
            num3 = b3.group(0)
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            if xpattern != '0':
                f.write('{0:.0f} dihedrals\n'.format(eval(b4.group(0)) * (x*y*z
                                                                          - int(x/2) * eval(axpattern)
                                                                          - int(eval(axpattern)/2) * (eval(axpattern)%2))))
            if ypattern != '0':
                f.write('{0:.0f} dihedrals\n'.format(eval(b4.group(0)) * (x*y*z
                                                                          - int(y/2) * eval(aypattern)
                                                                          - int(eval(aypattern)/2) * (eval(aypattern)%2))))
            if zpattern != '0':
                f.write('{0:.0f} dihedrals\n'.format(eval(b4.group(0)) * (x*y*z
                                                                          - int(z/2) * eval(azpattern)
                                                                          - int(eval(azpattern)/2) * (eval(azpattern)%2))))
            num4 = b4.group(0)
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            if xpattern != '0':
                f.write('{0:.0f} impropers\n'.format(eval(b5.group(0)) * (x*y*z
                                                                          - int(x/2) * eval(axpattern)
                                                                          - int(eval(axpattern)/2) * (eval(axpattern)%2))))
            if ypattern != '0':
                f.write('{0:.0f} impropers\n'.format(eval(b5.group(0)) * (x*y*z
                                                                          - int(y/2) * eval(aypattern)
                                                                          - int(eval(aypattern)/2) * (eval(aypattern)%2))))
            if zpattern != '0':
                f.write('{0:.0f} impropers\n'.format(eval(b5.group(0)) * (x*y*z
                                                                          - int(z/2) * eval(azpattern)
                                                                          - int(eval(azpattern)/2) * (eval(azpattern)%2))))
            num5 = b5.group(0)
        elif a6:
            x1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            xd = eval(x1[1])-eval(x1[0])
            if absolute == True:
                real[0] = xd-1/x*xgap
            f.write('{0:.5f} {1:.5f} xlo xhi\n'.format(xloc,\
                                                       xloc + (xd-real[0])*x + xgap*(x-1)))
            if absolute == True:
                real[0] = xd
        elif a7:
            y1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            yd = eval(y1[1])-eval(y1[0])
            if absolute == True:
                real[1] = yd-1/y*ygap
            f.write('{0:.5f} {1:.5f} ylo yhi\n'.format(yloc,\
                                                       yloc + (yd-real[1])*y + ygap*(y-1)))
            if absolute == True:
                real[1] = yd
        elif a8:
            z1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            zd = eval(z1[1])-eval(z1[0])
            if absolute == True:
                real[2] = zd-1/z*zgap
            f.write('{0:.5f} {1:.5f} zlo zhi\n'.format(zloc,\
                                                       zloc + (zd-real[2])*z + zgap*(z-1)))
            if absolute == True:
                real[2] = zd
        elif c1:
            w1 = 1
            f.write(line)
            f.write('\n')
        elif c2:
            w2 = 1
            f.write(line)
            f.write('\n')
        elif c3:
            w3 = 1
            f.write(line)
            f.write('\n')
        elif c4:
            w4 = 1
            f.write(line)
            f.write('\n')
        elif c5:
            w5 = 1
            f.write(line)
            f.write('\n')
        elif w1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(l1) >= 6:
                l.append(l1)
                if l1[0] == num:
                    w1 = 2
                    j = 1
        elif w1 == 2:
            if xpattern != '0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern == 'y':
                            if i3 == None:
                                i3 = 1
                            half = i1 % 2 / 2
                        for i3 in range(z):
                            if xpattern == 'y':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5 - half
                            if xpattern == 'z':
                                half = i1 % 2 / 2
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5 - half
                            for i4 in range(len(l)):
                                if xpattern == 'y':
                                    if half == 0.5 and i2 == y - 1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,\
                                                                                                       eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),\
                                                                                                       eval(l[i4][-3]) + (xd+xgap-real[0])*i1 - eval(x1[0]),\
                                                                                                       eval(l[i4][-2]) + (yd+ygap-real[1])*(i2+half) - eval(y1[0]),\
                                                                                                       eval(l[i4][-1]) + (zd+zgap-real[2])*i3 - eval(z1[0])))
                                        j = j + 1
                                elif xpattern == 'z':
                                    if half == 0.5 and i3 == z - 1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,\
                                                                                                       eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),\
                                                                                                       eval(l[i4][-3]) + (xd+xgap-real[0])*i1 - eval(x1[0]),\
                                                                                                       eval(l[i4][-2]) + (yd+ygap-real[1])*i2 - eval(y1[0]),\
                                                                                                       eval(l[i4][-1]) + (zd+zgap-real[2])*(i3+half) - eval(z1[0])))
                                        j = j + 1
            elif ypattern != '0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern == 'x':
                            if i3 == None:
                                i3 = 1
                            half = i2 % 2 / 2
                        for i3 in range(z):
                            if ypattern == 'x':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5 - half
                            if ypattern == 'z':
                                half = i2 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            for i4 in range(len(l)):
                                if ypattern == 'x':
                                    if half == 0.5 and i1 == x - 1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,\
                                                                                                       eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),\
                                                                                                       eval(l[i4][-3]) + (xd+xgap-real[0])*(i1+half) - eval(x1[0]),\
                                                                                                       eval(l[i4][-2]) + (yd+ygap-real[1])*i2 - eval(y1[0]),\
                                                                                                       eval(l[i4][-1]) + (zd+zgap-real[2])*i3 - eval(z1[0])))
                                        j = j + 1
                                elif ypattern == 'z':
                                    if half == 0.5 and i3 == z - 1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,\
                                                                                                       eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),\
                                                                                                       eval(l[i4][-3]) + (xd+xgap-real[0])*i1 - eval(x1[0]),\
                                                                                                       eval(l[i4][-2]) + (yd+ygap-real[1])*i2 - eval(y1[0]),\
                                                                                                       eval(l[i4][-1]) + (zd+zgap-real[2])*(i3+half) - eval(z1[0])))
                                        j = j + 1
            elif zpattern != '0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern == 'x':
                            if i2 == None:
                                i2 = 1
                            half = i3 % 2 / 2
                        for i2 in range(y):
                            if zpattern == 'x':
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5 - half
                            if zpattern == 'y':
                                half = i3 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            for i4 in range(len(l)):
                                if zpattern == 'x':
                                    if half == 0.5 and i1 == x - 1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,\
                                                                                                       eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),\
                                                                                                       eval(l[i4][-3]) + (xd-real[0]+xgap)*(i1+half) - eval(x1[0]),\
                                                                                                       eval(l[i4][-2]) + (yd-real[1]+ygap)*i2 - eval(y1[0]),\
                                                                                                       eval(l[i4][-1]) + (zd-real[2]+zgap)*i3 - eval(z1[0])))
                                        j = j + 1
                                elif zpattern == 'y':
                                    if half == 0.5 and i2 == y - 1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,\
                                                                                                       eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),\
                                                                                                       eval(l[i4][-3]) + (xd-real[0]+xgap)*i1 - eval(x1[0]),\
                                                                                                       eval(l[i4][-2]) + (yd-real[1]+ygap)*(i2+half) - eval(y1[0]),\
                                                                                                       eval(l[i4][-1]) + (zd-real[2]+zgap)*i3 - eval(z1[0])))
                                        j = j + 1
            l = []
            j = 1
            w1 = 0
            f.write('\n')
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if l1[0] == num2:
                    w2 = 2
                    j = 1
        elif w2 == 2:
            if xpattern != '0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern == 'y':
                            if i3 == None:
                                i3 = 1
                            half = i1%2/2
                        for i3 in range(z):
                            if xpattern == 'y':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5-half
                            if xpattern == 'z':
                                half = i1%2/2
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5-half
                            if half == 0.5 and i2 == y - 1 and xpattern == 'y':
                                continue
                            if half == 0.5 and i3 == z - 1 and xpattern == 'z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                   eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            if ypattern != '0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern == 'x':
                            if i3 == None:
                                i3 = 1
                            half = i2 % 2 / 2
                        for i3 in range(z):
                            if ypattern == 'x':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5 - half
                            if ypattern == 'z':
                                half = i2 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i1 == x - 1 and ypattern == 'x':
                                continue
                            if half == 0.5 and i3 == z - 1 and ypattern == 'z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                   eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            if zpattern != '0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern == 'x':
                            if i2 == None:
                                i2 = 1
                            half = i3 % 2 / 2
                        for i2 in range(y):
                            if zpattern == 'x':
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5 - half
                            if zpattern == 'y':
                                half = i3 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i1 == x - 1 and zpattern == 'x':
                                continue
                            if half == 0.5 and i2 == y - 1 and zpattern == 'y':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                   eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            l = []
            j = 1
            w2 = 0
            k = 0
            f.write('\n')
        elif w3 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 5:
                l.append(l1)
                if l1[0] == num3:
                    w3 = 2
                    j = 1
        elif w3 == 2:
            if xpattern != '0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern == 'y':
                            if i3 == None:
                                i3 = 1
                            half = i1 % 2 / 2
                        for i3 in range(z):
                            if xpattern == 'y':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5 - half
                            if xpattern == 'z':
                                half = i1 % 2 / 2
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i2 == y - 1 and xpattern == 'y':
                                continue
                            if half == 0.5 and i3 == z - 1 and xpattern == 'z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4}\n'.format(j, eval(l[i4][1]),\
                                                                       eval(l[i4][2]) + eval(num)*k,\
                                                                       eval(l[i4][3]) + eval(num)*k,\
                                                                       eval(l[i4][4]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            if ypattern != '0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern == 'x':
                            if i3 == None:
                                i3 = 1
                            half = i2 % 2 / 2
                        for i3 in range(z):
                            if ypattern == 'x':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5 - half
                            if ypattern == 'z':
                                half = i2 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i1 == x - 1 and ypattern == 'x':
                                continue
                            if half == 0.5 and i3 == z - 1 and ypattern == 'z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4}\n'.format(j, eval(l[i4][1]),\
                                                                       eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k,\
                                                                       eval(l[i4][4]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            if zpattern != '0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern == 'x':
                            if i2 == None:
                                i2 = 1
                            half = i3 % 2 / 2
                        for i2 in range(y):
                            if zpattern == 'x':
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5 - half
                            if zpattern == 'y':
                                half = i3 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i1 == x - 1 and zpattern == 'x':
                                continue
                            if half == 0.5 and i2 == y - 1 and zpattern == 'y':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4}\n'.format(j, eval(l[i4][1]),\
                                                                       eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k,\
                                                                       eval(l[i4][4]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            l = []
            j = 1
            k = 0
            w3 = 0
            f.write('\n')
        elif w4 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num4:
                    w4 = 2
                    j = 1
        elif w4 == 2:
            if xpattern != '0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern == 'y':
                            if i3 == None:
                                i3 = 1
                            half = i1 % 2 / 2
                        for i3 in range(z):
                            if xpattern == 'y':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5 - half
                            if xpattern == 'z':
                                half = i1 % 2 / 2
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i2 == y - 1 and xpattern == 'y':
                                continue
                            if half == 0.5 and i3 == z - 1 and xpattern == 'z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                                           eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k,\
                                                                           eval(l[i4][4]) + eval(num)*k, eval(l[i4][5]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            if ypattern != '0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern == 'x':
                            if i3 == None:
                                i3 = 1
                            half = i2 % 2 / 2
                        for i3 in range(z):
                            if ypattern == 'x':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5 - half
                            if ypattern == 'z':
                                half = i2 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i1 == x - 1 and ypattern == 'x':
                                continue
                            if half == 0.5 and i3 == z - 1 and ypattern == 'z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                                           eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k,\
                                                                           eval(l[i4][4]) + eval(num)*k, eval(l[i4][5]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            if zpattern != '0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern == 'x':
                            if i2 == None:
                                i2 = 1
                            half = i3 % 2 / 2
                        for i2 in range(y):
                            if zpattern == 'x':
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5 - half
                            if zpattern == 'y':
                                half = i3 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i1 == x - 1 and zpattern == 'x':
                                continue
                            if half == 0.5 and i2 == y - 1 and zpattern == 'y':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                                           eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k,\
                                                                           eval(l[i4][4]) + eval(num)*k, eval(l[i4][5]) + eval(num)*k))
                                j = j + 1
                            k = k + 1

            l = []
            j = 1
            k = 0
            w4 = 0
            f.write('\n')
        elif w5 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num5:
                    w5 = 2
                    j = 1
        elif w5 == 2:
            if xpattern != '0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern == 'y':
                            if i3 == None:
                                i3 = 1
                            half = i1 % 2 / 2
                        for i3 in range(z):
                            if xpattern == 'y':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5 - half
                            if xpattern == 'z':
                                half = i1%2/2
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i2 == y - 1 and xpattern == 'y':
                                continue
                            if half == 0.5 and i3 == z - 1 and xpattern == 'z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                                           eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k,\
                                                                           eval(l[i4][4]) + eval(num)*k, eval(l[i4][5]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            if ypattern != '0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern == 'x':
                            if i3 == None:
                                i3 = 1
                            half = i2 % 2 / 2
                        for i3 in range(z):
                            if ypattern == 'x':
                                if shuffle == 1 and i3 % 2 == 1:
                                    half = 0.5 - half
                            if ypattern == 'z':
                                half = i2 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i1 == x - 1 and ypattern == 'x':
                                continue
                            if half == 0.5 and i3 == z - 1 and ypattern == 'z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                                           eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k,\
                                                                           eval(l[i4][4]) + eval(num)*k, eval(l[i4][5]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            if zpattern != '0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern == 'x':
                            if i2 == None:
                                i2 = 1
                            half = i3 % 2 / 2
                        for i2 in range(y):
                            if zpattern == 'x':
                                if shuffle == 1 and i2 % 2 == 1:
                                    half = 0.5 - half
                            if zpattern == 'y':
                                half = i3 % 2 / 2
                                if shuffle == 1 and i1 % 2 == 1:
                                    half = 0.5 - half
                            if half == 0.5 and i1 == x - 1 and zpattern == 'x':
                                continue
                            if half == 0.5 and i2 == y - 1 and zpattern == 'y':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                                           eval(l[i4][2]) + eval(num)*k, eval(l[i4][3]) + eval(num)*k,\
                                                                           eval(l[i4][4]) + eval(num)*k, eval(l[i4][5]) + eval(num)*k))
                                j = j + 1
                            k = k + 1
            l = []
            j = 1
            k = 0
            w4 = 0
            f.write('\n')
        else:
            f.write(line)
    f.close()
    os.chdir(opath)

class mould:
    """
    The definition of molecule system.
    mould(polymer, loc = './')
    atoms: Your molecule system name on your .data file.
    loc: File Location. The default is your current location.
    You can get the further information by .cube and .brick.
    """
    
    def __init__(self, atoms, loc = './'):
        self.atoms = atoms
        self.loc = loc
    
    def cube(self, x, xgap, y, ygap, z, zgap,\
             absolute = False, startloc = [0, 0, 0], chain = False):
        """
    The method to create a box with replicated molecular, or create a crystal
    * * * * *
    * * * * *
    * * * * *
    * * * * *     Cube Structure
    
    cube(x, xgap, y, ygap, z, zgap, absolute = False, startloc = [0, 0, 0])
    x: number axis in x axis
    xgap: distance of nearest atoms between each molecule systems in x axis direction
    y: number axis in y axis
    ygap: distance of nearest atoms between each molecule systems in y axis direction
    z: number axis in z axis
    zgap: distance of nearest atoms between each molecule systems in z axis direction
    absolute: Use the distance between the molecule systems box instead of nearest atoms in xgap, ygap and zgap. The default is True.
    startloc: The location of the beginning the replicated system drawing. Mostly, it is minimum of x, y, z in molecule system.
    chain: To connect chain by this function in x axis. The default is false.
        TIPS: It can only alleviate the problem of connecting multiple polymer molecules, but to totally solve this problem, the best way is to get a single long polymer chain.
    Example:
        Input:
            from MCPoly.lmpset import mould
            atoms = mould('Poly1')
            atoms.cube(6, 5, 6, 5, 3, 5)
        Output in Poly1_663.data AND Poly1.data (original file):
            #File: Poly1                                    #File: Poly1

            972 atoms                                       9 atoms
            864 bonds                                       8 bonds
            1404 angles                                     13 angles
            1296 dihedrals                                  12 dihedrals
            432 impropers                                   4 impropers

            4 atom types                                    4 atom types
            4 bond types                                    4 bond types
            5 angle types                                   5 angle types
            4 dihedral types                                4 dihedral types
            1 improper types                                1 improper types
            
            0.00000 53.32948 xlo xhi  -> startloc[0]        -9.29324 -4.57166 xlo xhi
            0.00000 54.87382 ylo yhi  -> startloc[1]        -1.72888 3.25009 ylo yhi
            0.00000 26.46613 zlo zhi  -> startloc[2]        -2.04102 3.44769 zlo zhi
            ...                                             ...
            
            
            Atoms                                           Atoms

            1 1 1 -0.2893 1.62654 2.03208 1.98722           1 1 1 -0.28935 -7.66305 0.30313 -0.05343  
            2 1 1 0.1096 1.09404 2.95284 3.06549            2 1 1 0.10963 -8.19932 1.22343 1.02443
            3 1 2 0.0915 1.27330 2.34271 1.00000            3 1 2 0.0915 -8.01953 0.61332 -1.04153
            4 1 2 0.0915 2.72158 2.03208 1.98722            4 1 2 0.0915 -6.57166 0.30320 -0.05380
            5 1 3 -0.6883 1.56973 2.53667 4.33729           5 1 3 -0.6883 -7.72351 0.80779 2.29626
            6 1 2 0.0939 0.00000 2.94776 3.08107            6 1 2 0.0939 -9.29323 1.21888 1.04005
            
            ...                                             ...
        """
        return replication(self.atoms, x, xgap, y, ygap, z,\
                           zgap, absolute, startloc, self.loc, chain)
    
    def brick(self, x, xgap, y, ygap, z, zgap,\
              xpattern = '0', ypattern = '0', zpattern = '0',\
              shuffle = 0, absolute = False, startloc = [0, 0, 0]):
        """
    The method to create a box with replicated molecular in brick settlement.
    * * * * *
     * * * *
    * * * * *
     * * * *       Brick Structure
    
    brick(x, xgap, y, ygap, z, zgap, xpattern = '0', ypattern = '0', zpattern = '0', shuffle = 0, absolute = False, startloc = [0, 0, 0])
    x: number axis in x axis
    xgap: distance of nearest atoms between each molecule systems in x axis direction
    y: number axis in y axis
    ygap: distance of nearest atoms between each molecule systems in y axis direction
    z: number axis in z axis
    zgap: distance of nearest atoms between each molecule systems in z axis direction
    (x, y, z)pattern: To show the direction of you brick patterns. e.g. xpattern = 'y' mean brick patterns on xy plane and other horizental planes, with linear aligns in y direction. 
    shuffle: With two reversed brick pattern on alternate planes. The default is false.
    absolute: Use the distance between the molecule systems box instead of nearest atoms in xgap, ygap and zgap. The default is True.
    startloc: The location of the beginning the replicated system drawing. Mostly, it has minimum of x, y, z.
    TIPS: You can only use  one of (x, y, z)pattern.
    Example:
        Input:
            from MCPoly.lmpset import mould
            atoms = mould('Poly1')
            atoms.brick(6, 5, 6, 5, 4, 5, xpattern = 'y')
        Output in Poly1_664_brickxy.data AND Poly1.data (original file):
            #File: Poly1                                    #File: Poly1

            1188 atoms                                      9 atoms
            1056 bonds                                      8 bonds
            1716 angles                                     13 angles
            1584 dihedrals                                  12 dihedrals
            528 impropers                                   4 impropers

            4 atom types                                    4 atom types
            4 bond types                                    4 bond types
            5 angle types                                   5 angle types
            4 dihedral types                                4 dihedral types
            1 improper types                                1 improper types
            
            0.00000 53.32948 xlo xhi  -> startloc[0]        -9.29324 -4.57166 xlo xhi
            0.00000 54.87382 ylo yhi  -> startloc[1]        -1.72888 3.25009 ylo yhi
            0.00000 36.95484 zlo zhi  -> startloc[2]        -2.04102 3.44769 zlo zhi
            ...                                             ...
            
            Atoms                                           Atoms

            1 1 1 -0.2893 1.62654 2.03208 1.98722           1 1 1 -0.28935 -7.66305 0.30313 -0.05343  
            2 1 1 0.1096 1.09404 2.95284 3.06549            2 1 1 0.10963 -8.19932 1.22343 1.02443
            3 1 2 0.0915 1.27330 2.34271 1.00000            3 1 2 0.0915 -8.01953 0.61332 -1.04153
            4 1 2 0.0915 2.72158 2.03208 1.98722            4 1 2 0.0915 -6.57166 0.30320 -0.05380
            5 1 3 -0.6883 1.56973 2.53667 4.33729           5 1 3 -0.6883 -7.72351 0.80779 2.29626
            6 1 2 0.0939 0.00000 2.94776 3.08107            6 1 2 0.0939 -9.29323 1.21888 1.04005
            
            ...                                             ...
        """
        
        if xpattern != '0' and ypattern != '0':
            raise AssertionError('You can only use  one of (x, y, z)pattern.')
        if xpattern != '0' and zpattern != '0':
            raise AssertionError('You can only use  one of (x, y, z)pattern.')
        if ypattern != '0' and zpattern != '0':
            raise AssertionError('You can only use  one of (x, y, z)pattern.')
        return replication_brick(self.atoms, x, xgap, y, ygap, z, zgap,\
                                 xpattern, ypattern, zpattern, shuffle, absolute, startloc, self.loc)

    def DataToXyz(self, savename = ''):
        """
        The method to change LAMMPS Data File into XYZ File. Powered by ASE.
        DataToXyz(savename = '')
        savename: Name of the saved XYZ File. The default is name of origin LAMMPS Data File.
        Example:
            Input:
                from MCPoly.lmpset import mould
                atoms = mould('Poly1')
                atoms.DataToXyz()
            Output in Poly1_Chain.xyz:
                82
                Properties = species:S:1:pos:R:3 pbc = "F F F"
                H       34.26389000       7.81856000       4.72112000
                C       33.60940000       7.29160000       5.40582000
                H       33.81841000       6.25379000       5.63907000
                C       32.30303000       7.89744000       5.79799000
                H       31.98877000       7.45674000       6.75085000
                C       31.23940000       7.64473000       4.74277000
        
                ...
        """
        return cxyz(self.atoms, self.loc, savename)
    
    def DataToMolTxt(self, types = {0:0}, savename = ''):
        """
        The method to change LAMMPS Data File into LAMMPS Molecule Text File.
        DataToMolTxt(types = {0:0, ...}, savename = '')
        file: Your molecule system name on your .data file.
        types: Change the Atom Type number to suit the main LAMMPS input files.
        loc: File Location. The default is your current location.
        savename: Name of the saved XYZ File. The default is name of origin LAMMPS Data File.
        
        Example 1:
            Input:
                from MCPoly.lmpset import mould
                atoms = mould('Poly1')
                atoms.DataToMolTxt()
                
            Output in Poly1.txt:
                #Polymer: BioPolymer1

                84 atoms
                85 bonds
                158 angles
                217 dihedrals
                48 impropers
                
                Coords
                
                1 1.3822 0.027 0.8364
                    
                    ...
           
                Types
        
                1 1
                2 2
                3 1
                4 2
                5 2
                6 3
                7 2
                8 4
                9 5
                
                    ...
        
        Example 2:
            Input:
                from MCPoly.lmpset import mould
                atoms = mould('Poly1')
                atoms.DataToMolTxt(types = {1:11}, savename = 'Poly1_rami')
                
            Output in Poly1_rami:
                #Polymer: BioPolymer1
    
                84 atoms
                85 bonds
                158 angles
                217 dihedrals
                48 impropers
                
                Coords
                
                1 1.3822 0.027 0.8364
                    
                    ...
           
                Types
        
                1 11
                2 2
                3 11
                4 2
                5 2
                6 3
                7 2
                8 4
                9 5
                
                    ...
        """
        return ctxt(self.atoms, types, self.loc, savename)
    
    def infchain(self, bondtype):
        """
        The method to create a single molecule for periodical chain, specialised for polymers.
        infchain(bondtype)
        file: Your molecule system name on your .data file.
        bondtype: The bond type between the start and the end of the polymer. 
        loc: File Location. The default is your current location.
        Example:
            Input:
                from MCPoly.lmpset import mould
                atoms = mould('Poly1')
                atoms.infchain(11)
            Output in Poly1_Chain.data:
                82 atoms
                84 bonds
                154 angles
                211 dihedrals
                47 impropers
                
                9 atom types
                12 bond types
                21 angle types
                31 dihedral types
                3 improper types
                
                ...
                84 11 79 1
                
                ...
        """
        return inf(self.atoms, bondtype, self.loc)

    def chain(self, bondtype, degree):
        """
        The method to create a single molecule for a finite chain, specialised for polymers.
        chain(bondtype, degree)
        degree: Degree of polymerisation.
        bondtype: The bond type between the start and the end of the polymer. 
        Example:
            Input:
                from MCPoly.lmpset import mould
                atoms = mould('Poly1')
                atoms.chain('BioPolymer2', 11, 16)
            Output in Poly1_16x.data:
                1314 atoms
                1345 bonds
                2468 angles
                3382 dihedrals
                753 impropers

                9 atom types
                12 bond types
                21 angle types
                31 dihedral types
                3 improper types

                ...
                82 1 83 1

                ...
        """
        return normalchain(self.atoms, bondtype, degree, self.loc)
    
    def rebuild(self, substitute, x = [99999, 99999], y = [99999, 99999], z = [99999, 99999], subloc = './'):
        """
        The method to replace the geometry structure of LAMMPS Input File from an XYZ File.
        rebuild(substitute, x = [xmin, xmax], y = [ymin, ymax], z = [zmin, zmax], subloc = './')
        substitute: XYZ File used to replace.
        x, y, z: Caption of the box.
        subloc: XYZ File Location. The default is your current location.
        Example:
            Input:
                from MCPoly.lmpset import mould
                
                atoms = mould('Atoms1')
                atoms = atoms.cube(1, 0, 1, 0, 2, 0) # Create the file Atoms1_112.data, see in lmpset.cube.
                atoms.DataToXyz('Atoms1_112') # Create Atoms1_112.xyz, see in lmpset.DataToXyz
                # After changing the geometry structure by GaussView or Avogrado and save it as Atoms1_replace.xyz
                atoms.rebuild('Atoms1_replace')
                
            
            Atoms1_replace.xyz:
                106
                
                O          9.64263        4.67845        3.83021
                C          9.97863        4.59047        2.66070
                N         11.06345        5.23157        2.07395
                O          9.31893        3.82840        1.73881
                C         12.08548        6.09614        2.75528
                H         11.25704        4.90681        1.13587
                C          8.15829        2.96807        2.00836
                ......
            
            Output in Atoms1_112.data:
                ......
                
                Atoms
                
                1 1 1 -0.4804 9.64263 4.67845 3.83021
                2 1 2 0.7694 9.97863 4.59047 2.66070
                3 1 3 -1.1201 11.06345 5.23157 2.07395
                4 1 4 -0.3798 9.31893 3.82840 1.73881
                5 1 5 0.1753 12.08548 6.09614 2.75528
                6 1 6 0.5150 11.25704 4.90681 1.13587
                7 1 5 0.0718 8.15829 2.96807 2.00836
                ......
            
    """
        return rb(self.atoms, substitute, x = x, y = y, z = z,\
                  loc = self.loc, subloc = subloc)