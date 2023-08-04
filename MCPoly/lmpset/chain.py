import re
import os
import sys

mydir = os.path.dirname( __file__ )
lmpdir = os.path.join(mydir, '..', 'lmpset')
sys.path.append(lmpdir)
from infchain import infchain

def chain(file, bondtype, degree, loc='./'):
    """
    The method to create a single molecule for a finite chain,  specialised for polymers.
    chain(file, bondtype, degree, loc='./')
    file: Your molecule system name on your .data file.
    degree: Degree of polymerisation.
    bondtype: The bond type between the start and the end of the polymer. 
    loc: File Location. The default is your current location.
    Example:
        Input:
            from MCPoly.lmpset import chain
            chain('BioPolymer2', 11, 16)
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
    opath = os.getcwd()
    os.chdir(loc)
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    n1 = 0
    t2 = 0
    t3 = 0
    t4 = 0
    t5 = 0
    H = []
    lx = []
    lost2_u = []
    lost2_d = []
    lost3_u = []
    lost3_d = []
    lost4_u = []
    lost4_d = []
    lost5_u = []
    lost5_d = []
    Hloc = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    l = []
    real = [0, 0, 0]
    inp = open('{0}.data'.format(file), 'r')
    for line in inp:
        a1 = re.search(r' atoms', line)
        a2 = re.search(r' bonds', line)
        a3 = re.search(r' angles', line)
        a4 = re.search(r' dihedrals', line)
        a5 = re.search(r' impropers', line)
        a6 = re.search(r' xlo xhi', line)
        a7 = re.search(r' ylo yhi', line)
        a8 = re.search(r' zlo zhi', line)
        m1 = re.search(r'Masses', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            tipnum = eval(b1.group(0))
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            tipnum2 = eval(b2.group(0))
        elif a3:
            b3 = re.search(r'[0-9]+', line)
            tipnum3 = eval(b3.group(0))
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            tipnum4 = eval(b4.group(0))
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            tipnum5 = eval(b5.group(0))
        elif a6:
            b6 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            tipgapx = eval(b6[1]) - eval(b6[0])
            xlo = eval(b6[0])
            xhi = eval(b6[1])
        elif a7:
            b7 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            gapy = eval(b7[1]) - eval(b7[0])
            ylo = eval(b7[0])
        elif a8:
            b8 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            gapz = eval(b8[1]) - eval(b8[0])
            zlo = eval(b8[0])
        if n1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l1 == []:
                continue
            else:
                n1 = 2
        if n1 == 2:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l1 == []:
                n1 = 0
                continue
            if eval(l1[-1]) < 1.1:
                H.append(eval(l1[0]))
        if m1:
            n1 = 1
        c1 = re.search('Atoms', line)
        c2 = re.search('Bonds', line)
        c3 = re.search('Angles', line)
        c4 = re.search('Dihedrals', line)
        c5 = re.search('Impropers', line)
        if c1:
            w1 = 1
        elif c2:
            w2 = 1
        elif c3:
            w3 = 1
        elif c4:
            w4 = 1
        elif c5:
            w5 = 1
        elif w1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(l1) >= 6:
                l.append(l1)
                if eval(l1[0]) == tipnum:
                    w1 = 2
        elif w1 == 2:
            xmin = 99999
            xmax = -99999
            xmaxnum = 0
            xminnum = 0
            realxmin = 99999
            realxmax = -99999
            realxmaxnum = 0
            realxminnum = 0
            for i4 in range(len(l)):
                if eval(l[i4][2]) not in H:
                    continue
                if eval(l[i4][-3]) > xmax:
                    xmax = eval(l[i4][-3])
                    xmaxnum = eval(l[i4][0])
                    Hloc[1] = [eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])]
                if eval(l[i4][-3]) < xmin:
                    xmin = eval(l[i4][-3])
                    xminnum = eval(l[i4][0])
                    Hloc[0] = [eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])]
            for i4 in range(len(l)):
                if eval(l[i4][2]) in H:
                    continue
                if eval(l[i4][-3]) > realxmax:
                    realxmax = eval(l[i4][-3])
                    realxmaxnum = eval(l[i4][0])
                if eval(l[i4][-3]) < realxmin:
                    realxmin = eval(l[i4][-3])
                    realxminnum = eval(l[i4][0])
            w1 = 0
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if eval(l1[0]) == tipnum2:
                    w2 = 2
                    j = 1
        elif w2 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    if eval(l[i4][2]) == xminnum:
                        lost2_d = [eval(l[i4][-2]), eval(l[i4][-1])]
                    elif eval(l[i4][2]) == xmaxnum:
                        lost2_u = [eval(l[i4][-2]), eval(l[i4][-1])]
                    continue
                elif eval(l[i4][3]) in [xminnum, xmaxnum]:
                    if eval(l[i4][3]) == xminnum:
                        lost2_d = [eval(l[i4][-2]), eval(l[i4][-1])]
                    elif eval(l[i4][3]) == xmaxnum:
                        lost2_u = [eval(l[i4][-2]), eval(l[i4][-1])]
                    continue
                j = j + 1
            k = k + 1
            l = []
            t2 = j
            j = 1
            w2 = 0
            k = 0
        elif w3 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 5:
                l.append(l1)
                if eval(l1[0]) == tipnum3:
                    w3 = 2
                    j = 1
        elif w3 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    if eval(l[i4][2]) == xminnum:
                        lost3_d.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][2]) == xmaxnum:
                        lost3_u.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][3]) in [xminnum, xmaxnum]:
                    if eval(l[i4][3]) == xminnum:
                        lost3_d.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    if eval(l[i4][3]) == xmaxnum:
                        lost3_u.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][4]) in [xminnum, xmaxnum]:
                    if eval(l[i4][4]) == xminnum:
                        lost3_d.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    if eval(l[i4][4]) == xmaxnum:
                        lost3_u.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                j = j + 1
            k = k + 1
            l = []
            t3 = j
            j = 1
            k = 0
            w3 = 0
        elif w4 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if eval(l1[0]) == tipnum4:
                    w4 = 2
                    j = 1
        elif w4 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    if eval(l[i4][2]) == xminnum:
                        lost4_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][2]) == xmaxnum:
                        lost4_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][3]) in [xminnum, xmaxnum]:
                    if eval(l[i4][3]) == xminnum:
                        lost4_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][3]) == xmaxnum:
                        lost4_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][4]) in [xminnum, xmaxnum]:
                    if eval(l[i4][4]) == xminnum:
                        lost4_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][4]) == xmaxnum:
                        lost4_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][5]) in [xminnum, xmaxnum]:
                    if eval(l[i4][5]) == xminnum:
                        lost4_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][5]) == xmaxnum:
                        lost4_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                j = j + 1
            k = k + 1
            l = []
            t4 = j
            j = 1
            k = 0
            w4 = 0
        elif w5 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if eval(l1[0]) == tipnum5:
                    w5 = 2
                    j = 1
        elif w5 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    if eval(l[i4][2]) == xminnum:
                        lost5_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][2]) == xmaxnum:
                        lost5_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][3]) in [xminnum, xmaxnum]:
                    if eval(l[i4][3]) == xminnum:
                        lost5_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][3]) == xmaxnum:
                        lost5_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][4]) in [xminnum, xmaxnum]:
                    if eval(l[i4][4]) == xminnum:
                        lost5_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][4]) == xmaxnum:
                        lost5_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][5]) in [xminnum, xmaxnum]:
                    if eval(l[i4][5]) == xminnum:
                        lost5_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][5]) == xmaxnum:
                        lost5_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                j = j + 1
            k = k + 1
            l = []
            t5 = j
            j = 1
            k = 0
            w4 = 0
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    gap = 0
    l = []
    try:
        med = open('{0}_Chain.data'.format(file), 'r')
    except:
        infchain(file, bondtype, loc)
        med = open('{0}_Chain.data'.format(file), 'r')
    for line in med:
        a1 = re.search(r' atoms', line)
        a2 = re.search(r' bonds', line)
        a3 = re.search(r' angles', line)
        a4 = re.search(r' dihedrals', line)
        a5 = re.search(r' impropers', line)
        a6 = re.search(r' xlo xhi', line)
        
        m1 = re.search(r'Masses', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            num = eval(b1.group(0))
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            num2 = eval(b2.group(0))
        elif a3:
            b3 = re.search(r'[0-9]+', line)
            num3 = eval(b3.group(0))
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            num4 = eval(b4.group(0))
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            num5 = eval(b5.group(0))
        elif a6:
            b6 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            gapx = eval(b6[1]) - eval(b6[0])
        if n1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l1 == []:
                continue
            else:
                n1 = 2
        if n1 == 2:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l1 == []:
                n1 = 0
                continue
            if eval(l1[-1])<1.1:
                H.append(eval(l1[0]))
        if m1:
            n1 = 1
        c1 = re.search('Atoms', line)
        c2 = re.search('Bonds', line)
        c3 = re.search('Angles', line)
        c4 = re.search('Dihedrals', line)
        c5 = re.search('Impropers', line)
        if c1:
            w1 = 1
        elif c2:
            w2 = 1
        elif c3:
            w3 = 1
        elif c4:
            w4 = 1
        elif c5:
            w5 = 1
        elif w1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(l1) >= 6:
                l.append(l1)
                if eval(l1[0]) == tipnum:
                    w1 = 2
        elif w1 == 2:
            w1 = 0
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if eval(l1[0]) == tipnum2:
                    w2 = 2
                    j = 1
        elif w2 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                j = j + 1
            k = k + 1
            l = []
            t2 = j
            j = 1
            w2 = 0
            k = 0
        elif w3 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 5:
                l.append(l1)
                if eval(l1[0]) == tipnum3:
                    w3 = 2
                    j = 1
        elif w3 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                j = j + 1
            k = k + 1
            l = []
            t3 = j
            j = 1
            k = 0
            w3 = 0
        elif w4 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if eval(l1[0]) == tipnum4:
                    w4 = 2
                    j = 1
        elif w4 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][5]) in [xminnum, xmaxnum]:
                    continue
                j = j + 1
            k = k + 1
            l = []
            t4 = j
            j = 1
            k = 0
            w4 = 0
        elif w5 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if eval(l1[0]) == tipnum5:
                    w5 = 2
                    j = 1
        elif w5 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][5]) in [xminnum, xmaxnum]:
                    continue
                j = j + 1
            k = k + 1
            l = []
            t5 = j
            j = 1
            k = 0
            w4 = 0
    try:
        f = open('{0}_{1}x.data'.format(file, degree), 'x')
    except:
        f = open('{0}_{1}x.data'.format(file, degree), 'w')
    inp.close()
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    gap = 0
    l = []
    all_x = []
    inp = open('{0}.data'.format(file), 'r')
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
            f.write('{0} atoms\n'.format((tipnum-2)*degree + 2))
            #allnum = b1.group(0)
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            f.write('{0} bonds\n'.format(tipnum2*degree - (degree-1)))
            #allnum2 = b2.group(0)
        elif a3:
            f.write('{0} angles\n'.format(tipnum3 + num3*(degree-1)))
            #allnum3 = b3.group(0)
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            f.write('{0} dihedrals\n'.format(tipnum4 + num4*(degree-1)))
            #allnum4 = b4.group(0)
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            f.write('{0} impropers\n'.format(tipnum5 + num5*(degree-1)))
            #allnum5 = b5.group(0)
        elif a6:
            x1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            xd = eval(x1[1]) - eval(x1[0])
            f.write('{0:.5f} {1:.5f} xlo xhi\n'.format(0.0,\
                                                       tipgapx*(degree-1) + gapx))
        elif a7:
            y1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            yd = eval(y1[1]) - eval(y1[0])
            f.write('{0:.5f} {1:.5f} ylo yhi\n'.format(0.0, gapy))
        elif a8:
            z1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            zd = eval(z1[1]) - eval(z1[0])
            f.write('{0:.5f} {1:.5f} zlo zhi\n'.format(0.0, gapz))
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
                if eval(l1[0]) == tipnum:
                    w1 = 2
                    j = 1
        elif w1 == 2:
            distance = 0
            allnum = 0
            for i3 in range(degree):
                for i4 in range(eval(b1.group(0))):
                    if i3 != 0 and i3 != degree-1:
                        if i4 + 1 in [xminnum, xmaxnum]:
                            gap = gap + 1
                            j = j + 1
                            continue
                    elif i3 == 0:
                        if i4 + 1 == xmaxnum:
                            gap = gap + 1
                            j = j + 1
                            all_x.append(eval(l[i4][-3]))
                            continue
                    elif i3 == degree - 1:
                        if i4 + 1 == xminnum:
                            gap = gap + 1
                            j = j + 1
                            continue
                    #print(i3, j-gap, eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1]))
                    #if realxmin< = realxmax:
                    f.write('{0} {1} '.format(j-gap, eval(l[i4][1])))
                    f.write('{0} {1:.4f} '.format(eval(l[i4][2]), eval(l[i4][3])))
                    f.write('{0:.5f} '.format(eval(l[i4][-3]) + distance - xlo))
                    f.write('{0:.5f} '.format(eval(l[i4][-2]) - ylo))
                    f.write('{0:.5f}\n'.format(eval(l[i4][-1]) - zlo))
#                    else:
#                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j-gap,\
#eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),
#eval(l[i4][-3]) + distance-xhi, 
#eval(l[i4][-2])-ylo, eval(l[i4][-1])-zlo))
                    j = j + 1
                    if i3 == 0:
                        all_x.append(eval(l[i4][-3]))
                if i3 != 0 and i3 != degree - 1:
                    distance = distance + tipgapx
                    allnum = allnum + num - 2
                else:
                    distance = distance + gapx
                    allnum = allnum + num - 1
            l = []
            j = 1
            w1 = 0
            f.write('\n')
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if eval(l1[0]) == tipnum2:
                    lx = l
                    w2 = 2
                    j = 1
        elif w2 == 2:
            every = 0
            for i3 in range(degree):
                repeated = []
                #print(i3)
                for i4 in range(eval(b2.group(0))):
                    if i3 != 0 and i3 != degree - 1:
                        if eval(l[i4][2]) in [xminnum, xmaxnum]:
                            #print(l[i4])
                            continue
                        elif eval(l[i4][2])>min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        elif eval(l[i4][2])>max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) - 2)
                        if eval(l[i4][3]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            continue
                        elif eval(l[i4][3])>min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                        elif eval(l[i4][3])>max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) - 2)
                    elif i3 == 0:
                        if eval(l[i4][2]) == xmaxnum:
                            continue
                        elif eval(l[i4][2]) > xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                    elif i3 == degree-1:
                        if eval(l[i4][2]) == xminnum:
                            #print('{0} {1}\n'.format(eval(l[i4][2]), eval(l[i4][3])))
                            continue
                        elif eval(l[i4][2]) > xminnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xminnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                    f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                       eval(l[i4][2]) + every, eval(l[i4][3]) + every))
                    #print('{0} {1} {2} {3}'.format(j, eval(l[i4][1]), eval(l[i4][2]) + every, eval(l[i4][3]) + every))
                    j = j + 1
                if i3 != 0:
                    if realxminnum <= realxmaxnum:
                        f.write('{0} {1} {2} {3}\n'.format(j, bondtype, realxminnum + every,\
                                                           realxmaxnum + every - num - 1))
                    else:
                        if i3 == degree - 1:
                            f.write('{0} {1} {2} {3}\n'.format(j, bondtype, realxminnum + every,\
                                                               realxmaxnum + every - num - 1))
                        else:
                            f.write('{0} {1} {2} {3}\n'.format(j, bondtype, realxminnum + every - 1,\
                                                               realxmaxnum + every - num - 1))
                        #if all_x[eval(l[i4][2])-1]< = all_x[eval(l[i4][3])-1]:
                        #    if eval(l[i4][2])>eval(l[i4][3]) + tipnum*i3:
                        #        f.write('{0} {1} {2} {3}\n'.format(j, bondtype,\
                        #                                           realxminnum + every, realxmaxnum + every - num))
                        #    else:
                        #        f.write('{0} {1} {2} {3}\n'.format(j, bondtype,\
                        #                                           realxminnum + every-num, realxmaxnum + every))
                        #elif all_x[eval(l[i4][2])-1]>all_x[eval(l[i4][3])-1]:
                        #    if eval(l[i4][2])>eval(l[i4][3]) + tipnum*i3:
                        #        f.write('{0} {1} {2} {3}\n'.format(j, bondtype,\
                        #                                           realxminnum + every-num, realxmaxnum + every))
                        #    else:
                        #        f.write('{0} {1} {2} {3}\n'.format(j, bondtype,\
                        #                                           realxminnum + every, realxmaxnum + every - num))
                    #print('{0} {1} {2} {3}'.format(j, bondtype,\
                    #                               realxminnum + every-1, realxmaxnum + every - num - 1))
                    j = j + 1
                #l = lx
                if i3 == 0 or i3 == degree - 1:
                    every = every + tipnum - 1
                else:
                    every = every + tipnum - 2
                #print(l[-1], xmaxnum)
                for i4 in range(eval(b2.group(0))):
                    if i3 != 0 and i3 != degree - 1:
                        if (([eval(l[i4][2]), eval(l[i4][3])] == lost2_u
                                 or [eval(l[i4][2]), eval(l[i4][3])] == lost2_d) and
                            [eval(l[i4][2]), eval(l[i4][3])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3])])
                            continue
                        elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 2)
                        if eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        elif eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 2)
                    elif i3 == 0:
                        if (([eval(l[i4][2]), eval(l[i4][3])] == lost2_u
                                 or [eval(l[i4][2]), eval(l[i4][3])] == lost2_d) and
                            [eval(l[i4][2]), eval(l[i4][3])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3])])
                            continue
                        if eval(l[i4][2]) >= xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)  
                        if eval(l[i4][3]) >= xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                    elif i3 == degree-1:
                        if (([eval(l[i4][2]), eval(l[i4][3])] == lost2_u
                                 or [eval(l[i4][2]), eval(l[i4][3])] == lost2_d) and
                            [eval(l[i4][2]), eval(l[i4][3])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3])])
                            continue
                        if eval(l[i4][2]) >= xminnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xminnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                #if i3 == 1:
                    #print(l)
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
                if eval(l1[0]) == tipnum3:
                    w3 = 2
                    j = 1
            #lx = l
        elif w3 == 2:
            every = 0
            for i3 in range(degree):
                repeated = []
                for i4 in range(len(l)):
                    if i3 == 0:
                        if eval(l[i4][2]) == xmaxnum:
                            continue
                        elif eval(l[i4][2]) > xmaxnum:
                            l[i4][2] = str(eval(l[i4][2])-1) 
                        if eval(l[i4][3]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xmaxnum:
                            l[i4][3] = str(eval(l[i4][3])-1)
                        if eval(l[i4][4]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xmaxnum:
                            l[i4][4] = str(eval(l[i4][4])-1)  
                    elif i3 == degree-1:
                        if eval(l[i4][2]) == xminnum:
                            continue
                        elif eval(l[i4][2]) > xminnum:
                            l[i4][2] = str(eval(l[i4][2])-1) 
                        if eval(l[i4][3]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xmaxnum:
                            l[i4][3] = str(eval(l[i4][3])-1)
                        if eval(l[i4][4]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xminnum:
                            l[i4][4] = str(eval(l[i4][4])-1)           
                    else:
                        if eval(l[i4][2]) in [xminnum, xmaxnum]:
                            continue
                        elif eval(l[i4][2])>min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-1)
                        elif eval(l[i4][2])>max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-2)  
                        if eval(l[i4][3]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            continue
                        elif eval(l[i4][3])>min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-1)
                        elif eval(l[i4][3])>max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-2)  
                        if eval(l[i4][4]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            if eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            elif eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            continue
                        elif eval(l[i4][4])>min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4])-1)
                        elif eval(l[i4][4])>max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4])-2)  
                    f.write('{0} {1} {2} {3} {4}\n'.format(j, eval(l[i4][1]),\
                                                           eval(l[i4][2]) + every, eval(l[i4][3]) + every, eval(l[i4][4]) + every))
                    j = j + 1
                if i3 == 0 or i3 == degree-1:
                    every = every + tipnum-1
                else:
                    every = every + tipnum-2
                for i4 in range(len(l)):
                    if i3 == 0:
                        if (([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_u
                                 or [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_d) and
                            [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])])
                            continue
                        elif eval(l[i4][2]) >= xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1) 
                        if eval(l[i4][3]) >= xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)  
                    elif i3 == degree-1:
                        if (([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_u
                                 or [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_d) and
                            [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])])
                            continue
                        elif eval(l[i4][2]) >= xminnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1) 
                        if eval(l[i4][3]) >= xminnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xminnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)           
                    else:
                        if (([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_u
                                 or [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_d) and
                            [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])])
                            continue
                        if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 2)  
                        if eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        elif eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 2)  
                        if eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        elif eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 2)
                #l = lx
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
                if eval(l1[0]) == tipnum4:
                    w4 = 2
                    j = 1
            #lx = l
        elif w4 == 2:
            every = 0
            for i3 in range(degree):
                repeated = []
                for i4 in range(len(l)):
                    if i3 == 0:
                        if eval(l[i4][2]) == xmaxnum:
                            continue
                        elif eval(l[i4][2]) > xmaxnum:
                            l[i4][2] = str(eval(l[i4][2])-1)
                        if eval(l[i4][3]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xmaxnum:
                            l[i4][3] = str(eval(l[i4][3])-1)  
                        if eval(l[i4][4]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xmaxnum:
                            l[i4][4] = str(eval(l[i4][4])-1)
                        if eval(l[i4][5]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) >= xmaxnum:
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            continue
                        elif eval(l[i4][5]) > xmaxnum:
                            l[i4][5] = str(eval(l[i4][5])-1)
                    elif i3 == degree-1:
                        if eval(l[i4][2]) == xminnum:
                            continue
                        elif eval(l[i4][2]) > xminnum:
                            l[i4][2] = str(eval(l[i4][2])-1)  
                        if eval(l[i4][3]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xminnum:
                            l[i4][3] = str(eval(l[i4][3])-1)  
                        if eval(l[i4][4]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xminnum:
                            l[i4][4] = str(eval(l[i4][4])-1)
                        if eval(l[i4][5]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) >= xminnum:
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            continue
                        elif eval(l[i4][5]) > xminnum:
                            l[i4][5] = str(eval(l[i4][5])-1)
                    else:
                        if eval(l[i4][2]) in [xminnum, xmaxnum]:
                            continue
                        elif eval(l[i4][2])>min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-1)
                        elif eval(l[i4][2])>max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-2)  
                        if eval(l[i4][3]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            continue
                        elif eval(l[i4][3])>min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-1)
                        elif eval(l[i4][3])>max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-2)  
                        if eval(l[i4][4]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            if eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            elif eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            continue
                        elif eval(l[i4][4])>min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4])-1)
                        elif eval(l[i4][4])>max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4])-2)  
                        if eval(l[i4][5]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            if eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            elif eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            if eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            elif eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][4] = str(eval(l[i4][2]) + 2)
                            continue
                        elif eval(l[i4][5])>min([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5])-1)
                        elif eval(l[i4][5])>max([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5])-2)
                    f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                               eval(l[i4][2]) + every, eval(l[i4][3]) + every,\
                                                               eval(l[i4][4]) + every, eval(l[i4][5]) + every))
                    j = j + 1
                for i4 in range(len(l)):
                    if i3 == 0:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost4_u
                                 or [eval(l[i4][2]), eval(l[i4][3]),\
                                    eval(l[i4][4]), eval(l[i4][5])] in lost4_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) >= xmaxnum:
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                    elif i3 == degree - 1:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost4_u
                                 or [eval(l[i4][2]), eval(l[i4][3]),\
                                     eval(l[i4][4]), eval(l[i4][5])] in lost4_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= xminnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xminnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xminnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) >= xminnum:
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                    else:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] in lost4_u
                                or [eval(l[i4][2]), eval(l[i4][3]),\
                                    eval(l[i4][4]), eval(l[i4][5])] in lost4_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 2)  
                        if eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        elif eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 2)  
                        if eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        elif eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 2)  
                        if eval(l[i4][5]) >= min([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                        elif eval(l[i4][5]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5]) + 2)

                if i3 == 0 or i3 == degree - 1:
                    every = every + tipnum - 1
                else:
                    every = every + tipnum - 2
                #l = lx
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
                if eval(l1[0]) == tipnum5:
                    w5 = 2
                    j = 1
            #lx = l
        elif w5 == 2:
            every = 0
            for i3 in range(degree):
                repeated = []
                for i4 in range(len(l)):
                    if i3 == 0:
                        if eval(l[i4][2]) == xmaxnum:
                            continue
                        elif eval(l[i4][2]) > xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                        if eval(l[i4][4]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) - 1)
                        if eval(l[i4][5]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) >= xmaxnum:
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            continue
                        elif eval(l[i4][5]) > xmaxnum:
                            l[i4][5] = str(eval(l[i4][5]) - 1)
                    elif i3 == degree-1:
                        if eval(l[i4][2]) == xminnum:
                            continue
                        elif eval(l[i4][2]) > xminnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xminnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                        if eval(l[i4][4]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xminnum:
                            l[i4][4] = str(eval(l[i4][4])-1)
                        if eval(l[i4][5]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) >= xminnum:
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            continue
                        elif eval(l[i4][5]) > xminnum:
                            l[i4][5] = str(eval(l[i4][5])-1)
                    else:
                        if eval(l[i4][2]) in [xminnum, xmaxnum]:
                            continue
                        elif eval(l[i4][2])>min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-1)
                        elif eval(l[i4][2])>max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-2)  
                        if eval(l[i4][3]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            continue
                        elif eval(l[i4][3])>min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-1)
                        elif eval(l[i4][3])>max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-2)  
                        if eval(l[i4][4]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            if eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            elif eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            continue
                        elif eval(l[i4][4])>min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) - 1)
                        elif eval(l[i4][4])>max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) - 2)  
                        if eval(l[i4][5]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            if eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            elif eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            if eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            elif eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][4] = str(eval(l[i4][4]) + 2)
                            continue
                        elif eval(l[i4][5])>min([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5])-1)
                        elif eval(l[i4][5])>max([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5])-2)
                    f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                               eval(l[i4][2]) + every, eval(l[i4][3]) + every,\
                                                               eval(l[i4][4]) + every, eval(l[i4][5]) + every))
                    j = j + 1
                for i4 in range(len(l)):
                    if i3 == 0:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost5_u
                                  or [eval(l[i4][2]), eval(l[i4][3]),\
                                      eval(l[i4][4]), eval(l[i4][5])] in lost5_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) >= xmaxnum:
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                    elif i3 == degree - 1:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost5_u
                                 or [eval(l[i4][2]), eval(l[i4][3]),\
                                     eval(l[i4][4]), eval(l[i4][5])] in lost5_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= xminnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xminnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xminnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) >= xminnum:
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                    else:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost5_u
                                 or [eval(l[i4][2]), eval(l[i4][3]),\
                                     eval(l[i4][4]), eval(l[i4][5])] in lost5_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 2)  
                        if eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        elif eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 2)  
                        if eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        elif eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 2)
                        if eval(l[i4][5]) >= min([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                        elif eval(l[i4][5]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5]) + 2)
                if i3 == 0 or i3 == degree - 1:
                    every = every + tipnum - 1
                else:
                    every = every + tipnum - 2
                #lx = l
            k = k + 1
            l = []
            j = 1
            k = 0
            w4 = 0
            l = lx
            f.write('\n')
        else:
            f.write(line)
    f.close()
    inp.close()
    med.close()
    os.chdir(opath)