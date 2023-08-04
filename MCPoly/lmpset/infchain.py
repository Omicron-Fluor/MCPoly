import re
import os
import warnings
import math as m
from ase.io import read,  write

def infchain(file, bondtype, loc = './'):
    """
    The method to create a single molecule for periodical chain, specialised for polymers.
    infchain(file, bondtype, loc = './')
    file: Your molecule system name on your .data file.
    bondtype: The bond type between the start and the end of the polymer. 
    loc: File Location. The default is your current location.
    Example:
        Input:
            from MCPoly.lmpset import chain
            infchain('Poly1', 11)
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
    l = []
    real = [0, 0, 0]
    inp  =  open('{0}.data'.format(file), 'r')
    for line in inp:
        a1 = re.search(r' atoms', line)
        a2 = re.search(r' bonds', line)
        a3 = re.search(r' angles', line)
        a4 = re.search(r' dihedrals', line)
        a5 = re.search(r' impropers', line)
        m1 = re.search(r'Masses', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            num = b1.group(0)
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            num2 = b2.group(0)
        elif a3:
            b3 = re.search(r'[0-9]+', line)
            num3 = b3.group(0)
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            num4 = b4.group(0)
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            num5 = b5.group(0)
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
                if l1[0] == num:
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
                if eval(l[i4][-3]) < xmin:
                    xmin = eval(l[i4][-3])
                    xminnum = eval(l[i4][0])
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
                if l1[0] == num2:
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
                if l1[0] == num3:
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
                j = j+1
            k = k+1
            l = []
            t3 = j
            j = 1
            k = 0
            w3 = 0
        elif w4 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num4:
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
                j = j+1
            k = k+1
            l = []
            t4 = j
            j = 1
            k = 0
            w4 = 0
        elif w5 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num5:
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
    inp.close()
#    print(t2, t3, t4, t5)
    inp  =  open('{0}.data'.format(file), 'r')
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    gap = 0
    l = []
    try:
        f  =  open('{0}_Chain.data'.format(file), 'x')
    except:
        f  =  open('{0}_Chain.data'.format(file), 'w')
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
            f.write('{0} atoms\n'.format(eval(b1.group(0))-2))
            num = b1.group(0)
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            f.write('{0} bonds\n'.format(eval(b2.group(0))-1))
            num2 = b2.group(0)
        elif a3:
            f.write('{0} angles\n'.format(t3-1))
            num3 = b3.group(0)
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            f.write('{0} dihedrals\n'.format(t4-1))
            num4 = b4.group(0)
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            f.write('{0} impropers\n'.format(t5-1))
            num5 = b5.group(0)
        elif a6:
            x1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            xd = eval(x1[1])-eval(x1[0])
            f.write('{0:.5f} {1:.5f} xlo xhi\n'.format(xmin, xmax))
        elif a7:
            y1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            yd = eval(y1[1])-eval(y1[0])
            f.write('{0:.5f} {1:.5f} ylo yhi\n'.format(eval(y1[0]), eval(y1[1])))
        elif a8:
            z1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            zd = eval(z1[1])-eval(z1[0])
            f.write('{0:.5f} {1:.5f} zlo zhi\n'.format(eval(z1[0]), eval(z1[1])))
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
            for i4 in range(eval(b1.group(0))-2+1):
                if j in [xminnum, xmaxnum]:
                    gap = gap + 1
                    j = j + 1
                    continue
                f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j-gap,\
                                                                               eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),\
                                                                               eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])))
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
            for i4 in range(eval(b2.group(0))-1+1):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-1)
                elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-2)  
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-1)
                elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-2)  
                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                   eval(l[i4][2]), eval(l[i4][3])))
                j = j + 1
            k = k + 1
            l = []
            if realxmaxnum in [xminnum, xmaxnum]:
                continue
            elif realxmaxnum > min([xminnum, xmaxnum]):
                realxmaxnum = realxmaxnum-1
            elif realxmaxnum > max([xminnum, xmaxnum]):
                realxmaxnum = realxmaxnum-2  
            if realxminnum in [xminnum, xmaxnum]:
                continue
            elif realxminnum > min([xminnum, xmaxnum]):
                realxminnum = realxminnum-1
            elif realxminnum > max([xminnum, xmaxnum]):
                realxminnum = realxminnum-2  
            f.write('{0} {1} {2} {3}\n'.format(j, bondtype,\
                                               realxmaxnum, realxminnum))
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
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-1)
                elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-2)  
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-1)
                elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-2)  
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][4]) > min([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-1)
                elif eval(l[i4][4]) > max([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-2)  
                f.write('{0} {1} {2} {3} {4}\n'.format(j, eval(l[i4][1]),\
                                                       eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])))
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
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-1)
                elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-2)  
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-1)
                elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-2)  
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][4]) > min([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-1)
                elif eval(l[i4][4]) > max([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-2)  
                if eval(l[i4][5]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][5]) > min([xminnum, xmaxnum]):
                    l[i4][5] = str(eval(l[i4][5])-1)
                elif eval(l[i4][5]) > max([xminnum, xmaxnum]):
                    l[i4][5] = str(eval(l[i4][5])-2)
                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                           eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4]), eval(l[i4][5])))
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
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-1)
                elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-2)  
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-1)
                elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-2)  
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][4]) > min([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-1)
                elif eval(l[i4][4]) > max([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-2)  
                if eval(l[i4][5]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][5]) > min([xminnum, xmaxnum]):
                    l[i4][5] = str(eval(l[i4][5])-1)
                elif eval(l[i4][5]) > max([xminnum, xmaxnum]):
                    l[i4][5] = str(eval(l[i4][5])-2)
                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                           eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4]), eval(l[i4][5])))
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
    inp.close()
    os.chdir(opath)