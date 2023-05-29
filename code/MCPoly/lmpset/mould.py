import re
import os
import warnings
import math as m
from ase.io import read, write

def replication(filename,x,xgap,y,ygap,z,zgap,absolute=False,startloc=[0,0,0],loc='./'):
    if absolute==True:
        warnings.warn('Make sure you avoid all the overlaps of molecules!')
    opath=os.getcwd()
    os.chdir(loc)
    xloc=startloc[0]
    yloc=startloc[1]
    zloc=startloc[2]
    k=0
    w1=0
    w2=0
    w3=0
    w4=0
    w5=0
    l=[]
    real=[0,0,0]
    try:
        f = open('{0}_{1}{2}{3}.data'.format(filename,x,y,z),'x')
    except:
        f = open('{0}_{1}{2}{3}.data'.format(filename,x,y,z),'w')
    inp = open('{0}.data'.format(filename),'r')
    for line in inp:
        a1=re.search(r' atoms',line)
        a2=re.search(r' bonds',line)
        a3=re.search(r' angles',line)
        a4=re.search(r' dihedrals',line)
        a5=re.search(r' impropers',line)
        a6=re.search(r' xlo xhi',line)
        a7=re.search(r' ylo yhi',line)
        a8=re.search(r' zlo zhi',line)
        c1=re.search('Atoms',line)
        c2=re.search('Bonds',line)
        c3=re.search('Angles',line)
        c4=re.search('Dihedrals',line)
        c5=re.search('Impropers',line)
        if a1:
            b1=re.search(r'[0-9]+',line)
            f.write('{0} atoms\n'.format(eval(b1.group(0))*x*y*z))
            num=b1.group(0)
        elif a2:
            b2=re.search(r'[0-9]+',line)
            f.write('{0} bonds\n'.format(eval(b2.group(0))*x*y*z))
            num2=b2.group(0)
        elif a3:
            b3=re.search(r'[0-9]+',line)
            f.write('{0} angles\n'.format(eval(b3.group(0))*x*y*z))
            num3=b3.group(0)
        elif a4:
            b4=re.search(r'[0-9]+',line)
            f.write('{0} dihedrals\n'.format(eval(b4.group(0))*x*y*z))
            num4=b4.group(0)
        elif a5:
            b5=re.search(r'[0-9]+',line)
            f.write('{0} impropers\n'.format(eval(b5.group(0))*x*y*z))
            num5=b5.group(0)
        elif a6:
            x1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            xd=eval(x1[1])-eval(x1[0])
            if absolute==True:
                real[0]=xd-1/x*xgap
            f.write('{0:.5f} {1:.5f} xlo xhi\n'.format(xloc,xloc+(xd-real[0])*x+xgap*(x-1)))
            if absolute==True:
                real[0]=xd
        elif a7:
            y1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            yd=eval(y1[1])-eval(y1[0])
            if absolute==True:
                real[1]=yd-1/y*ygap
            f.write('{0:.5f} {1:.5f} ylo yhi\n'.format(yloc,yloc+(yd-real[1])*y+ygap*(y-1)))
            if absolute==True:
                real[1]=yd
        elif a8:
            z1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            zd=eval(z1[1])-eval(z1[0])
            if absolute==True:
                real[2]=zd-1/z*zgap
            f.write('{0:.5f} {1:.5f} zlo zhi\n'.format(zloc,zloc+(zd-real[2])*z+zgap*(z-1)))
            if absolute==True:
                real[2]=zd
        elif c1:
            w1=1
            f.write(line)
            f.write('\n')
        elif c2:
            w2=1
            f.write(line)
            f.write('\n')
        elif c3:
            w3=1
            f.write(line)
            f.write('\n')
        elif c4:
            w4=1
            f.write(line)
            f.write('\n')
        elif c5:
            w5=1
            f.write(line)
            f.write('\n')
        elif w1==1:
            l1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            if len(l1)>=6:
                l.append(l1)
                if l1[0]==num:
                    w1=2
                    j=1
        elif w1==2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,eval(l[i4][1]),eval(l[i4][2]),eval(l[i4][3]),eval(l[i4][-3])+(xd-real[0]+xgap)*i1-eval(x1[0]),eval(l[i4][-2])+(yd-real[1]+ygap)*i2-eval(y1[0]),eval(l[i4][-1])+(zd-real[2]+zgap)*i3-eval(z1[0])))
                            j=j+1
            l=[]
            j=1
            w1=0
            f.write('\n')
        elif w2==1:
            l1=re.findall(r'[0-9]+',line)
            if len(l1)==4:
                l.append(l1)
                if l1[0]==num2:
                    w2=2
                    j=1
        elif w2==2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            f.write('{0} {1} {2} {3}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k))
                            j=j+1
                        k=k+1
            l=[]
            j=1
            w2=0
            k=0
            f.write('\n')
        elif w3==1:
            l1=re.findall(r'[0-9]+',line)
            if len(l1)==5:
                l.append(l1)
                if l1[0]==num3:
                    w3=2
                    j=1
        elif w3==2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            f.write('{0} {1} {2} {3} {4}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k))
                            j=j+1
                        k=k+1
            l=[]
            j=1
            k=0
            w3=0
            f.write('\n')
        elif w4==1:
            l1=re.findall(r'[0-9]+',line)
            if len(l1)==6:
                l.append(l1)
                if l1[0]==num4:
                    w4=2
                    j=1
        elif w4==2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k,eval(l[i4][5])+eval(num)*k))
                            j=j+1
                        k=k+1
            l=[]
            j=1
            k=0
            w4=0
            f.write('\n')
        elif w5==1:
            l1=re.findall(r'[0-9]+',line)
            if len(l1)==6:
                l.append(l1)
                if l1[0]==num5:
                    w5=2
                    j=1
        elif w5==2:
            for i1 in range(x):
                for i2 in range(y):
                    for i3 in range(z):
                        for i4 in range(len(l)):
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k,eval(l[i4][5])+eval(num)*k))
                            j=j+1
                        k=k+1
            l=[]
            j=1
            k=0
            w4=0
            f.write('\n')
        else:
            f.write(line)
    f.close()
    os.chdir(opath)

def replication_brick(filename,x,xgap,y,ygap,z,zgap,xpattern='0',ypattern='0',zpattern='0',shuffle=0,absolute=False,startloc=[0,0,0],loc='./'):
    if xpattern=='0' and ypattern=='0' and zpattern=='0':
        warnings.warn("You didn't set any patterns!")
    xloc=startloc[0]
    yloc=startloc[1]
    zloc=startloc[2]
    opath=os.getcwd()
    os.chdir(loc)
    k=0
    w1=0
    w2=0
    w3=0
    w4=0
    w5=0
    half=0.0
    l=[]
    real=[0,0,0]
    pattern=''
    if xpattern=='y':
        pattern='xy'
    if xpattern=='z':
        pattern='xz'
    if ypattern=='x':
        pattern='yx'
    if ypattern=='z':
        pattern='yz'
    if zpattern=='x':
        pattern='zx'
    if zpattern=='y':
        pattern='zy'
    try:
        f = open('{0}_{1}{2}{3}_brick{4}.data'.format(filename,x,y,z,pattern),'x')
    except:
        f = open('{0}_{1}{2}{3}_brick{4}.data'.format(filename,x,y,z,pattern),'w')
    inp = open('{0}.data'.format(filename),'r')
    for line in inp:
        i1,i2,i3=None,None,None
        a1=re.search(r' atoms',line)
        a2=re.search(r' bonds',line)
        a3=re.search(r' angles',line)
        a4=re.search(r' dihedrals',line)
        a5=re.search(r' impropers',line)
        a6=re.search(r' xlo xhi',line)
        a7=re.search(r' ylo yhi',line)
        a8=re.search(r' zlo zhi',line)
        c1=re.search('Atoms',line)
        c2=re.search('Bonds',line)
        c3=re.search('Angles',line)
        c4=re.search('Dihedrals',line)
        c5=re.search('Impropers',line)
        if a1:
            b1=re.search(r'[0-9]+',line)
            if xpattern!='0':
                if xpattern=='y':
                    axpattern='z'
                else:
                    axpattern='y'
                f.write('{0:.0f} atoms\n'.format(eval(b1.group(0))*(x*y*z-int(x/2)*eval(axpattern)-int(eval(axpattern)/2)*(eval(axpattern)%2))))
            if ypattern!='0':
                if ypattern=='z':
                    aypattern='x'
                else:
                    aypattern='z'
                f.write('{0:.0f} atoms\n'.format(eval(b1.group(0))*(x*y*z-int(y/2)*eval(aypattern)-int(eval(aypattern)/2)*(eval(aypattern)%2))))
            if zpattern!='0':
                if zpattern=='y':
                    azpattern='x'
                else:
                    azpattern='y'
                f.write('{0:.0f} atoms\n'.format(eval(b1.group(0))*(x*y*z-int(z/2)*eval(azpattern)-int(eval(azpattern)/2)*(eval(azpattern)%2))))
            num=b1.group(0)
        elif a2:
            b2=re.search(r'[0-9]+',line)
            if xpattern!='0':
                f.write('{0:.0f} bonds\n'.format(eval(b2.group(0))*(x*y*z-int(x/2)*eval(axpattern)-int(eval(axpattern)/2)*(eval(axpattern)%2))))
            if ypattern!='0':
                f.write('{0:.0f} bonds\n'.format(eval(b2.group(0))*(x*y*z-int(y/2)*eval(aypattern)-int(eval(aypattern)/2)*(eval(aypattern)%2))))
            if zpattern!='0':
                f.write('{0:.0f} bonds\n'.format(eval(b2.group(0))*(x*y*z-int(z/2)*eval(azpattern)-int(eval(azpattern)/2)*(eval(azpattern)%2))))
            num2=b2.group(0)
        elif a3:
            b3=re.search(r'[0-9]+',line)
            if xpattern!='0':
                f.write('{0:.0f} angles\n'.format(eval(b3.group(0))*(x*y*z-int(x/2)*eval(axpattern)-int(eval(axpattern)/2)*(eval(axpattern)%2))))
            if ypattern!='0':
                f.write('{0:.0f} angles\n'.format(eval(b3.group(0))*(x*y*z-int(y/2)*eval(aypattern)-int(eval(aypattern)/2)*(eval(aypattern)%2))))
            if zpattern!='0':
                f.write('{0:.0f} angles\n'.format(eval(b3.group(0))*(x*y*z-int(z/2)*eval(azpattern)-int(eval(azpattern)/2)*(eval(azpattern)%2))))
            num3=b3.group(0)
        elif a4:
            b4=re.search(r'[0-9]+',line)
            if xpattern!='0':
                f.write('{0:.0f} dihedrals\n'.format(eval(b4.group(0))*(x*y*z-int(x/2)*eval(axpattern)-int(eval(axpattern)/2)*(eval(axpattern)%2))))
            if ypattern!='0':
                f.write('{0:.0f} dihedrals\n'.format(eval(b4.group(0))*(x*y*z-int(y/2)*eval(aypattern)-int(eval(aypattern)/2)*(eval(aypattern)%2))))
            if zpattern!='0':
                f.write('{0:.0f} dihedrals\n'.format(eval(b4.group(0))*(x*y*z-int(z/2)*eval(azpattern)-int(eval(azpattern)/2)*(eval(azpattern)%2))))
            num4=b4.group(0)
        elif a5:
            b5=re.search(r'[0-9]+',line)
            if xpattern!='0':
                f.write('{0:.0f} impropers\n'.format(eval(b5.group(0))*(x*y*z-int(x/2)*eval(axpattern)-int(eval(axpattern)/2)*(eval(axpattern)%2))))
            if ypattern!='0':
                f.write('{0:.0f} impropers\n'.format(eval(b5.group(0))*(x*y*z-int(y/2)*eval(aypattern)-int(eval(aypattern)/2)*(eval(aypattern)%2))))
            if zpattern!='0':
                f.write('{0:.0f} impropers\n'.format(eval(b5.group(0))*(x*y*z-int(z/2)*eval(azpattern)-int(eval(azpattern)/2)*(eval(azpattern)%2))))
            num5=b5.group(0)
        elif a6:
            x1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            xd=eval(x1[1])-eval(x1[0])
            if absolute==True:
                real[0]=xd-1/x*xgap
            f.write('{0:.5f} {1:.5f} xlo xhi\n'.format(xloc,xloc+(xd-real[0])*x+xgap*(x-1)))
            if absolute==True:
                real[0]=xd
        elif a7:
            y1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            yd=eval(y1[1])-eval(y1[0])
            if absolute==True:
                real[1]=yd-1/y*ygap
            f.write('{0:.5f} {1:.5f} ylo yhi\n'.format(yloc,yloc+(yd-real[1])*y+ygap*(y-1)))
            if absolute==True:
                real[1]=yd
        elif a8:
            z1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            zd=eval(z1[1])-eval(z1[0])
            if absolute==True:
                real[2]=zd-1/z*zgap
            f.write('{0:.5f} {1:.5f} zlo zhi\n'.format(zloc,zloc+(zd-real[2])*z+zgap*(z-1)))
            if absolute==True:
                real[2]=zd
        elif c1:
            w1=1
            f.write(line)
            f.write('\n')
        elif c2:
            w2=1
            f.write(line)
            f.write('\n')
        elif c3:
            w3=1
            f.write(line)
            f.write('\n')
        elif c4:
            w4=1
            f.write(line)
            f.write('\n')
        elif c5:
            w5=1
            f.write(line)
            f.write('\n')
        elif w1==1:
            l1=re.findall(r'\-?[0-9]+\.?[0-9]*',line)
            if len(l1)>=6:
                l.append(l1)
                if l1[0]==num:
                    w1=2
                    j=1
        elif w1==2:
            if xpattern!='0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern=='y':
                            if i3==None:
                                i3=1
                            half=i1%2/2
                        for i3 in range(z):
                            if xpattern=='y':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if xpattern=='z':
                                half=i1%2/2
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            for i4 in range(len(l)):
                                if xpattern=='y':
                                    if half==0.5 and i2==y-1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,eval(l[i4][1]),eval(l[i4][2]),eval(l[i4][3]),eval(l[i4][-3])+(xd+xgap-real[0])*i1-eval(x1[0]),eval(l[i4][-2])+(yd+ygap-real[1])*(i2+half)-eval(y1[0]),eval(l[i4][-1])+(zd+zgap-real[2])*i3-eval(z1[0])))
                                        j=j+1
                                elif xpattern=='z':
                                    if half==0.5 and i3==z-1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,eval(l[i4][1]),eval(l[i4][2]),eval(l[i4][3]),eval(l[i4][-3])+(xd+xgap-real[0])*i1-eval(x1[0]),eval(l[i4][-2])+(yd+ygap-real[1])*i2-eval(y1[0]),eval(l[i4][-1])+(zd+zgap-real[2])*(i3+half)-eval(z1[0])))
                                        j=j+1
            elif ypattern!='0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern=='x':
                            if i3==None:
                                i3=1
                            half=i2%2/2
                        for i3 in range(z):
                            if ypattern=='x':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if ypattern=='z':
                                half=i2%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            for i4 in range(len(l)):
                                if ypattern=='x':
                                    if half==0.5 and i1==x-1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,eval(l[i4][1]),eval(l[i4][2]),eval(l[i4][3]),eval(l[i4][-3])+(xd+xgap-real[0])*(i1+half)-eval(x1[0]),eval(l[i4][-2])+(yd+ygap-real[1])*i2-eval(y1[0]),eval(l[i4][-1])+(zd+zgap-real[2])*i3-eval(z1[0])))
                                        j=j+1
                                elif ypattern=='z':
                                    if half==0.5 and i3==z-1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,eval(l[i4][1]),eval(l[i4][2]),eval(l[i4][3]),eval(l[i4][-3])+(xd+xgap-real[0])*i1-eval(x1[0]),eval(l[i4][-2])+(yd+ygap-real[1])*i2-eval(y1[0]),eval(l[i4][-1])+(zd+zgap-real[2])*(i3+half)-eval(z1[0])))
                                        j=j+1
            elif zpattern!='0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern=='x':
                            if i2==None:
                                i2=1
                            half=i3%2/2
                        for i2 in range(y):
                            if zpattern=='x':
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            if zpattern=='y':
                                half=i3%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            for i4 in range(len(l)):
                                if zpattern=='x':
                                    if half==0.5 and i1==x-1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,eval(l[i4][1]),eval(l[i4][2]),eval(l[i4][3]),eval(l[i4][-3])+(xd-real[0]+xgap)*(i1+half)-eval(x1[0]),eval(l[i4][-2])+(yd-real[1]+ygap)*i2-eval(y1[0]),eval(l[i4][-1])+(zd-real[2]+zgap)*i3-eval(z1[0])))
                                        j=j+1
                                elif zpattern=='y':
                                    if half==0.5 and i2==y-1:
                                        continue
                                    else:
                                        f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j,eval(l[i4][1]),eval(l[i4][2]),eval(l[i4][3]),eval(l[i4][-3])+(xd-real[0]+xgap)*i1-eval(x1[0]),eval(l[i4][-2])+(yd-real[1]+ygap)*(i2+half)-eval(y1[0]),eval(l[i4][-1])+(zd-real[2]+zgap)*i3-eval(z1[0])))
                                        j=j+1
            l=[]
            j=1
            w1=0
            f.write('\n')
        elif w2==1:
            l1=re.findall(r'[0-9]+',line)
            if len(l1)==4:
                l.append(l1)
                if l1[0]==num2:
                    w2=2
                    j=1
        elif w2==2:
            if xpattern!='0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern=='y':
                            if i3==None:
                                i3=1
                            half=i1%2/2
                        for i3 in range(z):
                            if xpattern=='y':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if xpattern=='z':
                                half=i1%2/2
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            if half==0.5 and i2==y-1 and xpattern=='y':
                                continue
                            if half==0.5 and i3==z-1 and xpattern=='z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k))
                                j=j+1
                            k=k+1
            if ypattern!='0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern=='x':
                            if i3==None:
                                i3=1
                            half=i2%2/2
                        for i3 in range(z):
                            if ypattern=='x':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if ypattern=='z':
                                half=i2%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            if half==0.5 and i1==x-1 and ypattern=='x':
                                continue
                            if half==0.5 and i3==z-1 and ypattern=='z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k))
                                j=j+1
                            k=k+1
            if zpattern!='0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern=='x':
                            if i2==None:
                                i2=1
                            half=i3%2/2
                        for i2 in range(y):
                            if zpattern=='x':
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            if zpattern=='y':
                                half=i3%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            if half==0.5 and i1==x-1 and zpattern=='x':
                                continue
                            if half==0.5 and i2==y-1 and zpattern=='y':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k))
                                j=j+1
                            k=k+1
            l=[]
            j=1
            w2=0
            k=0
            f.write('\n')
        elif w3==1:
            l1=re.findall(r'[0-9]+',line)
            if len(l1)==5:
                l.append(l1)
                if l1[0]==num3:
                    w3=2
                    j=1
        elif w3==2:
            if xpattern!='0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern=='y':
                            if i3==None:
                                i3=1
                            half=i1%2/2
                        for i3 in range(z):
                            if xpattern=='y':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if xpattern=='z':
                                half=i1%2/2
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            if half==0.5 and i2==y-1 and xpattern=='y':
                                continue
                            if half==0.5 and i3==z-1 and xpattern=='z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k))
                                j=j+1
                            k=k+1
            if ypattern!='0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern=='x':
                            if i3==None:
                                i3=1
                            half=i2%2/2
                        for i3 in range(z):
                            if ypattern=='x':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if ypattern=='z':
                                half=i2%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            if half==0.5 and i1==x-1 and ypattern=='x':
                                continue
                            if half==0.5 and i3==z-1 and ypattern=='z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k))
                                j=j+1
                            k=k+1
            if zpattern!='0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern=='x':
                            if i2==None:
                                i2=1
                            half=i3%2/2
                        for i2 in range(y):
                            if zpattern=='x':
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            if zpattern=='y':
                                half=i3%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            if half==0.5 and i1==x-1 and zpattern=='x':
                                continue
                            if half==0.5 and i2==y-1 and zpattern=='y':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k))
                                j=j+1
                            k=k+1
            l=[]
            j=1
            k=0
            w3=0
            f.write('\n')
        elif w4==1:
            l1=re.findall(r'[0-9]+',line)
            if len(l1)==6:
                l.append(l1)
                if l1[0]==num4:
                    w4=2
                    j=1
        elif w4==2:
            if xpattern!='0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern=='y':
                            if i3==None:
                                i3=1
                            half=i1%2/2
                        for i3 in range(z):
                            if xpattern=='y':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if xpattern=='z':
                                half=i1%2/2
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            if half==0.5 and i2==y-1 and xpattern=='y':
                                continue
                            if half==0.5 and i3==z-1 and xpattern=='z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k,eval(l[i4][5])+eval(num)*k))
                                j=j+1
                            k=k+1
            if ypattern!='0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern=='x':
                            if i3==None:
                                i3=1
                            half=i2%2/2
                        for i3 in range(z):
                            if ypattern=='x':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if ypattern=='z':
                                half=i2%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            if half==0.5 and i1==x-1 and ypattern=='x':
                                continue
                            if half==0.5 and i3==z-1 and ypattern=='z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k,eval(l[i4][5])+eval(num)*k))
                                j=j+1
                            k=k+1
            if zpattern!='0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern=='x':
                            if i2==None:
                                i2=1
                            half=i3%2/2
                        for i2 in range(y):
                            if zpattern=='x':
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            if zpattern=='y':
                                half=i3%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            if half==0.5 and i1==x-1 and zpattern=='x':
                                continue
                            if half==0.5 and i2==y-1 and zpattern=='y':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k,eval(l[i4][5])+eval(num)*k))
                                j=j+1
                            k=k+1

            l=[]
            j=1
            k=0
            w4=0
            f.write('\n')
        elif w5==1:
            l1=re.findall(r'[0-9]+',line)
            if len(l1)==6:
                l.append(l1)
                if l1[0]==num5:
                    w5=2
                    j=1
        elif w5==2:
            if xpattern!='0':
                for i1 in range(x):
                    for i2 in range(y):
                        if xpattern=='y':
                            if i3==None:
                                i3=1
                            half=i1%2/2
                        for i3 in range(z):
                            if xpattern=='y':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if xpattern=='z':
                                half=i1%2/2
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            if half==0.5 and i2==y-1 and xpattern=='y':
                                continue
                            if half==0.5 and i3==z-1 and xpattern=='z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k,eval(l[i4][5])+eval(num)*k))
                                j=j+1
                            k=k+1
            if ypattern!='0':
                for i2 in range(y):
                    for i1 in range(x):
                        if ypattern=='x':
                            if i3==None:
                                i3=1
                            half=i2%2/2
                        for i3 in range(z):
                            if ypattern=='x':
                                if shuffle==1 and i3%2==1:
                                    half=0.5-half
                            if ypattern=='z':
                                half=i2%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            if half==0.5 and i1==x-1 and ypattern=='x':
                                continue
                            if half==0.5 and i3==z-1 and ypattern=='z':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k,eval(l[i4][5])+eval(num)*k))
                                j=j+1
                            k=k+1
            if zpattern!='0':
                for i3 in range(z):
                    for i1 in range(x):
                        if zpattern=='x':
                            if i2==None:
                                i2=1
                            half=i3%2/2
                        for i2 in range(y):
                            if zpattern=='x':
                                if shuffle==1 and i2%2==1:
                                    half=0.5-half
                            if zpattern=='y':
                                half=i3%2/2
                                if shuffle==1 and i1%2==1:
                                    half=0.5-half
                            if half==0.5 and i1==x-1 and zpattern=='x':
                                continue
                            if half==0.5 and i2==y-1 and zpattern=='y':
                                continue
                            for i4 in range(len(l)):
                                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j,eval(l[i4][1]),eval(l[i4][2])+eval(num)*k,eval(l[i4][3])+eval(num)*k,eval(l[i4][4])+eval(num)*k,eval(l[i4][5])+eval(num)*k))
                                j=j+1
                            k=k+1
            l=[]
            j=1
            k=0
            w4=0
            f.write('\n')
        else:
            f.write(line)
    f.close()
    os.chdir(opath)

class mould:
    """
    The definition of molecule system.
    mould(polymer, loc='./')
    atoms: Your molecule system name on your .data file.
    loc: File Location. The default is your current location.
    You can get the further information by .cube and .brick.
    """
    
    def __init__(self, atoms, loc='./'):
        self.atoms=atoms
        self.loc=loc
    
    def cube(self,x,xgap,y,ygap,z,zgap,absolute=False,startloc=[0,0,0]):
        """
    The method to create a box with replicated molecular, or create a crystal
    * * * * *
    * * * * *
    * * * * *
    * * * * *     Cube Structure
    
    cube(x,xgap,y,ygap,z,zgap,absolute=False,startloc=[0,0,0])
    x: number axis in x axis
    xgap: distance of nearest atoms between each molecule systems in x axis direction
    y: number axis in y axis
    ygap: distance of nearest atoms between each molecule systems in y axis direction
    z: number axis in z axis
    zgap: distance of nearest atoms between each molecule systems in z axis direction
    absolute: Use the distance between the molecule systems box instead of nearest atoms in xgap, ygap and zgap. The default is True.
    startloc: The location of the beginning the replicated system drawing. Mostly, it is minimum of x, y, z in molecule system.
    Example:
        Input:
            from MCPoly.lmpset import mould
            atoms=mould('Poly1')
            atoms.cube(6,5,6,5,3,5)
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
        return replication(self.atoms,x,xgap,y,ygap,z,zgap,absolute,startloc,self.loc)
    
    def brick(self,x,xgap,y,ygap,z,zgap,xpattern='0',ypattern='0',zpattern='0',shuffle=0,absolute=False,startloc=[0,0,0]):
        """
    The method to create a box with replicated molecular in brick settlement.
    * * * * *
     * * * *
    * * * * *
     * * * *       Brick Structure
    
    brick(x,xgap,y,ygap,z,zgap,xpattern='0',ypattern='0',zpattern='0',shuffle=0,absolute=False,startloc=[0,0,0])
    x: number axis in x axis
    xgap: distance of nearest atoms between each molecule systems in x axis direction
    y: number axis in y axis
    ygap: distance of nearest atoms between each molecule systems in y axis direction
    z: number axis in z axis
    zgap: distance of nearest atoms between each molecule systems in z axis direction
    (x,y,z)pattern: To show the direction of you brick patterns. e.g. xpattern='y' mean brick patterns on xy plane and other horizental planes, with linear aligns in y direction. 
    shuffle: With two reversed brick pattern on alternate planes. The default is false.
    absolute: Use the distance between the molecule systems box instead of nearest atoms in xgap, ygap and zgap. The default is True.
    startloc: The location of the beginning the replicated system drawing. Mostly, it has minimum of x, y, z.
    TIPS: You can only use  one of (x,y,z)pattern.
    Example:
        Input:
            from MCPoly.lmpset import mould
            atoms=mould('Poly1')
            atoms.brick(6,5,6,5,4,5,xpattern='y')
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
        
        if xpattern!='0' and ypattern!='0':
            raise AssertionError('You can only use  one of (x,y,z)pattern.')
        if xpattern!='0' and zpattern!='0':
            raise AssertionError('You can only use  one of (x,y,z)pattern.')
        if ypattern!='0' and zpattern!='0':
            raise AssertionError('You can only use  one of (x,y,z)pattern.')
        return replication_brick(self.atoms,x,xgap,y,ygap,z,zgap,xpattern,ypattern,zpattern,shuffle,absolute,startloc,self.loc)