import glob
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.optimize import curve_fit

from matplotlib import rcParams
config = {
        "font.family":'serif',
        "font.size": 10,
        "mathtext.fontset" :'stix',
        #"font.serif": ['SimSun'],
        }
rcParams.update(config)
def power_law(xi,A,B):
    return A*xi**B

file=[]
nf=0
files = sorted(glob.glob('rcorp20.dat'))
for ff in files:
    file.append(ff)
    nf=nf+1
alldata=[]
for i in range(nf):
    data = np.loadtxt(fname = file[i])
    alldata.append(data)
print(alldata)

x=[]
y=[]
dy=[]
tcor=[]
prob=np.zeros(nf)
size=np.zeros(nf)
for i in range(nf):
    data = np.loadtxt(fname = file[i])
    xx=[]; yy=[]; ddy=[]; tt=[]
    for j in range(len(data)):
        xx.append(data[j][0])
        yy.append(data[j][1])
        ddy.append(data[j][2])
        tt.append(data[j][4])
    x.append(xx)
    y.append(yy)
    dy.append(ddy)
    tcor.append(tt)
    #size[i]=data[j][4]
    #prob[i]=data[j][4]
#print(alldata[1])
#print(x[1],size[1],prob[1])


#nump=len(pfile)


#A1, B1 = curve_fit(power_law, x1, y1)[0]
#---find same prob---
fig=plt.figure(figsize=(18,6))
marker = itertools.cycle(('p', '^', '.', 'o', '*'))
color1 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange'))
color2 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange'))
j=0
for j in range(nf):
    for kk in range(100):
        if(x[j][kk] == x[j][23]):
            #print(1)
            #print(min(x[j]),max(x[j]))
            AA, BB = curve_fit(power_law, x[j], y[j])[0]
            xf=np.arange((x[j][kk])/1.2,max(x[j])*1.2,2.0)
            AA=x[j][kk]*y[j][kk]
            BB=-1
            yf=AA*xf**BB#*1.1
            break

    #plt.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker), label='$P_{rm}=1-\\frac{20}{L}$',
    #                     linestyle='',color=next(color1),capsize=2)
    left,bottom,width,height=0.1,0.1,0.8,0.8   #表示从整张图的左边10%，下边10%的位置开始，宽度、高度分别是是整个图的80%
    ax1=fig.add_axes([left,bottom,width,height])
    ax1.plot(x[j],y[j],'r',marker='o',label='$P_{rm}=1-\\frac{20}{L}$')   #'r'表示红色的线
    #ax1.plot(xf,yf,linestyle='-.',color=next(color2),label='$L^{-1}$')
    kk=0
    for a, b, c in zip(x[j], y[j], tcor[j]):
        kk+=1
        ax1.text(a, b+0.2*(-1)**kk*b, '%.2f' % (1-c*1./a), ha='center', va='bottom', fontsize=6, color='blue')

    left,bottom,width,height=0.2,0.2,0.2,0.2   #表示从整张图的左边10%，下边10%的位置开始，宽度、高度分别是是整个图的80%
    ax2=fig.add_axes([left,bottom,width,height])
    x1=[]
    y1=[]
    xx1=[]
    yy1=[]
    for kk in range(21):
        x1.append(x[j][kk])
        y1.append(y[j][kk])
    for zz in range(7,21):
        xx1.append(x[j][zz])
        yy1.append(y[j][zz])
    for zz in range(len(xx1)):
        #AA, BB = curve_fit(power_law, xx1, yy1)[0]
        xf=np.arange(xx1[0],96,2.0)
        BB=-0.5
        AA=xx1[9]**-BB*yy1[9]
        yf=AA*xf**BB#*1.1
        break
    ax2.plot(x1,y1,'b',marker='o')
    ax2.plot(xf,yf,linestyle='-.',color=next(color2),label='$L^{%.2f}$' % BB)

    left,bottom,width,height=0.6,0.6,0.2,0.2   #表示从整张图的左边10%，下边10%的位置开始，宽度、高度分别是是整个图的80%
    ax3=fig.add_axes([left,bottom,width,height])
    x2=[]
    y2=[]
    for kk in range(20,28):
        x2.append(x[j][kk])
        y2.append(y[j][kk])
    for kk in range(len(x2)):
            #AA, BB = curve_fit(power_law, x[j], y[j])[0]
        #AA, BB = curve_fit(power_law, x2, y2)[0]
        xf=np.arange(96,1024,2.0)
        BB=-1
        AA=x2[3]**-BB*y2[3]
        yf=AA*xf**BB#*1.1
        break
    ax3.plot(x2,y2,'g',marker='o')
    ax3.plot(xf,yf,linestyle='-.',color='g',label='$L^{%.2f}$' % BB)
    

#plt.yscale('log')
#plt.xscale('log')
#plt.ylim([2*10**-5,2*10**-3])
ax1.set_ylim([2*10**-5,2*10**-3])
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax1.legend(prop={'size':8},loc=3)
ax2.legend(prop={'size':8},loc=3)
ax3.legend(prop={'size':8},loc=3)
ax1.set_xlabel(r'$L$')
ax1.set_ylabel(r'$C(L)$')
ax2.set_xlabel(r'$L$')
ax2.set_ylabel(r'$C(L)$')
ax3.set_xlabel(r'$L$')
ax3.set_ylabel(r'$C(L)$')

#plt.title('end-to-end spins correlation,SDRG')
plt.xlabel('$L$')
plt.ylabel('$C(L)$')
plt.legend()
plt.savefig('sdrgl20Cor.png')
plt.savefig('sdrgl20Cor.pdf',bbox_inches="tight")
plt.show()
            
    
#for i in range(nf):
#    if(
#    plt.figure(figsize=(6,6))
#    plt.axes( xlim=(min(x[i])/1.2,max(x[i])*1.2),ylim=(min(y[i])/2.,max(y[i])*5/4.) )
