import glob
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.optimize import curve_fit

def power_law(xi,A,B):
    return A*xi**B

file=[]
nf=0
files = sorted(glob.glob('rcoree*.dat'))
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
prob=np.zeros(nf)
size=np.zeros(nf)
for i in range(nf):
    data = np.loadtxt(fname = file[i])
    xx=[]; yy=[]; ddy=[]
    for j in range(len(data)):
        xx.append(data[j][0])
        yy.append(data[j][1])
        ddy.append(data[j][2])
    x.append(xx)
    y.append(yy)
    dy.append(ddy)
    #size[i]=data[j][4]
    #prob[i]=data[j][4]
#print(alldata[1])
#print(x[1],size[1],prob[1])


#nump=len(pfile)


#A1, B1 = curve_fit(power_law, x1, y1)[0]
#---find same prob---
#plt.figure(figsize=(6,6))
marker = itertools.cycle(('p', '^', '.', 'o', '*'))
color1 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange'))
color2 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange'))
j=0
for j in range(nf):
    for kk in range(100):
        if(x[j][kk] == 16):
            #print(1)
            #print(min(x[j]),max(x[j]))
            AA, BB = curve_fit(power_law, x[j], y[j])[0]
            xf=np.arange((x[j][kk])/1.2,max(x[j])*1.2,2.0)
            AA=x[j][kk]*y[j][kk]
            BB=-1
            yf=AA*xf**BB#*1.1
            break

    #plt.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker), label='SDRG ly=%2i' % size[j],
    #                     linestyle='',color=next(color1),capsize=2)
    plt.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker), label='C(L)',
                         linestyle='',color=next(color1),capsize=2)
    plt.plot(xf,yf,linestyle='-.',color=next(color2),label='y=%5.3fx^%5.3f' % (AA,BB))
plt.yscale('log')
plt.xscale('log')

plt.title('end-to-end spins correlation,SDRG')
plt.xlabel('$L$')
plt.ylabel('$C(L)$')
plt.legend()
plt.savefig('sdrgCor.png')
plt.savefig('sdrgCor.pdf')
plt.show()
            
    
#for i in range(nf):
#    if(
#    plt.figure(figsize=(6,6))
#    plt.axes( xlim=(min(x[i])/1.2,max(x[i])*1.2),ylim=(min(y[i])/2.,max(y[i])*5/4.) )
