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
#files = sorted(glob.glob('p*f*_cor.dat'))
files = sorted(glob.glob('./mp*128locsus.dat'))
for ff in files:
    file.append(ff)
    nf=nf+1
    print(ff)
alldata=[]
for i in range(nf):
    data = np.loadtxt(fname = file[i])
    alldata.append(data)
print(alldata)

print(file)
x=[]
y=[]
dy=[]
prob=np.zeros(nf)
size=np.zeros(nf)
for i in range(nf):
    data = np.loadtxt(fname = file[i])
    xx=[]; yy=[]; ddy=[]
    for j in range(len(data)):
        xx.append(1./data[j][0])
        yy.append(data[j][1])
        ddy.append(data[j][2])
    x.append(xx)  #T=1/beta
    y.append(yy)
    dy.append(ddy)
    size[i]=data[j][4]
    prob[i]=data[j][5]
#print(alldata[1])
#print(x[1],size[1],prob[1])


pfile = np.loadtxt( fname = 'lseed.in' ) #pnum.in write the total number of l seed
print(pfile)
numl=len(pfile)
numl=pfile.size
print('numl',numl)
print(np.shape(pfile))


#A1, B1 = curve_fit(power_law, x1, y1)[0]
#---find same prob---
plt.figure(figsize=(8,4))
marker = itertools.cycle(('p', '^', 'd', 'o', 'D'))
color1 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange','purple','olive','tomato','deepskyblue','orchid'))
color11 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange','purple','olive','tomato','deepskyblue','orchid'))
markerface = itertools.cycle(('black', 'red', 'blue', 'green', 'orange','purple','olive','tomato','deepskyblue','orchid','none','none','none','none','none','none'))
#color2 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange'))
#plt.errorbar(x[0],y[0],yerr=dy[0],marker=next(marker),markerfacecolor=next(markerface), label='16' ,linestyle='-',color=next(color1),capsize=2)
#               #plt.plot(xf,yf,linestyle='-.',color=next(color2))#,label='y=%5.3fx^%5.3f' % (AA,BB))
#plt.yscale('log')
#plt.xscale('log')
print(prob)
for i in range(numl):
    for j in range(nf):
        if(numl == 1):
            if(size[j] == pfile):
                print(size[j])
                for kk in range(100):
                    if(x[j][kk] == 16):
                        #print(1)
                        #print(min(x[j]),max(x[j]))
                        AA, BB = curve_fit(power_law, x[j], y[j])[0]
                        xf=np.arange((x[j][kk])/1.2,max(x[j])*1.2,2.0)
                        AA=x[j][kk]*y[j][kk]
                        BB=-1
                        yf=AA*xf**BB
                        break
        else:
         #   if(size[j] == pfile[i]):
            if(prob[j] == pfile[i]):
                print(prob[j])
            #    for kk in range(100):
            #        if(x[j][kk] == 64):
            #            #print(1)
            #            #print(min(x[j]),max(x[j]))
            #            AA, BB = curve_fit(power_law, x[j], y[j])[0]
            #            xf=np.arange((x[j][kk])/1.2,max(x[j])*1.2,2.0)
            #            AA=x[j][kk]*y[j][kk]
            #            BB=-1
            #            yf=AA*xf**BB
            #            break
            
               # if(prob[j] == 1.0):
               #     for kk in range(100):
               #         if(x[j][kk] == 1./64.):
               #             temp=1./64.
               #             xaa=y[j][kk]
               #             conf=xaa*temp*(np.log(temp))**2
               #             xf=np.arange((x[j][kk])/5.,max(x[j])*1.2,0.0001)
               #             yf=conf/(xf*(np.log(xf)**2))
               #             print(conf)
               #             break
               #     plt.plot(xf,yf,linestyle='-.',color='black',label='$\\frac{1}{Tln^2|T|}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
               #     print(prob[j])

                for kk in range(100):
                    if(x[j][kk] == 1./128.):
                         temp=1./128.
                         xaa=y[j][kk]
                         conf=xaa*temp*(np.log(temp))**2
                         xf=np.arange((x[j][kk])/5.,max(x[j])*1.2,0.0001)
                         yf=conf/(xf*(np.log(xf)**2))
                         print(conf)
                         break
                if(prob[j] == 1.0):
                    plt.plot(xf,yf,linestyle='-.',color=next(color11),label='$\\frac{1}{Tln^2|T|}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
                else:
                    plt.plot(xf,yf,linestyle='-.',color=next(color11))#,label='y=%5.3fx^%5.3f' % (AA,BB))
                #plt.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(markerface), label='$L=$%3.3f' % size[j],linestyle='-',color=next(color1),capsize=2)
                plt.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(markerface), label='$P_{rm}=$%1.3f' % prob[j],linestyle='-',color=next(color1),capsize=2)
#                plt.plot(xf,yf,linestyle='-.',color=next(color2))#,label='y=%5.3fx^%5.3f' % (AA,BB))
                plt.yscale('log')
                plt.xscale('log')
                #plt.xlim(8e-3,0.2)
                #plt.ylim(0.2,5)
    #plt.title('ladder(obc),random J,rm size=%1.2f' % (pfile[i]))
#
#fitx = np.arange(60,130,0.1)
#A=10**(-2)*60.
#plt.plot(fitx,A*fitx**(-1),linestyle='-.',color='black',label='$L^{-1}$')
plt.title('2-leg ladder(obc) $\chi_{loc}$, L=128 ,random coupling,SSE')
#plt.xlabel(chr(946))
plt.xlabel('$T$')
plt.ylabel('$\chi_{loc}$')
plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)
plt.tight_layout() #避免被裁匡
#plt.savefig('eqCor080.png')
plt.savefig('diagram_mlocsus.pdf')
plt.show()
            
    
#for i in range(nf):
#    if(
#    plt.figure(figsize=(6,6))
#    plt.axes( xlim=(min(x[i])/1.2,max(x[i])*1.2),ylim=(min(y[i])/2.,max(y[i])*5/4.) )
