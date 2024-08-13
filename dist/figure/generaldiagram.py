import glob
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.optimize import curve_fit

def power_law(xi,A,B):
    return A*xi**B

# Define the function to fit
def func(x, c, z):
    return c * x**(-1 + 1/z)
#files = sorted(glob.glob('p*f*_cor.dat'))
#files = sorted(glob.glob('./p*128locsus.dat'))
def data_to_array(files,file,nf):
    file=[]
    nf=0
    for ff in files:
        file.append(ff)
        nf=nf+1
        print(ff)
    alldata=[]
    for i in range(nf):
        data = np.loadtxt(fname = file[i])
        alldata.append(data)
    #print(alldata)
    return data,file,nf


def array_to_xy(file,nf,x,y,dy,prob,size):
    x=[]
    y=[]
    dy=[]
    prob=np.zeros(nf)
    size=np.zeros(nf)
    for i in range(nf):
        data = np.loadtxt(fname = file[i])
        xx=[]; yy=[]; ddy=[]
        for j in range(len(data)):
            xx.append(data[j][5])
            yy.append(data[j][1])
        x.append(xx)  #T=1/beta
        y.append(yy)
        size[i]=data[j][0]
        prob[i]=data[j][6]
    return x,y,dy,size,prob



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
#markerface = itertools.cycle(('black', 'red', 'blue', 'green', 'orange','purple','olive','tomato','deepskyblue','orchid','none','none','none','none','none','none'))
markerface = itertools.cycle(('none','none','none','none','none','none'))
#color2 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange'))
#plt.errorbar(x[0],y[0],yerr=dy[0],marker=next(marker),markerfacecolor=next(markerface), label='16' ,linestyle='-',color=next(color1),capsize=2)
#               #plt.plot(xf,yf,linestyle='-.',color=next(color2))#,label='y=%5.3fx^%5.3f' % (AA,BB))
#plt.yscale('log')
#plt.xscale('log')
#print(prob)
def orderdraw(pfile,size,prob,x,y,dy,legendname):
    for i in range(numl):
        for j in range(nf):
            if(numl == 999999):
                break
            else:
                if(prob[j] == pfile[i]):
                    #plt.errorbar(x[j],y[j],yerr=0,marker=next(marker),markerfacecolor=next(markerface), label=legendname % (prob[j]),linestyle='-',color=next(color1),capsize=0)
                    plt.errorbar(x[j],y[j],yerr=0,marker='',markerfacecolor=next(markerface), label=legendname % (prob[j]),linestyle='-',color=next(color1),capsize=0)

                    #plt.yscale('log')
                    #plt.xscale('log')
    return

#-------------part 1----------------#
files = sorted(glob.glob('./ladderevol*.dat'))
file=[]
nf=0
data,file,nf = data_to_array(files,file,nf)
x=[]
y=[]
dy=[]
prob=np.zeros(nf)
size=np.zeros(nf)
x,y,dy,size,prob = array_to_xy(file,nf,x,y,dy,prob,size)
legendname='$P_{rm}=$%1.2f'
orderdraw(pfile,size,prob,x,y,dy,legendname)
#plt.xlim(280,360)
#-----------------------------------#
##-------------part 2----------------#
#files = sorted(glob.glob('./cleanp0128locsus.dat'))
#file=[]
#nf=0
#data,file,nf = data_to_array(files,file,nf)
#x=[]
#y=[]
#dy=[]
#prob=np.zeros(nf)
#size=np.zeros(nf)
#x,y,dy,size,prob = array_to_xy(file,nf,x,y,dy,prob,size)
##legendname='$clean \ P_{rm}=$%1.1f'
#legendname='$clean \ P_{rm}=$%1.1f ; $z_{pure}=$%1.3f'
#orderdraw(pfile,size,prob,x,y,dy,legendname)
#-----------------------------------#
#fitx = np.arange(60,130,0.1)
#A=10**(-2)*60.
#plt.plot(fitx,A*fitx**(-1),linestyle='-.',color='black',label='$L^{-1}$')
plt.title('2-leg ladder(obc) $\chi_{loc} \propto T^{-1+1/z}$, L=128, uniform disorder')
#plt.xlabel(chr(946))
plt.xlabel('$T$')
plt.ylabel('$\chi_{loc}$')
plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)
plt.tight_layout() #避免被裁匡
#plt.savefig('newlocsysdiagram.png')
#plt.savefig('newlocsysdiagram.pdf')
plt.show()
            
    
#for i in range(nf):
#    if(
#    plt.figure(figsize=(6,6))
#    plt.axes( xlim=(min(x[i])/1.2,max(x[i])*1.2),ylim=(min(y[i])/2.,max(y[i])*5/4.) )
