import glob
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.optimize import curve_fit
from matplotlib import rcParams
#config = {
#        #"font.family":'Times New Roman',
#        "font.family":'serif',
#        "font.size": 10,
#        "mathtext.fontset" :'stix',
#        #"font.serif": ['SimSun'],
#        }
#rcParams.update(config)

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
            if(data[j][0] > 100):
                xx.append(1./data[j][0])
                yy.append(data[j][1])
                ddy.append(data[j][2])
        x.append(xx)  #T=1/beta
        y.append(yy)
        dy.append(ddy)
        size[i]=data[j][4]
        prob[i]=data[j][5]
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
color11 = itertools.cycle(('green','black', 'red', 'blue', 'green', 'orange','purple','olive','tomato','deepskyblue','orchid'))
markerface = itertools.cycle(('black', 'red', 'blue', 'green', 'orange','purple','olive','tomato','deepskyblue','orchid','none','none','none','none','none','none'))
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
                if(size[j] == pfile):
                    print(size[j])
                    for kk in range(100):
                        if(x[j][kk] == 16):
                            AA, BB = curve_fit(power_law, x[j], y[j])[0]
                            xf=np.arange((x[j][kk])/1.2,max(x[j])*1.2,2.0)
                            AA=x[j][kk]*y[j][kk]
                            BB=-1
                            yf=AA*xf**BB
                            break
            else:
                if(prob[j] == pfile[i]):
                    #print(prob[j])
                    popt, pcov = curve_fit(func, x[j], y[j])
                    c_fit = popt[0]
                    z_fit = popt[1]
                    for kk in range(100):
                        #if(x[j][kk] == 1./512.):
                        #    AA, BB = curve_fit(power_law, x[j], y[j])[0]
                        #    xf=np.arange((x[j][kk])/5.,max(x[j])*1.2,0.0001)
                        #    AA=x[j][kk]*y[j][kk]
                        #    BB=-1
                        #    yf=AA*xf**BB
                        #    break
                        if(x[j][kk] == 1./1024.):
                             temp=x[j][kk]
                             xaa=y[j][kk]
                             conf=xaa*temp*(np.log(temp))**2
                             xf=np.arange((x[j][kk])/5.,max(x[j])*1.2,0.0001)
                             yf=conf/(xf*(np.log(xf)**2))
                             print(conf)
                             break
                    #if(prob[j] > 0.7 and prob[j] < 1.1):
                    plt.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(markerface), label=legendname % (prob[j],z_fit),linestyle='-',color=next(color1),capsize=2)
                    if(prob[j] > 0.9):
                        plt.plot(xf,yf,linestyle='-.',color=next(color11),label='$\\frac{1}{Tln^2|T|}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
                    #elif(prob[j] > 0.6):
                    #    plt.plot(xf,yf,linestyle='-.',color=next(color11))#,label='$\\frac{1}{Tln^2|T|}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
#                        plt.plot(xf,yf,linestyle='-.',color=next(color11),label='$\\frac{1}{Tln^2|T|}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
                        #for kk in range(100):
                        #    if(x[j][kk] == 1./128.):
                        #         temp=1./128.
                        #         #xaa=y[j][kk]
                        #         #conf=xaa*temp*(np.log(temp))**2
                        #         conf=y[j][kk]*x[j][kk]
                        #         xf=np.arange((x[j][kk])/5.,max(x[j])*1.2,0.0001)
                        #         yf=conf/xf
                        #         print(conf)
                        #         break
                        #plt.plot(xf,yf,linestyle='-.',color='green',label='$\\frac{1}{T}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
                    #else:
                    #    plt.plot(xf,yf,linestyle='-.',color=next(color11))#,label='y=%5.3fx^%5.3f' % (AA,BB))
                    
                    #plt.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(markerface), label=legendname % (prob[j],z_fit),linestyle='-',color=next(color1),capsize=2)

                    plt.yscale('log')
                    plt.xscale('log')
    return

def noporderdraw(pfile,size,x,y,dy,legendname):
    j=0
    popt, pcov = curve_fit(func, x[j], y[j])
    c_fit = popt[0]
    z_fit = popt[1]
    for kk in range(100):
        if(x[j][kk] == 1./256.):
            temp=1./256.
            xaa=y[j][kk]
            conf=xaa*temp*(np.log(temp))**2
            xf=np.arange((x[j][kk])/5.,max(x[j])*1.2,0.0001)
            yf=conf/(xf*(np.log(xf)**2))
            print(conf)
            break
    plt.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(markerface), label=legendname % (z_fit),linestyle='-',color=next(color1),capsize=2)

    plt.yscale('log')
    plt.xscale('log')
    return
#-------------part 1----------------#
files = sorted(glob.glob('./p*128unsus.dat'))
file=[]
nf=0
data,file,nf = data_to_array(files,file,nf)
x=[]
y=[]
dy=[]
prob=np.zeros(nf)
size=np.zeros(nf)
x,y,dy,size,prob = array_to_xy(file,nf,x,y,dy,prob,size)
legendname='$P_{rm}=$%1.1f ; $z\'=$%1.3f'
orderdraw(pfile,size,prob,x,y,dy,legendname)
#-----------------------------------#
#-------------part 2----------------#
#-------------part 2----------------#
files = sorted(glob.glob('./cleanun*.dat'))
file=[]
nf=0
data,file,nf = data_to_array(files,file,nf)
x=[]
y=[]
dy=[]
prob=np.zeros(nf)
size=np.zeros(nf)
x,y,dy,size,prob = array_to_xy(file,nf,x,y,dy,prob,size)
legendname='$clean \ P_{rm}=$%1.1f'
legendname='clean disorder connectivities ; $z\'=$%1.3f'
noporderdraw(pfile,size,x,y,dy,legendname)
#-----------------------------------#
files = sorted(glob.glob('./cross4un*.dat'))
file=[]
nf=0
data,file,nf = data_to_array(files,file,nf)
x=[]
y=[]
dy=[]
prob=np.zeros(nf)
size=np.zeros(nf)
x,y,dy,size,prob = array_to_xy(file,nf,x,y,dy,prob,size)
#legendname='$clean \ P_{rm}=$%1.1f'
legendname='$cross-$%i ; $z\'=$%1.3f'
orderdraw(pfile,size,prob,x,y,dy,legendname)
#-----------------------------------#
#fitx = np.arange(60,130,0.1)
#A=10**(-2)*60.
#plt.plot(fitx,A*fitx**(-1),linestyle='-.',color='black',label='$L^{-1}$')
#plt.title('2-leg ladder(obc) $\chi_{loc} \propto T^{-1+1/z}$, L=128, uniform disorder')
#plt.xlabel(chr(946))
plt.xlabel(r'$T$')
plt.ylabel(r'$\chi_{u}$')
plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)
plt.tight_layout() #避免被裁匡
plt.savefig('newunsus.png')
plt.savefig('newunsus.pdf')
plt.show()
            
    
#for i in range(nf):
#    if(
#    plt.figure(figsize=(6,6))
#    plt.axes( xlim=(min(x[i])/1.2,max(x[i])*1.2),ylim=(min(y[i])/2.,max(y[i])*5/4.) )
