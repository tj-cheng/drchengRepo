import glob
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.optimize import curve_fit

#plt.rc('font',family='Times New Roman')
 # 修改公式中默认字体
from matplotlib import rcParams
config = {
        "font.family":'serif',
        "font.size": 10,
        "mathtext.fontset" :'stix',
        #"font.serif": ['SimSun'],
        }
rcParams.update(config)
#rcParams['mathtext.default'] = 'regular'
#plt.figure(figsize=(7,5))
#plt.figure(num=None, figsize=(2.8,1.7), dpi=300)
plt.figure(num=None, figsize=(4,3), dpi=300)
#:plt.subplots(nrows=3, ncols=3, figsize=(6, 6))
#:plt.subplots_adjust(left=0.2, right=0.85, top=0.9, bottom=0.15, wspace=0.01, hspace=0.1)
 # 经测试, 发现坐标轴刻度字体大小采用16, Label字体大小采用22, legend大小采用`x-large`, 线宽采用2比较合适, 即使在文章排版后经过缩放也能保证看得清.
#plt.xticks(fontsize=16); plt.yticks(fontsize=16); plt.tick_params(labelsize=16)
plt.xticks(fontsize=8, fontname = 'Times New Roman')
plt.yticks(fontsize=8, fontname = 'Times New Roman')
plt.tick_params(labelsize=8)
 # 修改X轴刻度显示间隔为0.5的倍数(例如0.5, 1.0, 1.5, 2.0, ...)
from matplotlib.pyplot import MultipleLocator
plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))

def power_law(xi,A,B):
    return A*xi**B
file=[]
nf=0
#files = sorted(glob.glob(r'hL*_cor.dat'))
files = sorted(glob.glob(r'rcorpL*.dat'))
for ff in files:
    file.append(ff)
    nf=nf+1
alldata=[]
for i in range(nf):
    data = np.loadtxt(fname = file[i])
    alldata.append(data)
#print(alldata)
print(file)

x=[]
y=[]
dy=[]
prob=np.zeros(nf)
size=np.zeros(nf)
tcor=[]
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
    #size[i]=data[j][2]
#    prob[i]=data[j][3]
#print(alldata[1])
#print(x[1],size[1],prob[1])


#---find same prob---
#Lfile=np.loadtxt( fname= 'lseed.in' )
#pfile = np.loadtxt( fname = 'pseed.in' )
#np=len(pfile)
#nump=pfile.size
marker = itertools.cycle(('p', 'o', 's', 'D', '*','p','^','h'))
color1 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange','brown','pink'))
color11 = itertools.cycle(('red', 'blue', 'green', 'orange','brown','pink'))
label1 = itertools.cycle(('1AA','2AA','3AA','1AB','2AB','3AB','4AB'))
mface  = itertools.cycle(('none','none','none','green','orange','brown','pink'))
#plt.figure(figsize=(6,6))
order=[]
for j in range(nf):
    plt.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(mface),
                 linestyle='',color=next(color1),label='$P_{rm}=1-\\frac{20}{L}$',capsize=2,markersize=2)
                 #linestyle='',color=next(color1),label='$Type-%i$' % (tcor[j][0]),capsize=2)
    for kk in range(100):
        if(x[j][kk] == x[j][0]):
            #print(1)
            #print(min(x[j]),max(x[j]))
            AA, BB = curve_fit(power_law, x[j], y[j])[0]
            xf=np.arange((x[j][kk])/1.2,max(x[j])*1.2,2.0)
            AA=x[j][kk]**1*y[j][kk]
            BB=-1
            yf=AA*xf**BB
            break
    plt.plot(xf,yf,linestyle='-.',color=next(color11),label='$L^{-1}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
    for a, b, c in zip(x[j], y[j], tcor[j]):
        plt.text(a, 1.2*b, '$P_{rm}=$%.3f' % c, ha='center', va='bottom', fontsize=6, color='blue')
        #for kk in range(100):
        #    if(x[j][kk] == 68):
        #        #print(1)
        #        #print(min(x[j]),max(x[j]))
        #        AA, BB = curve_fit(power_law, x[j], y[j])[0]
        #        xf=np.arange((x[j][kk])/1.2,max(x[j])*1.2,2.0)
        #        AA=x[j][kk]**2*y[j][kk]
        #        BB=-2
        #        yf=AA*xf**BB
        #        break
        #plt.plot(xf,yf,linestyle='-.',color=next(color11),label='$r^{-2}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
    #for kk in range(100):
    #    if(x[j][kk] == 64):
    #        #print(1)
    #        #print(min(x[j]),max(x[j]))
    #        AA, BB = curve_fit(power_law, x[j], y[j])[0]
    #        xf=np.arange((x[j][kk])/1.2,max(x[j])*1.2,2.0)
    #        AA=x[j][kk]**2*y[j][kk]
    #        BB=-2
    #        yf=AA*xf**BB
    #        break
    #plt.plot(xf,yf,linestyle='-.',color=next(color11),label='$L^{-2}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
    plt.xscale('log')
    plt.yscale('log')
#fitx = np.arange(0,5,0.1)
#plt.plot(fitx,(-fitx/4.)*np.exp((-fitx**2)/2),linestyle='-.',color='black',
         #label='$P(x)=\\frac{x}{4}exp(-x^2/2)$')
#D=0.4598
#A=0.629
#B=1./0.1128
#plt.plot(fitx,(fitx*A)*np.exp(-(fitx**2)/B-D*(fitx)),linestyle='-.',color='black')#,
#          label='$P(x)=%1.1f*|x|exp(-|x|^2/%1.1f-%1.3f|x|)$' % (A,B,D))
#fitx = np.arange(50,150,0.1)
#D=0.7
##A=10**(-2)*60.
#A=11**(-2)*60.
#plt.plot(fitx,A*fitx**(-1),linestyle='-.',color='black',label='$L^{-1}$')
#A=1**(-2)*40
#plt.plot(fitx,A*fitx**(-2),linestyle='-.',color='red',label='$L^{-2}$')
#plt.title('cross4,random J')
#plt.xlabel('$r$', fontdict={'family' : 'Times New Roman', 'size':8})
plt.xlabel(r'$L$')
plt.ylabel(r'$C(L)$')
plt.legend(prop={'size':8})
#plt.legend(fontsize='x-large', loc=0, frameon=False, bbox_to_anchor=(0.575, 0.38))
plt.tight_layout() #避免被裁匡
plt.savefig('pair_cor.png')
plt.savefig('pair_cor.pdf',bbox_inches="tight")
plt.show()


            
    
#for i in range(nf):
#    if(
#    plt.figure(figsize=(6,6))
#    plt.axes( xlim=(min(x[i])/1.2,max(x[i])*1.2),ylim=(min(y[i])/2.,max(y[i])*5/4.) )
