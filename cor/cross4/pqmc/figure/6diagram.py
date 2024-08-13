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
        #"font.family":'Times New Roman',
        "font.size": 10,
        "mathtext.fontset" :'stix',
        #"font.serif": ['SimSun'],
        }
rcParams.update(config)
#rcParams['mathtext.default'] = 'regular'
#plt.figure(figsize=(7,5))
#plt.figure(num=None, figsize=(2.8,1.7), dpi=300)
#plt.figure(num=None, figsize=(4,3), dpi=300)
fig,axes=plt.subplots(nrows=2, ncols=4, figsize=(8, 6))
plt.subplots_adjust(left=0.2, right=0.85, top=0.9, bottom=0.15, wspace=0.01, hspace=0.1)
ax1=axes[0,0]
ax2=axes[0,1]
ax3=axes[0,2]
ax4=axes[0,3]
ax5=axes[1,0]
ax6=axes[1,1]
ax7=axes[1,2]
fig.delaxes(axes[1,3]) # The indexing is zero-based here
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax4.set_xscale('log')
ax4.set_yscale('log')
ax5.set_xscale('log')
ax5.set_yscale('log')
ax6.set_xscale('log')
ax6.set_yscale('log')
ax7.set_xscale('log')
ax7.set_yscale('log')


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
files = sorted(glob.glob(r'rcoreetype*.dat'))
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
        xx.append(data[j][4])
        yy.append(data[j][1])
        ddy.append(data[j][2])
        tt.append(data[j][5])
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
color11 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange','brown','pink'))
#color11 = itertools.cycle(('red', 'blue', 'green', 'orange','brown','pink'))
label1 = itertools.cycle(('1AA','2AA','3AA','1AB','2AB','3AB','4AB'))
mface  = itertools.cycle(('none','none','none','green','orange','brown','pink'))
#plt.figure(figsize=(6,6))
order=[]
for j in range(nf):
    if(tcor[j][0]==0):
        ax1.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(mface),
                     linestyle='',color=next(color1),label='$Type-%i$' % (tcor[j][0]),capsize=2)
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
        ax1.plot(xf,yf,linestyle='-.',color=next(color11),label='$r^{-1}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
    if(tcor[j][0]==1):
        ax2.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(mface),
                     linestyle='',color=next(color1),label='$Type-%i$' % (tcor[j][0]),capsize=2)
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
        ax2.plot(xf,yf,linestyle='-.',color=next(color11),label='$r^{-1}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
    if(tcor[j][0]==2):
        ax3.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(mface),
                     linestyle='',color=next(color1),label='$Type-%i$' % (tcor[j][0]),capsize=2)
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
        ax3.plot(xf,yf,linestyle='-.',color=next(color11),label='$r^{-1}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
    if(tcor[j][0]==3):
        ax4.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(mface),
                     linestyle='',color=next(color1),label='$Type-%i$' % (tcor[j][0]),capsize=2)
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
        ax4.plot(xf,yf,linestyle='-.',color=next(color11),label='$r^{-1}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
    if(tcor[j][0]==4):
        for kk in range(100):
            #if(x[j][kk] == x[j][2]):
            if(x[j][kk] == x[j][0]):
                AA, BB = curve_fit(power_law, x[j], y[j])[0]
                xf=np.arange((x[j][kk])/1.2,max(x[j])*1.2,2.0)
                AA=x[j][kk]**1*y[j][kk]
                BB=-1
                yf=AA*xf**BB
                break
        ax5.plot(xf,yf,linestyle='-.',color=next(color11),label='$r^{-1}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
        ax5.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(mface),
                     linestyle='',color=next(color1),label='$Type-%i$' % (tcor[j][0]),capsize=2)
    if(tcor[j][0]==5):
        ax6.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(mface),
                     linestyle='',color=next(color1),label='$Type-%i$' % (tcor[j][0]),capsize=2)
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
        ax6.plot(xf,yf,linestyle='-.',color=next(color11),label='$r^{-1}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
    if(tcor[j][0]==6):
        ax7.errorbar(x[j],y[j],yerr=dy[j],marker=next(marker),markerfacecolor=next(mface),
                     linestyle='',color=next(color1),label='$Type-%i$' % (tcor[j][0]),capsize=2)
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
        ax7.plot(xf,yf,linestyle='-.',color=next(color11),label='$r^{-1}$')#,label='y=%5.3fx^%5.3f' % (AA,BB))
    plt.xscale('log')
    plt.yscale('log')
#fitx = np.arange(0,5,0.1)
#plt.plot(fitx,(-fitx/4.)*np.exp((-fitx**2)/2),linestyle='-.',color='black',
         #label='$P(x)=\\frac{x}{4}exp(-x^2/2)$')

ax1.set_xlim([5, 220])
ax1.set_ylim([10**-4.2, 10**-1])
ax2.set_xlim([5, 220])
ax2.set_ylim([10**-4.2, 10**-1])
ax3.set_xlim([5, 220])
ax3.set_ylim([10**-4.2, 10**-1])
ax4.set_xlim([5, 220])
ax4.set_ylim([10**-4.2, 10**-1])
ax5.set_xlim([5, 220])
ax5.set_ylim([10**-4.2, 10**-1])
ax6.set_xlim([5, 220])
ax6.set_ylim([10**-4.2, 10**-1])
ax7.set_xlim([5, 220])
ax7.set_ylim([10**-4.2, 10**-1])
ax1.set_xlabel(r'$r$')
ax1.set_ylabel(r'$C(r)$')
ax2.set_xlabel(r'$r$')
ax2.set_ylabel(r'$C(r)$')
ax3.set_xlabel(r'$r$')
ax3.set_ylabel(r'$C(r)$')
ax4.set_xlabel(r'$r$')
ax4.set_ylabel(r'$C(r)$')
ax5.set_xlabel(r'$r$')
ax5.set_ylabel(r'$C(r)$')
ax6.set_xlabel(r'$r$')
ax6.set_ylabel(r'$C(r)$')
ax7.set_xlabel(r'$r$')
ax7.set_ylabel(r'$C(r)$')
ax1.legend(prop={'size':8},loc=3)
ax2.legend(prop={'size':8},loc=3)
ax3.legend(prop={'size':8},loc=3)
ax4.legend(prop={'size':8},loc=3)
ax5.legend(prop={'size':8},loc=3)
ax6.legend(prop={'size':8},loc=3)
ax7.legend(prop={'size':8},loc=3)
#plt.legend(fontsize='x-large', loc=0, frameon=False, bbox_to_anchor=(0.575, 0.38))
plt.tight_layout() #避免被裁匡
plt.savefig('cross4_eecor.png')
plt.savefig('new_cross4_eecor.pdf',bbox_inches="tight")
plt.show()


            
    
#for i in range(nf):
#    if(
#    plt.figure(figsize=(6,6))
#    plt.axes( xlim=(min(x[i])/1.2,max(x[i])*1.2),ylim=(min(y[i])/2.,max(y[i])*5/4.) )
