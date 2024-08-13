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
#plt.figure(figsize=(5,5))
#plt.figure(num=None, figsize=(2.8,1.7), dpi=300)
#plt.figure(num=None, figsize=(4,3), dpi=300)
#:plt.subplots(nrows=3, ncols=3, figsize=(6, 6))
#:plt.subplots_adjust(left=0.2, right=0.85, top=0.9, bottom=0.15, wspace=0.01, hspace=0.1)
 # 经测试, 发现坐标轴刻度字体大小采用16, Label字体大小采用22, legend大小采用`x-large`, 线宽采用2比较合适, 即使在文章排版后经过缩放也能保证看得清.
#plt.xticks(fontsize=16); plt.yticks(fontsize=16); plt.tick_params(labelsize=16)
#plt.xticks(fontsize=8, fontname = 'Times New Roman')
#plt.yticks(fontsize=8, fontname = 'Times New Roman')
#plt.tick_params(labelsize=8)
 # 修改X轴刻度显示间隔为0.5的倍数(例如0.5, 1.0, 1.5, 2.0, ...)
from matplotlib.pyplot import MultipleLocator
#plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))

ychain = np.arange(1,19)
xchain = np.arange(1,17)
y1=[]
y2=[]
x1=[]
x2=[]
for i in range(0,18):
    y1.append(4)
#plt.plot(y1,ychain,linestyle='-.',color='black')#,label='tes')
for i in range(0,18):
    y2.append(13)
#plt.plot(y2,ychain,linestyle='-.',color='black')#,label='tes')
for i in range(0,16):
    x1.append(5)
#plt.plot(xchain,x1,linestyle='-.',color='black')#,label='tes')
for i in range(0,16):
    x2.append(14)
#plt.plot(xchain,x2,linestyle='-.',color='black')#,label='tes')

def draw_cross(ax):
    ax.plot(y1,ychain,linestyle='-.',color='black')#,label='tes')
    ax.plot(y2,ychain,linestyle='-.',color='black')#,label='tes')
    ax.plot(xchain,x1,linestyle='-.',color='black')#,label='tes')
    ax.plot(xchain,x2,linestyle='-.',color='black')#,label='tes')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')


# 創建2x2的subplot佈局
#fig, axs = plt.subplots(2, 2, figsize=(6, 6))
fig,axes=plt.subplots(nrows=2, ncols=5, figsize=(10,4), sharex=True, sharey=True)
#plt.subplots_adjust(left=0.2, right=0.85, top=0.9, bottom=0.15, wspace=0.01, hspace=0.1)

# 在每個subplot中畫井字號
for ax in axes.flat:
    draw_cross(ax)
    ax.set_box_aspect(1)
    #ax.set(adjustable='box-forced', aspect='equal')

def draw_end(ax,cx1,cy1,name,cc):
    for i in range(len(cx1)):
        ax.scatter(cx1[i],cy1[i],color=cc)
        if(i==1):
            ax.plot(cx1[i],cy1[i],color=cc,label=name)
        else:
            ax.plot(cx1[i],cy1[i],color=cc)
    ax.legend(prop={'size':8})


color1 = itertools.cycle(('black', 'red', 'blue', 'green', 'orange','brown','pink'))
ax1=axes[0,0]
cc1=next(color1)
xx = [[1,4],[1,4],[13,16],[13,16]]
yy = [[14,18],[5,1],[18,14],[1,5]]
draw_end(ax1,xx,yy,"Type-0",cc1)

ax2=axes[0,1]
cc2=next(color1)
xx = [[1,1],[16,16]]
yy = [[5,14],[5,14]]
draw_end(ax2,xx,yy,"Type-1",cc2)
ax2=axes[0,2]
xx = [[1,16],[1,16]]
yy = [[5,5],[14,14]]
draw_end(ax2,xx,yy,"Type-1",cc2)


ax3=axes[0,3]
cc3=next(color1)
xx = [[4,16],[1,13],[4,16],[1,13]]
yy = [[18,14],[14,18],[1,5],[5,1]]
draw_end(ax3,xx,yy,"Type-2",cc3)
ax3=axes[0,4]
xx = [[1,4],[1,4],[13,16],[13,16]]
yy = [[14,1],[5,18],[18,5],[1,14]]
draw_end(ax3,xx,yy,"Type-2",cc3)

ax4=axes[1,0]
cc4=next(color1)
xx = [[4,13],[4,13]]
yy = [[1,1],[18,18]]
draw_end(ax4,xx,yy,"Type-3",cc4)
ax4=axes[1,1]
xx = [[4,4],[13,13]]
yy = [[1,18],[1,18]]
draw_end(ax4,xx,yy,"Type-3",cc4)

ax5=axes[1,2]
cc5=next(color1)
xx = [[1,16],[1,16]]
yy = [[14,5],[5,14]]
draw_end(ax5,xx,yy,"Type-4",cc5)

ax6=axes[1,3]
cc6=next(color1)
xx = [[1,13],[1,13],[4,16],[4,16]]
yy = [[5,18],[14,1],[18,5],[1,14]]
draw_end(ax6,xx,yy,"Type-5",cc6)

ax7=axes[1,4]
cc7=next(color1)
xx = [[4,13],[4,13]]
yy = [[18,1],[1,18]]
draw_end(ax7,xx,yy,"Type-6",cc7)
#for i in range(len(cx1)):
#    if(i==1):
#        ax1.plot(cx1[i],cy1[i],color='r',label='Type-0')
#    else:
#        ax1.plot(cx1[i],cy1[i],color='r')
#ax1.legend(prop={'size':8})
#
#ax2=axes[0,1]
#cx2 = [[1,1],[16,16]]
#cy2 = [[5,14],[5,14]]
#for i in range(len(cx2)):
#    if(i==1):
#        ax2.plot(cx2[i],cy2[i],color='r',label='Type-1')
#    else:
#        ax2.plot(cx2[i],cy2[i],color='r')
#ax2.legend(prop={'size':8})

#plt.legend(fontsize='x-large', loc=0, frameon=False, bbox_to_anchor=(0.575, 0.38))
plt.tight_layout() #避免被裁匡
plt.savefig('cross4lattice.png')
plt.savefig('cross4lattice.pdf',bbox_inches="tight")
plt.show()


            
    
#for i in range(nf):
#    if(
#    plt.figure(figsize=(6,6))
#    plt.axes( xlim=(min(x[i])/1.2,max(x[i])*1.2),ylim=(min(y[i])/2.,max(y[i])*5/4.) )
