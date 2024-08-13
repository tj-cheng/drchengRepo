import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#data8 = np.loadtxt(fname = 'tsor8.dat')
#data4 = np.loadtxt(fname = 'tsor4.dat')
data1 = np.loadtxt(fname = 'cm1.dat')
data2 = np.loadtxt(fname = 'cm2.dat')
#data3 = np.loadtxt(fname = 'cx.dat')
#data4 = np.loadtxt(fname = 'cy.dat')
data3 = np.loadtxt(fname = 'cm3.dat')
data4 = np.loadtxt(fname = 'cm4.dat')
data5 = np.loadtxt(fname = 'cp1.dat')
data6 = np.loadtxt(fname = 'cp2.dat')
data7 = np.loadtxt(fname = 'cp3.dat')
x1=[]; x2=[]; x3=[]; x4=[]; x5=[]; x6=[]; x7=[]
y1=[]; y2=[]; y3=[]; y4=[]; y5=[]; y6=[]; y7=[]
dy1=[]; dy2=[]; dy3=[]; dy4=[]; dy5=[]; dy6=[]; dy7=[]
def power_law(xi,A,B):
    return A*xi**B
for i in range(len(data1)):
    x1.append(data1[i][0])
    y1.append(data1[i][1])
    dy1.append(data1[i][2])
    x2.append(data2[i][0])
    y2.append(data2[i][1])
    dy2.append(data2[i][2])
    x3.append(data3[i][0])
    y3.append(data3[i][1])
    dy3.append(data3[i][2])
    x4.append(data4[i][0])
    y4.append(data4[i][1])
    dy4.append(data4[i][2])
    x5.append(data5[i][0])
    y5.append(data5[i][1])
    dy5.append(data5[i][2])
    x6.append(data6[i][0])
    y6.append(data6[i][1])
    dy6.append(data6[i][2])
    x7.append(data7[i][0])
    y7.append(data7[i][1])
    dy7.append(data7[i][2])


A1, B1 = curve_fit(power_law, x1, y1)[0]
A2, B2 = curve_fit(power_law, x2, y2)[0]
A3, B3 = curve_fit(power_law, x3, y3)[0]
A4, B4 = curve_fit(power_law, x4, y4)[0]
A5, B5 = curve_fit(power_law, x5, y5)[0]
A6, B6 = curve_fit(power_law, x6, y6)[0]
A7, B7 = curve_fit(power_law, x7, y7)[0]
#xx=np.arange(min(Ndat)-10,max(Ndat)+10,1.0)

xf1=np.arange(min(x1)/1.2,max(x1)*1.2,2.0)
yf1=A1*xf1**B1
xf2=np.arange(min(x2)/1.2,max(x2)*1.2,2.0)
yf2=A2*xf2**B2
xf3=np.arange(min(x3)/1.2,max(x3)*1.2,2.0)
yf3=A3*xf3**B3
xf4=np.arange(min(x4)/1.2,max(x4)*1.2,2.0)
yf4=A4*xf4**B4
xf5=np.arange(min(x5)/1.2,max(x5)*1.2,2.0)
yf5=A5*xf5**B5
xf6=np.arange(min(x6)/1.2,max(x6)*1.2,2.0)
yf6=A6*xf6**B6
xf7=np.arange(min(x7)/1.2,max(x7)*1.2,2.0)
yf7=A7*xf7**B7

plt.figure(figsize=(6,6))


ax=plt.axes( xlim=(min(x1)/1.2,max(x1)*1.2),ylim=(min(y7)/10,max(y1)*10) )
p1 = plt.errorbar(x1,y1,yerr=dy1,fmt='o',ecolor='black',
        elinewidth=1,ms=5,mfc='black',mec='black', capsize=3)#, label='1AB')
p2 = plt.errorbar(x2,y2,yerr=dy2,fmt='v',ecolor='red',
        elinewidth=1,ms=5,mfc='red',mec='red', capsize=3)#, label='2AB ')
p3 = plt.errorbar(x3,y3,yerr=dy3,fmt='s',ecolor='blue',
        elinewidth=1,ms=5,mfc='blue',mec='blue', capsize=3)#, label='3AB ')
p4 = plt.errorbar(x4,y4,yerr=dy4,fmt='d',ecolor='green',
        elinewidth=1,ms=5,mfc='green',mec='green', capsize=3)#, label='4AB ')
p5 = plt.errorbar(x5,y5,yerr=dy5,fmt='^',ecolor='sienna',
        elinewidth=1,ms=5,mfc='sienna',mec='sienna', capsize=3)#, label='1AA ')
p6 = plt.errorbar(x6,y6,yerr=dy6,fmt='<',ecolor='purple',
        elinewidth=1,ms=5,mfc='purple',mec='purple', capsize=3)#, label='2AA ')
p7 = plt.errorbar(x7,y7,yerr=dy7,fmt='>',ecolor='tomato',
        elinewidth=1,ms=5,mfc='tomato',mec='tomato', capsize=3)#, label='3AA ')
l1=plt.legend([p1,p2,p3,p4,p5,p6,p7],['1AB','2AB',
    '3AB','4AB','1AA','2AA','3AA'],loc='upper right')
#salmon
r1 = ax.plot(xf1,yf1,color='black',linestyle='-.',label='y=%5.3fx^%5.3f' % (A1,B1))
r2 = ax.plot(xf2,yf2,color='red',linestyle='-.',label='y=%5.3fx^%5.3f' % (A2,B2))
r3 = ax.plot(xf3,yf3,color='blue',linestyle='-.',label='y=%5.3fx^%5.3f' % (A3,B3))
r4 = ax.plot(xf4,yf4,color='green',linestyle='-.',label='y=%5.3fx^%5.3f' % (A4,B4))
r5 = ax.plot(xf5,yf5,color='sienna',linestyle='-.',label='y=%5.3fx^%5.3f' % (A5,B5))
r6 = ax.plot(xf6,yf6,color='purple',linestyle='-.',label='y=%5.3fx^%5.3f' % (A6,B6))
r7 = ax.plot(xf7,yf7,color='tomato',linestyle='-.',label='y=%5.3fx^%5.3f' % (A7,B7))
plt.legend(loc='lower left')
plt.gca().add_artist(l1)

plt.xlabel("L")
plt.ylabel("C(L)")
plt.title('Spin correlations (crossing 4)')
plt.yscale('log')
plt.xscale('log')
#plt.xscale('log')
#plt.legend()
plt.savefig('tjcrossing4.pdf')
plt.show()
print('plot is done')
