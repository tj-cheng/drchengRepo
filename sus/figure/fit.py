import numpy as np
import glob
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define the function to fit
def func(x, c, z):
    return c * x**(-1 + 1/z)

file=[]
nf=0
files = sorted(glob.glob('./p*128locsus.dat'))
for ff in files:
    file.append(ff)
    nf=nf+1
    print(ff)
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
        #if(data[j][0] > 100):
        if(data[j][0] > 100):
            xx.append(1./data[j][0])
            yy.append(data[j][1])
            ddy.append(data[j][2])
    x.append(xx)  #T=1/beta
    y.append(yy)
    dy.append(ddy)
    size[i]=data[j][4]
    prob[i]=data[j][5]

print(x)
print(y)
x = np.array(x)
y = np.array(y)
# Generate some sample data
#x = np.array([1, 2, 3, 4, 5])
#y = np.array([0.5, 0.3, 0.2, 0.15, 0.12])
pfile = np.loadtxt( fname = 'lseed.in' ) #pnum.in write the total number of l seed
print(pfile)
numl=len(pfile)
numl=pfile.size
print('numl',numl)
print(np.shape(pfile))

# Fit the function to the data
#:for i in range(nf):
#:    popt, pcov = curve_fit(func, x[i], y[i])
#:    c_fit = popt[0]
#:    z_fit = popt[1]
#:    print("Fitted values: c = {:.2f}, z = {:.2f}".format(c_fit, z_fit))


# Extract the fitted parameters

zz=[]
pp=[]
for i in range(numl):
    for j in range(nf):
        if(numl == 1):
            exit
        else:
            if(prob[j] == pfile[i]):
                print(prob[j])
                popt, pcov = curve_fit(func, x[j], y[j])
                c_fit = popt[0]
                z_fit = popt[1]
                pp.append(prob[j])
                zz.append(z_fit)
                print("Fitted values: c = {:.2f}, z = {:.2f}".format(c_fit, z_fit))

plt.figure(figsize=(8,4))
plt.title('2-leg ladder(obc) fit $\chi_{loc}=cT^{-1+1/z}$')
plt.plot(pp,zz,linestyle='-',marker='o',color='black')#,label='y=%5.3fx^%5.3f' % (AA,BB))
plt.yscale('log')
#plt.xscale('log')
plt.xlabel('$P_{rm}$')
plt.ylabel('$z$')
plt.savefig('diagram_fit_locsus.pdf')
plt.show()

