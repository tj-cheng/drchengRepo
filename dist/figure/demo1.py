import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import glob
from matplotlib import rcParams

config = {
        "font.family":'Times New Roman',
        "font.size": 10,
        "mathtext.fontset" :'stix',
        #"font.serif": ['SimSun'],
        }
rcParams.update(config)
# 示例数据
#data = [
#    [1, [(30, 0.5), (60, 0.8), (90, 0.2)]],
#    [2, [(30, 0.7), (60, 0.4), (90, 0.9)]],
#    [3, [(30, 0.2), (60, 0.6), (90, 0.3)]],
#]
files = sorted(glob.glob(r'p*dd100*.dat'), key=lambda x: int(x.split('p')[1].split('dd')[0]))
files2 = sorted(glob.glob(r'p*dd150*.dat'), key=lambda x: int(x.split('p')[1].split('dd')[0]))
files3 = sorted(glob.glob(r'p*dd200*.dat'), key=lambda x: int(x.split('p')[1].split('dd')[0]))
for i in range(len(files)):
    print(i, files[i])
aa=0.5
bb=0.55
#print(files)
#file_path = 'p0dd100evo.dat'  # 替换为您的文件路径

data = []
wspin = []
cc='cubehelix'
#cc='twilight_shifted'
#cc='cividis'

def readdata(data,open_path,wspin):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        angle_values = [int(val.split()[0]) for val in lines]
        #print(angle_values)
        wspin_values = [float(val.split()[1]) for val in lines]
        #print(wspin_values)
        radius_values = [float(val.split()[2]) for val in lines]
        #print(radius_values)
        angle_value_pairs = [(angle, wspin_values[j]) for j, angle in enumerate(range(0, 360))]  # 修正此处的索引
        data.append([(1-radius_values[0])*1, angle_value_pairs])
        wspin.append(wspin_values)
    return data,wspin

fig, axes = plt.subplots(nrows=1, ncols=3, subplot_kw=dict(projection='polar'))
ax1=axes[0]
ax2=axes[1]
ax3=axes[2]

for i in range(len(files)):
    file_path=files[i]
#file_path = 'p0dd100evo.dat'  # 替换为您的文件路径
    readdata(data,file_path,wspin)
#print(wspin)
#print(np.max(wspin,axis=1))
#print(np.max(wspin))
#print(data)
# 创建一个图形
#fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# 提取数据
radii = [item[0] for item in data]
#radii = [item[2] for item in data]
angle_values = [item[1] for item in data]

print('ra=',radii)
# 获取最大的角度数
max_angle = max(max(angle_values, key=lambda x: len(x)), key=lambda x: x[0])[0]

# 创建一个矩阵，用于存储值
matrix = np.zeros((len(radii), max_angle + 1))

# 将数据填充到矩阵中
for i, (_, values) in enumerate(data):
    for angle, value in values:
        matrix[i, angle] = value
vmin,vmax = np.min(wspin),np.max(wspin)
#print(radii)
#print(matrix[18])
#print(np.arange(0,max_angle+1,1))
#c1 = ax1.pcolormesh(np.radians(np.arange(0, max_angle + 1, 1)), radii, matrix, cmap=cc,vmin=aa,vmax=bb,shading='auto')
c1 = ax1.pcolormesh(np.radians(np.arange(0, max_angle + 1, 1)), radii, matrix, cmap=cc,vmin=aa,vmax=bb)
########################################
data=[]
wspin = []
for i in range(len(files2)):
    file_path=files2[i]
    readdata(data,file_path,wspin)
radii = [item[0] for item in data]
angle_values = [item[1] for item in data]
max_angle = max(max(angle_values, key=lambda x: len(x)), key=lambda x: x[0])[0]
matrix = np.zeros((len(radii), max_angle + 1))
for i, (_, values) in enumerate(data):
    for angle, value in values:
        matrix[i, angle] = value
vmin,vmax = np.min(wspin),np.max(wspin)
c2 = ax2.pcolormesh(np.radians(np.arange(0, max_angle + 1, 1)), radii, matrix, cmap=cc,vmin=aa,vmax=bb)
########################################
########################################
data=[]
wspin = []
for i in range(len(files3)):
    file_path=files3[i]
    readdata(data,file_path,wspin)
radii = [item[0] for item in data]
angle_values = [item[1] for item in data]
max_angle = max(max(angle_values, key=lambda x: len(x)), key=lambda x: x[0])[0]
matrix = np.zeros((len(radii), max_angle + 1))
for i, (_, values) in enumerate(data):
    for angle, value in values:
        matrix[i, angle] = value
vmin,vmax = np.min(wspin),np.max(wspin)
c3 = ax3.pcolormesh(np.radians(np.arange(0, max_angle + 1, 1)), radii, matrix, cmap=cc,vmin=aa,vmax=bb)
########################################
#
px=0*3.14/180
tt=[]
tx=[]
for i in range(11):
    tt.append(i*0.1)
    tx.append(px)


#ax1.set_rticks(tt)
ax1.set_rticks([])
ax2.set_rticks([])
ax3.set_rticks([])
#ax3.set_rticks([0,0.5,1])
ax1.text(20*3.14/180, 0.5, '$r=1-P_{rm}$', ha='center', va='center',color='pink',size=12)
ax1.plot(tx, tt, 'x', markersize=2, color='pink')
ax2.text(20*3.14/180, 0.5, '$r=1-P_{rm}$', ha='center', va='center',color='pink',size=12)
ax2.plot(tx, tt, 'x', markersize=2, color='pink')
ax3.text(20*3.14/180, 0.5, '$r=1-P_{rm}$', ha='center', va='center',color='pink',size=12)
ax3.plot(tx, tt, 'x', markersize=2, color='pink')
#ax2.plot(320*3.14/180, 0.3+0.01, 'x', markersize=5, color='blue')
#ax3.plot(320*3.14/180, 0.3+0.01, 'x', markersize=5, color='blue')
# 添加颜色条
# 设置箭头的弧度和长度
angles = np.radians([0, 90, 180, 270])
ax1.set_xticks(angles)
ax2.set_xticks(angles)
ax3.set_xticks(angles)
# 绘制箭头
#ax1.arrow(0, 0, angle, radius, head_width=0.1, head_length=0.1, fc='red', ec='red')

#ax.clim=(0.3,0.5)
cbar = plt.colorbar(c1, ax=ax1, orientation='horizontal')
cbar = plt.colorbar(c2, ax=ax2, orientation='horizontal')
cbar = plt.colorbar(c3, ax=ax3, orientation='horizontal')
#cbar = plt.colorbar(c, ax=ax, ticks=np.linspace(0.3, 0.5, 10))
# 添加标题
ax1.set_title(r'$dd=1$')
ax2.set_title(r'$dd=1.5$')
ax3.set_title(r'$dd=2.0$')
plt.tight_layout() #避免被裁匡

#plt.savefig('spareladder_nop_eecor.pdf',bbox_inches="tight")
# 显示图形
plt.show()

