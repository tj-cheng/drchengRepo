import numpy as np
import matplotlib.pyplot as plt

fig=plt.figure()

#big axes
x=[1,2,3,4,5,6,7]
y=[1,3,4,2,9,5,6]
left,bottom,width,height=0.1,0.1,0.8,0.8   #表示从整张图的左边10%，下边10%的位置开始，宽度、高度分别是是整个图的80%
ax1=fig.add_axes([left,bottom,width,height])
ax1.plot(x,y,'r')   #'r'表示红色的线
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title("big axes")

#在左上角画小图，y值是大图的2倍
x1=x
y1=[]
for a in y:
    b=a*2
    y1.append(b)
left,bottom,width,height=0.2,0.6,0.2,0.2   #表示从整张图的左边10%，下边10%的位置开始，宽度、高度分别是是整个图的80%
ax2=fig.add_axes([left,bottom,width,height])
ax2.plot(x1,y1,'b')
ax2.set_title('little1')

plt.show()
