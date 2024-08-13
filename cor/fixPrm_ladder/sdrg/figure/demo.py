import matplotlib.pyplot as plt
import numpy as np
import random
#ax = plt.subplot(1,1,1)
#ax = plt.subplot()
#plt.annotate("",
#            xy=(0,0),
#            xytext=(1,0),
#            size=20,va="center",ha="center",
#            arrowprops=dict(color='black',
#                            arrowstyle="-",
#                            connectionstyle="arc3,rad=0.4",
#                            )
#            )
#datajj = np.loadtxt(fname = 'L20p0f90bondcoup.dat')
rsing = np.loadtxt(fname = 'psinglet.dat')
rdsing = np.loadtxt(fname = 'pdash.dat')
#lsing = np.loadtxt(fname = 'L20p0f90lsinglet.dat')
#lspin = np.loadtxt(fname = 'Ll20p0f90spin.dat')
rspin = np.loadtxt(fname = 'pspin.dat')
spos = np.loadtxt(fname = 'ppos.dat')
#plt.figure(figsize=(np.int32(len(spos)/2),2))
plt.figure(figsize=(16,8))
jj=[]
xb=[];yb=[]
xir=[];xfr=[];yir=[];yfr=[]
xil=[];xfl=[];yil=[];yfl=[]
xid=[];xfd=[];yid=[];yfd=[]
#for i in range(len(datajj)):
#    xb.append([datajj[i][0],datajj[i][2]])
#    yb.append([datajj[i][1],datajj[i][3]])
#    jj.append(datajj[i][4])
for i in range(len(rsing)):
    xir.append(rsing[i][0])
    yir.append(rsing[i][1])
    xfr.append(rsing[i][2])
    yfr.append(rsing[i][3])
    #xil.append(lsing[i][0])
    #yil.append(lsing[i][1])
    #xfl.append(lsing[i][2])
    #yfl.append(lsing[i][3])
for i in range(len(rdsing)):
    xid.append(rdsing[i][0])
    yid.append(rdsing[i][1])
    xfd.append(rdsing[i][2])
    yfd.append(rdsing[i][3])
#jjsum=np.sum(jj)
#jj=jj/jjsum
#jjmax=np.max(jj)
#jj=jj/jjmax
#----------------I set lspn is the same as rspn-----------#
#spi=[1,1,1,1,1,1,1,1]
#xi=[0,1,2,4,5,6,8,9]
#yi=[0,0,0,0,0,0,1,1]
#xf=[3,11,10,15,7,14,12,13]
#yf=[0,1,1,1,0,1,1,1]
mk=5
def getrad(xi,xf,spos):
    radd=[]
    for i in range(len(xi)):
        delta=np.abs(xi[i]-xf[i])
        #radd.append(-0.05*delta+0.45)
        #radd.append(0.4*0.5**(delta/5.))
        radd.append(0.4*0.5**((delta-1)/len(xi)))
    return radd
def getrad_dash(xi,xf,spos):
    radd_dash=[]
    for i in range(len(xi)):
        delta=np.abs(xi[i]-xf[i])
        #radd.append(-0.05*delta+0.45)
        #radd.append(0.4*0.5**(delta/5.))
        radd_dash.append(0.5*0.3**((delta-1)/len(xi)))
    return radd_dash
def drawpicture1(ab,bb,spi,ai,bi,af,bf,adi,bdi,adf,bdf,radd,radd_dash,y0,pos):
    #for i in range(len(ab)):
    #    for j in range(2):
    #        bb[i][j]=bb[i][j]+y0   #跟排列vr跟vl順序而變
    #for i in range(len(ab)):
    #    plt.plot(ab[i],bb[i],color="royalblue",linestyle="-",linewidth=str(3*jj[i]+1))
    
    for i in range(len(spi)):
        if( spi[i] > 0):
            plt.plot(pos[i][0],pos[i][1]+y0,color="black",marker='o',markerfacecolor='white',markersize=mk,markeredgewidth=2)
        else:
            plt.plot(pos[i][0],pos[i][1]+y0,color="black",marker='o',markerfacecolor='red',markersize=mk,markeredgewidth=2)

    for i in range(len(ai)):
        plt.annotate("",
                    xy=(ai[i],bi[i]),
                    xytext=(af[i],bf[i]),
                    size=30,va="center",ha="center",
                    arrowprops=dict(color='black',
                                    arrowstyle="-",
                                    connectionstyle="arc3,"+"rad=" + str(radd[i]) ,
                                    )
                    )
    for i in range(len(adi)):
        plt.annotate("",
                    xy=(adi[i],bdi[i]),
                    xytext=(adf[i],bdf[i]),
                    size=30,va="center",ha="center",
                    arrowprops=dict(color='pink',
                                    arrowstyle="-",
                                    connectionstyle="arc3,"+"rad=" + str(radd_dash[i]) ,
                                    )
                    )
    return
#def drawpicture2(ab,bb,spi,jj,ai,bi,ci,di,af,bf,cf,df,radd1,radd2,y0,pos):
#    #print('radd1',radd1)
#    #print('radd2',radd2)
#    for i in range(len(ab)):
#        for j in range(2):
#            bb[i][j]=bb[i][j]+y0-2   #跟排列vr跟vl順序而變
#    for i in range(len(ab)):
#        plt.plot(ab[i],bb[i],color="royalblue",linestyle="-",linewidth=str(3*jj[i]+1))
#    
#    for i in range(len(spi)):
#        if( spi[i] > 0):
#            plt.plot(pos[i][0],pos[i][1]+y0,color="black",marker='o',markerfacecolor='white',markersize=mk,markeredgewidth=2)
#        else:
#            plt.plot(pos[i][0],pos[i][1]+y0,color="black",marker='o',markerfacecolor='red',markersize=mk,markeredgewidth=2)
#
#    for i in range(len(ai)):
#        plt.annotate("",
#                    xy=(ai[i],bi[i]),
#                    xytext=(af[i],bf[i]),
#                    size=30,va="center",ha="center",
#                    arrowprops=dict(color='black',
#                                    arrowstyle="-",
#                                    connectionstyle="arc3,"+"rad=" + str(radd1[i]) ,
#                                    )
#                    )
#        plt.annotate("",
#                    xy=(cf[i],df[i]),
#                    xytext=(ci[i],di[i]),
#                    size=30,va="center",ha="center",
#                    arrowprops=dict(color='lightcoral',
#                                    arrowstyle="-",
#                                    connectionstyle="arc3,"+"rad=" + str(radd2[i]) ,
#                                    )
#                    )
#    return
radd=getrad(xir,xfr,spos)
radd_dash=getrad_dash(xid,xfd,spos)
y0=0
drawpicture1(xb,yb,rspin,xir,yir,xfr,yfr,xid,yid,xfd,yfd,radd,radd_dash,y0,spos)
#dashx=np.arange(0,np.int32(len(spos)/2),1)
#dashy=[]
#for i in range(len(dashx)):
#    dashy.append(y0+1.5)
#plt.plot(dashx,dashy,color='green',linestyle="-.")
#y0=y0+2
#for i in range(len(xir)):
#    yir[i]=yir[i]+2
#    yfr[i]=yfr[i]+2
#radd=getrad(xir,xfr,spos)
#drawpicture1(xb,yb,rspin,jj,xir,yir,xfr,yfr,radd,y0,spos)
#dashx=np.arange(0,np.int32(len(spos)/2),1)
#dashy=[]
#for i in range(len(dashx)):
#    dashy.append(y0+1.5)
#y0=y0+4
#for i in range(len(xir)):
#    yir[i]=yir[i]+y0-2
#    yfr[i]=yfr[i]+y0-2
#    yil[i]=yil[i]+y0
#    yfl[i]=yfl[i]+y0
#raddr=getrad(xir,xfr,spos)
#raddl=getrad(xil,xfl,spos)
#drawpicture2(xb,yb,lspin,jj,xil,yil,xir,yir,xfl,yfl,xfr,yfr,raddl,raddr,y0,spos)
#plt.plot(dashx,dashy,color='green',linestyle="-.")
plt.ylim(-1,y0+2)
plt.title('open chain SDRG')
plt.savefig('picture.pdf')
#plt.savefig('l40p0f00case2.png')
plt.show()
