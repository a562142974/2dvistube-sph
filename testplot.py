from matplotlib import pyplot as plt
import numpy as np
from sklearn import neighbors



def plotpositions(filename):
    #读取常量
    f=open(filename)


    xl,xr,yd,yu,h,N_particle,Neis,Nbx,Nby=[float(t) for t in f.readline().strip().split()]
    N_particle=int(N_particle)
    Neis=int(Neis)
    Nbx=int(Nbx)
    Nby=int(Nby)

    #读取所有点
    AP=[]
    for i in range(N_particle):
        AP.append(f.readline().strip().split())

    AP=np.array(AP,float)

    #读取邻居点
    neighbors=[]
    for i in range(Neis):
        neighbors.append(f.readline().strip().split())
    neighbors=np.array(neighbors,float)

    #读取中心点
    xc,yc=f.readline().strip().split()
    xc=float(xc)
    yc=float(yc)

    #读取网格线
    xb=[]
    for i in range(Nbx):
        xb=xb+f.readline().strip().split()
    xb=np.array(xb,float)

    yb=[]
    for i in range(Nby):
        yb=yb+f.readline().strip().split()
    yb=np.array(yb,float)

    #画图



    #画背景网格
    plt.xticks(xb)
    plt.yticks(yb)
    plt.grid()#画背景网格


    #画所有点
    plt.scatter(AP[:,0],AP[:,1],c='blue')
    #标签
    for i in range(N_particle):
        plt.text(AP[i,0],AP[i,1],str(int(AP[i,2]))+","+str(int(AP[i,3])))
    #画所有邻居点
    plt.scatter(neighbors[:,0],neighbors[:,1],c='red')
    #画圆
    theta=np.linspace(0,2*np.pi,500);
    plt.plot(xc+h*np.cos(theta),yc+h*np.sin(theta))
    plt.show()

plotpositions('positions.txt')