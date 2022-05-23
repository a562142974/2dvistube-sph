from fileinput import filename
from matplotlib import pyplot as plt
import numpy as np
from sklearn import neighbors
from matplotlib import cm

filename='rho0.3.txt'
f=open(filename,'r')

nx,ny=f.readline().strip().split()
nx=int(nx)
ny=int(ny)

x=np.array(f.readline().strip().split(),np.double)
y=np.array(f.readline().strip().split(),np.double)

X, Y = np.meshgrid(x, y)

Z=np.array(f.readline().strip().split(),np.double)

for t in range(ny-1):
    temp=np.array(f.readline().strip().split(),np.double)
    Z=np.vstack((Z,temp))

norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
fig, ax = plt.subplots()
cset1 = ax.contourf(X, Y, Z)
ax.set_aspect(1)
plt.show()