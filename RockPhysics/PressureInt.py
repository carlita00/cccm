#!/usr/bin/env python
import Gassman as Gassman
import numpy as np
import math
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.mlab as mll



class sampling:
  def __init__(self,nx,dx,ox):
    self.nx = nx
    self.dx = dx
    self.ox = ox

  def geti(self,x):
    return int((x-self.ox)/self.dx)

  def checki(self,i):
    if i< 0:
      a = False
    elif i > self.nx:
      a = False
    else:
      a = True
    return a



def geti(x,ox,dx):
  i = int((x-ox)/dx)
  return i


def gridding(g3,g4,g5,g6,vpvsfile,ipfile,n3,n2,pref,p):

  pressure2interp = np.zeros((n3*n2*2,3))
  ii = 0
  pconfinment=22.5e6
  for i3 in range (n3):
    for i2 in range(n2):
      pressure2interp[ii][0] =  g4[i3][i2][320]
      pressure2interp[ii][1] =  g3[i3][i2][320]
      pressure2interp[ii][2] =  ((pconfinment-p[i3])-(pconfinment-p[pref]))/1.0e6 
      ii+=1

      pressure2interp[ii][0] =  g6[i3][i2][320]
      pressure2interp[ii][1] =  g5[i3][i2][320]
      pressure2interp[ii][2] =  ((pconfinment-p[i3])-(pconfinment-p[pref]))/1.0e6 
      ii+=1

  np.save('pressure_model_xzy.npy',pressure2interp)


def maps(vpvsfile,ipfile,n3,n2,pref,p):

  xline,iline,x,y,vpvs_inversion = np.loadtxt(vpvsfile,unpack=True)
  xline,iline,x,y,zp_inversion = np.loadtxt(ipfile,unpack=True)
  pressure2interp = np.load('pressure_model_xzy.npy')

  ####Gridding pressure###########################
  # Size of regular grid
  ny, nx = 100, 101

  # Generate a regular grid to interpolate the data.
  ox = -40.
  oy = -25.
  dx = (20-ox)/nx
  dy = (20-oy)/ny

  xi = np.linspace(ox, 20., nx)
  yi = np.linspace(oy, 20., ny)
  xi, yi = np.meshgrid(xi, yi)

  # Interpolate using delaunay triangularization 
  zi = mll.griddata(pressure2interp[:,0],pressure2interp[:,1],pressure2interp[:,2],xi,yi)
  fig =  figure(33, figsize=(6, 5))
  ax = plt.subplot(1,1,1)
  p = plt.pcolormesh(xi,yi,zi)
  cbar=plt.colorbar(p)

  n = zp_inversion.size

  dp = np.zeros(n) # differencial pressure

  for i in range(n):
    ix = geti(zp_inversion[i]*100,ox,dx) 
    iy = geti(vpvs_inversion[i]*100,oy,dy)
    if ix > nx or iy>ny or ix<0 or iy<0:
      dp[i] = 0.0 
    else:  
      if(zi.mask[iy,ix]):
        dp[i] = 0.0
      else:
        dp[i] = zi[iy,ix]
    #print(dp[i],ix,iy,zp_inversion[i]*100,vpvs_inversion[i]*100,dx,dy)

  # Size of regular grid
  ny, nx = 500, 500

  # Generate a regular grid to interpolate the data.
  xi = np.linspace(min(x), max(x), nx)
  yi = np.linspace(min(y), max(y), ny)
  xi, yi = np.meshgrid(xi, yi)

  # Interpolate using delaunay triangularization 
  zzi = mll.griddata(x,y,dp,xi,yi)

  fig =  figure(40, figsize=(10, 8))
  ax = plt.subplot(1,1,1)
  p = plt.pcolormesh(xi,yi,zzi,cmap='RdBu')
  cbar.set_clim(-5,5)
  cmap = mpl.colors.ListedColormap(['b','w','r'])
  cbar=plt.colorbar(p)
  '''
  name,xwell,ywell,types = np.loadtxt('./Logs/wells',dtype=None,unpack=True)

  n = xwell.size

  for i in range(n):
    if(types[i]=='injector'):
      ax.scatter(xwell[i],ywell[i] ,c= 'r',marker='1',s=100)
    if(types[i]=='producer'):
      ax.scatter(xwell[i],ywell[i] ,c= 'g',marker='1',s=50)
  '''

  fig.savefig('TW_pressure.pdf')
plt.show()

