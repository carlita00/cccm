#!/usr/bin/env python
import Gassman as Gassman
import numpy as np
import math
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.mlab as mll



tB=False
tF=False
tW=True
pL=False
minX=2275650
maxX=2285410
minY=636808
maxY=645065

if tB==True:
  vpvsfile='./Maps/VpVs_TB_norm'
  ipfile='./Maps/Zp_TB_norm'
  filename='TB_pressure.pdf'
  wellfile='./Logs/wellsTB'
if tW==True:
  vpvsfile='./Maps/VpVs_TW_norm2'
  ipfile='./Maps/Zp_TW_norm2'
  filename='TW_pressure.png'
  wellfile='./Logs/wellsTW'

if tF==True:
  vpvsfile='./Maps/VpVs_TW_norm'
  ipfile='./Maps/Zp_TW_norm'
  filename='TB_pressure.pdf'
  wellfile='./Logs/wellsTB'

if pL==True:
  vpvsfile='./Maps/VpVs_PL_norm'
  ipfile='./Maps/Zp_PL_norm'
  filename='PL_pressure.png'
  wellfile='./Logs/wellsPL'


xline,iline,x,y,vpvs_inversion = np.loadtxt(vpvsfile,unpack=True)
xline,iline,x,y,zp_inversion = np.loadtxt(ipfile,unpack=True)
pressure2interp = np.load('pressure_model_xzy.npy')
location = np.genfromtxt(wellfile,dtype=None)


def geti(x,ox,dx):
  i = int((x-ox)/dx)
  return i


plt.rcParams.update({'font.size': 15,'legend.fontsize': 12})


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
fig =  figure(1, figsize=(6, 5))
ax = plt.subplot(1,1,1)
p = plt.pcolormesh(xi,yi,zi,cmap='RdBu_r')
cbar=plt.colorbar(p)
ax.set_xlabel('Impedance % difference')
ax.set_ylabel('Vp/Vs ratio % difference')
#ax.set_title('Saturation ')
p = plt.axis([-20, 20, -20, 20])
cbar.set_label('Pore pressure change (MPa)')

fig.savefig('pressureInt.png', bbox_inches='tight')


n = zp_inversion.size

dp = np.zeros(n) # differencial pressure

todel=[]
for i in range(n):
  ix = geti(zp_inversion[i]*100,ox,dx) 
  iy = geti(vpvs_inversion[i]*100,oy,dy)
  print(ix,iy)
  ii=0
  if ix > nx-1 or iy>ny-1 or ix<0 or iy<0:
    dp[i] = 0.0 
    todel += [i]
  else:  
    if(zi.mask[iy,ix]):
      dp[i] = 0.0
      todel += [i]
    else:
      dp[i] = zi[iy,ix]

ndel = len(todel)
dp = np.delete(dp,todel)
x  = np.delete(x,todel)
y  = np.delete(y,todel)


  #print(dp[i],ix,iy,zp_inversion[i]*100,vpvs_inversion[i]*100,dx,dy)

# Size of regular grid
ny, nx = 500, 500

# Generate a regular grid to interpolate the data.
xi = np.linspace(min(x), max(x), nx)
yi = np.linspace(min(y), max(y), ny)
rangex =max(xi)-min(xi)
rangey =max(yi)-min(yi)

ratio=rangex/rangey


xi, yi = np.meshgrid(xi, yi)

# Interpolate using delaunay triangularization 
zzi = mll.griddata(x,y,dp,xi,yi)

fig =  figure(2)
ax = plt.subplot(1,1,1)

aspectratio=1.0
ratio_default=(ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0])
ax.set_aspect(ratio_default*aspectratio)
p = plt.pcolormesh(xi,yi,zzi,cmap='RdBu_r')
#cmap = mpl.colors.ListedColormap(['b','w','r'])
cbar=plt.colorbar(p)
cbar.set_label('Pore pressure change (MPa)')
cbar.set_clim(-5,5)
#p = plt.axis([min(x), max(x), min(y), max(y)])
p = plt.axis([minX, maxX, minY, maxY])



n = location.size

for i in range(n):
  if(location[i][3]=='injectorc'):
    ax.scatter(location[i][1],location[i][2] ,c= 'r',marker='^',s=200, lw = 2)
    ax.annotate(location[i][0], xy=(location[i][1], location[i][2]), xycoords='data',
    xytext=(-20, -20), textcoords='offset points', fontsize=10)
  if(location[i][3]=='injectorw'):
    ax.scatter(location[i][1],location[i][2] ,c= 'b',marker='^',s=200, lw = 2)
    ax.annotate(location[i][0], xy=(location[i][1], location[i][2]), xycoords='data',
    xytext=(-20, -20), textcoords='offset points', fontsize=10)


  if(location[i][3]=='producer'):
    ax.scatter(location[i][1],location[i][2] ,c= 'g',marker='o',s=100, lw = 2)
    ax.annotate(location[i][0], xy=(location[i][1], location[i][2]), xycoords='data',
    xytext=(-20, -20), textcoords='offset points', fontsize=10)
  if(location[i][3]=='producerw'):
    ax.scatter(location[i][1],location[i][2] ,c= 'b',marker='o',s=100, lw = 2)
    ax.annotate(location[i][0], xy=(location[i][1], location[i][2]), xycoords='data',
    xytext=(-20, -20), textcoords='offset points', fontsize=10)
ax.set_ylabel('Y ')
ax.set_xlabel('X ')
ax.set_yticklabels([])
ax.set_xticklabels([])

fig.savefig(filename)


N_inv = zp_inversion.size
#N_well = xwell.size

discrimination = np.zeros(N_inv)


for i in range(N_inv):
  if zp_inversion[i] <=0.0 and vpvs_inversion[i] >= 0.0:
    ax.scatter(zp_inversion[i]*100,vpvs_inversion[i]*100 ,c= 'g',edgecolor='none')
    discrimination[i]= 0.5
  elif zp_inversion[i] <=0.0 and vpvs_inversion[i] <= 0.0:
    ax.scatter(zp_inversion[i]*100,vpvs_inversion[i]*100 ,c= 'r',edgecolor='none')
    discrimination[i]= 1.5

  elif zp_inversion[i] >=0.0 and vpvs_inversion[i] >= 0.0:
    ax.scatter(zp_inversion[i]*100,vpvs_inversion[i]*100 ,c= 'b',edgecolor='none')
    discrimination[i]= 2.5

  elif zp_inversion[i] >=0.0 and vpvs_inversion[i] <= 0.0:
    ax.scatter(zp_inversion[i]*100,vpvs_inversion[i]*100 ,c= 'y',edgecolor='none')
    discrimination[i]= 3.5

  else:
    ax.scatter(zp_inversion[i]*100,vpvs_inversion[i]*100 ,c= 'w',edgecolor='none')
    discrimination[i]= 4.5

# Size of regular grid
ny, nx = 500, 500

xline,iline,x,y,vpvs_inversion = np.loadtxt(vpvsfile,unpack=True)
xline,iline,x,y,zp_inversion = np.loadtxt(ipfile,unpack=True)
pressure2interp = np.load('pressure_model_xzy.npy')

# Generate a regular grid to interpolate the data.
xi = np.linspace(min(x), max(x), nx)
yi = np.linspace(min(y), max(y), ny)
xi, yi = np.meshgrid(xi, yi)

# Interpolate using delaunay triangularization 
zi = mll.griddata(x,y,discrimination,xi,yi)

#fig =  figure(3, figsize=(6*ratio*1.3, 6))
fig =  figure(3)
ax = plt.subplot(1,1,1)
aspectratio=1.0
ratio_default=(ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0])
ax.set_aspect(ratio_default*aspectratio)
cmap = mpl.colors.ListedColormap(['Gold','Crimson','RoyalBlue','Indigo'])
bounds=[0,1,2,3,4]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
ax.set_xlabel('Y (m)')
ax.set_ylabel('X (m)')
p = plt.pcolormesh(xi,yi,zi,cmap=cmap)
#ax = plt.scatter(x,y,c=vpvs_inversion)
cbar=plt.colorbar(p, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0, 1, 2,3,4])
#p = plt.axis([min(x), max(x), min(y), max(y)])
p = plt.axis([minX, maxX, minY, maxY])
bounds = np.arange(5)
vals = bounds[:-1]
cbar.set_ticks(vals + .5)
cbar.set_ticklabels(['P+', 'CO2+', 'Brine+', 'P-'])
cbar.set_clim(0,4)




n = location.size

for i in range(n):
  if(location[i][3]=='injectorc'):
    ax.scatter(location[i][1],location[i][2] ,c= 'r',marker='^',s=200, lw = 2)
    ax.annotate(location[i][0], xy=(location[i][1], location[i][2]), xycoords='data',
    xytext=(-20, -20), textcoords='offset points', fontsize=10)
  if(location[i][3]=='injectorw'):
    ax.scatter(location[i][1],location[i][2] ,c= 'b',marker='^',s=200, lw = 2)
    ax.annotate(location[i][0], xy=(location[i][1], location[i][2]), xycoords='data',
    xytext=(-20, -20), textcoords='offset points', fontsize=10)


  if(location[i][3]=='producer'):
    ax.scatter(location[i][1],location[i][2] ,c= 'g',marker='o',s=100, lw = 2)
    ax.annotate(location[i][0], xy=(location[i][1], location[i][2]), xycoords='data',
    xytext=(-20, -20), textcoords='offset points', fontsize=10)
  if(location[i][3]=='producerw'):
    ax.scatter(location[i][1],location[i][2] ,c= 'b',marker='o',s=100, lw = 2)
    ax.annotate(location[i][0], xy=(location[i][1], location[i][2]), xycoords='data',
    xytext=(-20, -20), textcoords='offset points', fontsize=10)
ax.set_ylabel('Y ')
ax.set_xlabel('X ')
ax.set_yticklabels([])
ax.set_xticklabels([])

fig.savefig('change'+filename)


plt.show()

