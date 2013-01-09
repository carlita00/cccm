#!/usr/bin/env python
import Gassman as Gassman
import numpy as np
import math
from pylab import *
import numpy as np
import matplotlib.pyplot as plt

#---------Plotting----------------
class Plot:
  def __init__(self,kdry_l,gdry_l,kdry_c,gdry_c,kdryConstant,gdryConstant,phi,ks):
    self._kdry_l = kdry_l
    self._kdry_c = kdry_c
    self._gdry_l = gdry_l
    self._gdry_c = gdry_c
    self._kdryConstant = kdryConstant
    self._gdryConstant = gdryConstant
    self._ks = ks
    self._phi = phi
    

  def plotDry(self,kbri,koil,pref,ksat_well,gdry,phit,SO,vsh_iw,facies,vsh):

    #phid,phit,ksat_well,vsh_iw,g,SO,depth,well,facies,sid = np.loadtxt(fileName,unpack=True)
    kdry_well = Gassman.kdry(ksat_well,phit,SO,vsh,kbri,koil)

    n1 = len(ksat_well)
    dry_ratio = np.zeros(n1)
    #for i in range(n1):  
      #print(vsh[i])
    fig = figure(3, figsize=(10, 8))
    ax = plt.subplot(1,1,1)
    ax.set_xlabel('porosity $\phi$')
    ax.set_ylabel('Dry Bulk modulus $K$')
    #plot(self._phi[:],self._kdry_l[1,:],'b')
    #plot(self._phi[:],self._kdry_l[2,:],'b')
    plot(self._phi[:],self._kdry_l[pref,:],'b')
    plot(self._phi[:],self._kdry_c[:],'b')
    #plot(self._phi[:],self._kdryConstant[2,:],'b')
    #plot(self._phi[:],self._kdryConstant[9,:],'b')
    plot(self._phi[:],self._kdryConstant[18,:],'b')
    ylim(1.e08,4e10) 
    xlim(0,0.5) 
    p=ax.scatter(phit,kdry_well,c=vsh_iw)
    #p=ax.scatter(phit,ksat_well,c='r')
    p.set_clim([0,100])

    cbar= plt.colorbar(p)
    cbar.set_label('GR')
    plt.rcParams.update({'font.size': 20})

    fig = figure(4, figsize=(10, 8))
    bx = plt.subplot(1,1,1)
    bx.set_xlabel('porosity $\phi$')
    bx.set_ylabel('Dry Shear modulus $G$')
    plot(self._phi[:],self._gdry_l[pref,:],'b')
    plot(self._phi[:],self._gdry_c[:],'b')
    #plot(self._phi[:],self._gdryConstant[2,:],'b')
    #plot(self._phi[:],self._gdryConstant[9,:],'b')
    #plot(self._phi[:],self._gdryConstant[14,:],'b')
    ylim(1,8e10) 
    xlim(0,0.5) 

    p=bx.scatter(phit,gdry,c=vsh_iw)
    cbar= plt.colorbar(p)
    cbar.set_label('GR')
    p.set_clim([0,100])

    plt.rcParams.update({'font.size': 20})

  def plotDryratio(self,kbri,koil,ksat_well,gdry,phit,SO,vsh_iw,facies,depth,dmin,dmax,vsh):
    #phid,phit,ksat_well,vsh_iw,g,SO,depth,well,facies,sid = np.loadtxt(fileName,unpack=True)
    kdry_well= Gassman.kdry(ksat_well,phit,SO,vsh,kbri,koil)
    n1 = len(ksat_well)
    dry_ratio = np.zeros(n1)
    for i in range(n1):  
      #print(vsh[i],depth[i])  
      dry_ratio[i]=(kdry_well[i]+4.0/3.0*gdry[i])/gdry[i] 
      dry_ratio[i]=np.power(dry_ratio[i],0.5)
    fig = figure(5, figsize=(10, 8))
    bx = plt.subplot(1,1,1)
    bx.set_xlabel('Vp/Vs ratio')
    bx.set_ylabel('Depth (ft)')
    p=bx.scatter(dry_ratio,depth,c=SO)
    p.set_clim([0,0.5])
    cbar= plt.colorbar(p)
    cbar.set_label('Oil saturation')
    ylim(dmax,dmin) 
    xlim(1,2.4) 

    fig = figure(6, figsize=(10, 8))
    bx = plt.subplot(1,1,1)
    bx.set_xlabel('Vp/Vs ratio')
    bx.set_ylabel('Depth (ft)')
    plot(dry_ratio,depth,c='r')
    plot(SO,depth)
    xlim(1,2.4) 
    ylim(dmax,dmin)
    return kdry_well

  def plotfacies(self,ksat_well,gdry,phit,vsh_iw,facies,depth,kdry_well,pref):
    fig = figure(6, figsize=(8, 6))
    bx = plt.subplot(1,1,1)
    bx.set_xlabel('porosity $\phi$',fontsize=20)
    bx.set_ylabel('Dry Bulk modulus $K$ (Pa)',fontsize=20)
    axis([0.0, 0.4, 0, 4*1e10])
    colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.]
    n = len(depth)
    for i in range(n):
      for j in range(len(colors)):
        if vsh_iw[i] <30:
          if facies[i] == colors[j]:
            c=cm.hsv((j+1)/float(len(colors)),1)
            p=bx.scatter(phit[i],kdry_well[i],color=c)
            #print(facies[i])
    plot(self._phi[:],self._kdry_l[pref,:],'b')
    plot(self._phi[:],self._kdry_c[:],'b')
    plot(self._phi[:],self._kdryConstant[2,:],'b')
    plot(self._phi[:],self._kdryConstant[9,:],'b')
    plot(self._phi[:],self._kdryConstant[18,:],'b')

  def plottemplate(self,par,pref,vpvs,ip,so):
    fig = figure(7, figsize=(10, 8))
    bx = plt.subplot(1,1,1)
    bx.set_xlabel('Ip')
    bx.set_ylabel('vp/vs')
    a= par['vpvsBrineCo2']
    b= par['ipBrineCo2']
    c= par['vpvsBrineOil']
    d= par['ipBrineOil']


    p=bx.scatter(b[pref][0][:],a[pref][0][:],c='b')
    p=bx.scatter(b[pref][9][:],a[pref][9][:],c='r')
    p=bx.scatter(d[pref][3][:],c[pref][3][:],c='g')
    #p=bx.scatter(d[pref][9][:],c[pref][9][:],c=so)
    p=bx.scatter(ip,vpvs,c=so)


    cbar= plt.colorbar(p)
    cbar.set_label('Oil Saturation')
    #plot(self._phi[:],self._gdry_c[:],'b')
    ylim(1,3) 
    #xlim(0,0.5) 
    p.set_clim([0,0.3])

    plt.show()




plt.show()

