#!/usr/bin/env python
import Gassman as Gassman
import numpy as np
import math
from pylab import *
import numpy as np
import matplotlib.pyplot as plt


def logread140(fileName):
  depth,caliper,rhoC,depthP,bvi,bvw,cbw,phie,phit,gr,slowp,so,slows,vsh1,vsh2,boilD,dmrp, \
  dtc,dts,rhoD,facies,phit,rS,rL = np.loadtxt(fileName,unpack=True)
  return slowp,slows,dtc,dts,rhoD,so,phit,facies,depth,gr,bvi,cbw,vsh2,caliper,rS,rL

def logread159(fileName):
  depth,caliper,rhoC,depthP,bvi,bvw,cbw,phie,phit,gr,slowp,so,somril,slows,s1,s2,cal2, \
  dmroil,dmrp,dtc,dts,rhoD,t2k,facies,ml,rs,rl = np.loadtxt(fileName,unpack=True)
  return slowp,slows,dtc,dts,rhoD,so,phit,facies,depth,gr,bvi,cbw,somril,caliper,ml,rs,rl

def mergelogs(kSatD2,muDryD2,phit2,so2,gr2,facies2,vsh2,vpvs2,vs2,ip2, \
              depth2,cal2,rho2,kDry2,kSatD,muDryD,phit,so,gr1,facies1,vsh,vpvs,vs1,ip,depth1,cal,rho1,kDry1):
  depth = np.append(depth1,depth2,axis=0) 
  ksat = np.append(kSatD,kSatD2,axis=0) 
  kdry = np.append(kDry1,kDry2,axis=0) 
  mudry = np.append(muDryD,muDryD2,axis=0)
  phit = np.append(phit,phit2,axis=0)
  so = np.append(so,so2,axis=0)
  gr = np.append(gr1,gr2,axis=0)
  facies = np.append(facies1,facies2,axis=0)
  vsh = np.append(vsh,vsh2,axis=0)
  vpvs = np.append(vpvs,vpvs2,axis=0)
  ip = np.append(ip,ip2,axis=0)
  vs = np.append(vs1,vs2,axis=0)
  fig = figure(21, figsize=(10, 6))
  bx = plt.subplot(1,1,1)
  #bx.set_title('Paluxy and Tuscaloosa washover/transitioal facies')
  #bx.set_title('Tuscaloosa fluvial facies')

  bx.set_xlabel('Acoustic Impedance (kg/m$^2$s)*1e7')
  bx.set_ylabel(r'$V_{p}/V_{s}$')

  p=bx.scatter(ip/1e7,vpvs,c=so,edgecolor='none')
    
  cbar= plt.colorbar(p)
  cbar.set_label('Oil saturation')
  #plot(self._phi[:],self._gdry_c[:],'b')
  ylim(1.4,2.6) 
  xlim(0.4,1.1) 
  p.set_clim([0,0.5])

  fig.savefig('qc1.pdf')

  return depth,ksat,kdry,mudry,phit,so,gr,facies,vsh,vpvs,vs,ip

  
def plotqc(kSatD2,muDryD2,phit2,so2,gr2,facies2,vsh2,vpvs2,ip2,depth2,cal2,rho2,rhoCal2,vp2,vpCal2, \
           kSatD ,muDryD ,phit ,so ,gr1,facies1,vsh ,vpvs ,ip ,depth1,cal1,rho1,rhoCal1,vp1,vpCal1,rS,rL,ml,rs,rl):

  plt.rcParams.update({'font.size': 16,'legend.fontsize': 14})

  #-----Figure 1 de la tesis------#
  fig = figure(1, figsize=(9, 8))
  bx = plt.subplot(1,2,1)
  bx.set_xlabel('Gamma Ray')
  bx.set_ylabel('Depth (ft)')
  p=bx.plot(gr2,depth2,c='k')
  bx.xaxis.set_major_locator(MaxNLocator(4))
  bx.set_ylim(3440,3114) 
  bx.set_title('Well 159-2')

  colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0]
  n = len(depth2)
  for i in range(n):
    if depth2[i]>3186 and depth2[i]<3301:
      for j in range(len(colors)):
        if facies2[i] == colors[j]:
          if colors[j] >=5 and colors[j]<=9:
            j = 4
            c1='k'
            p=bx.scatter(gr2[i],depth2[i],c=c1)
          elif colors[j] >=50 and colors[j]<=90:
            j = 13
            c1='k'
            p=bx.scatter(gr2[i],depth2[i],c=c1)
          else:
            c1=cm.hsv((j+1)/float(len(colors)),1)
            p=bx.scatter(gr2[i],depth2[i],c=c1,edgecolor='none')
  bx.set_xlim(0,200)


  bx = plt.subplot(1,2,2)
  bx.set_xlabel('Gamma Ray')
  #bx2 = bx.twinx()

  #p=bx.plot(gr1,depth1,c='k')
  bx.xaxis.set_major_locator(MaxNLocator(4))
  p=bx.plot(gr1,depth1,c='k')
  bx.set_ylim(3440,3170) 
  bx.set_title('Well 140-1')
  #bx.set_yticklabels([])


  colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0]
  n = len(depth1)
  for i in range(n):
    if depth1[i]>3207 and depth1[i]<3440:

      for j in range(len(colors)):
        if facies1[i] == colors[j]:
          if colors[j] >=5 and colors[j]<=9:
            j = 4
            c1='k'
            p=bx.scatter(gr1[i],depth1[i],c=c1)
          elif colors[j] >=50 and colors[j]<=90:
            j = 13
            c1='k'
            p=bx.scatter(gr1[i],depth1[i],c=c1)
          else:
            c1=cm.hsv((j+1)/float(len(colors)),1)
            p=bx.scatter(gr1[i],depth1[i],c=c1,edgecolor='none')
  bx.set_xlim(0,200)

  p1 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((1)/float(len(colors)),1))
  p2 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((2)/float(len(colors)),1))
  p3 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((3)/float(len(colors)),1))
  p4 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((4)/float(len(colors)),1))
  p5 = Rectangle((0, 0), 0.2, 0.2, fc='k')

  p20 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((11)/float(len(colors)),1))
  p30 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((12)/float(len(colors)),1))
  p40 = Rectangle((0, 0), 0.2, 0.2, fc='k')

  labels = ('Tusc: Beach / Barrier bar' , 'Tusc: Washover', 'Tusc: Transitional','Tusc: Fluvial', \
          'Tusc: Poor rock', 'Paluxy: Distributary sand','Paluxy: mediocre distributary sand', \
          'Paluxy: Poor rock')
#  legend = plt.legend([p1,p2,p3,p4,p5,p20,p30,p40],labels,bbox_to_anchor=(0.5,-0.05), loc=2, borderaxespad=0., labelspacing=0.1)
  #box = bx.get_position()
  #bx.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])  
  legend = plt.legend([p1,p2,p3,p4,p5,p20,p30,p40],labels,bbox_to_anchor=(-0.1,-0.09), loc=9, fancybox=True, shadow=True, ncol=2)
#  fig.legend = plt.legend([p1,p2,p3,p4,p5,p20,p30,p40],labels,loc=9)
  
  ltext = gca().get_legend().get_texts()
  #fig.savefig('gr_facies.pdf')
  fig.savefig('gr_facies.pdf',bbox_extra_artists=(legend,), bbox_inches='tight')

  fig = figure(20, figsize=(8, 6))
  fig.suptitle("159-2")
  bx = plt.subplot(1,2,2)

  bx.set_xlabel('Resistivity (ohm-m)')
  #p=bx.plot(gr1,depth1,c='k')
  bx.xaxis.set_major_locator(MaxNLocator(4))
  bx.plot(rs,depth2,c='k',label='Res 10 inch')
  bx.plot(rl,depth2,c='b',label=" Res 90 inch ")
  bx.plot(ml,depth2,c='r',label=" Microlog")

  bx.set_xlim(0.,20)
  bx.set_yscale('log')  
  bx.set_ylim(3439,3114) 

  bx.set_yticklabels([])
  bx.legend()

 


  bx = plt.subplot(1,2,1)
  bx.set_xlabel('Gamma Ray')
  bx.set_ylabel('Depth (ft)')
  p=bx.plot(gr2,depth2,c='k')
  bx.xaxis.set_major_locator(MaxNLocator(4))
  bx.set_ylim(3439,3114) 


  colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0]
  n = len(depth2)
  for i in range(n):
    if depth2[i]>3186 and depth2[i]<3301:
      for j in range(len(colors)):
        if facies2[i] == colors[j]:
          if colors[j] >=5 and colors[j]<=9:
            j = 4
            c1='k'
            p=bx.scatter(gr2[i],depth2[i],c=c1)
          elif colors[j] >=50 and colors[j]<=90:
            j = 13
            c1='k'
            p=bx.scatter(gr2[i],depth2[i],c=c1)
          else:
            c1=cm.hsv((j+1)/float(len(colors)),1)
            p=bx.scatter(gr2[i],depth2[i],c=c1,edgecolor='none')
  bx.set_xlim(0,200)



  fig.savefig('invasion.pdf')

class Logs:
  def __init__(self,fileName,dmin,dmax,rhos):
    self._fileName = fileName
    self._dmin = dmin
    self._dmax = dmax
    self._rhos = rhos

  def parameters(self,slowp,slows,dtc,dts,rhoD,so,phit,facies,depth,gr,bvi,cbw,cal):
    #---Parameters calculation and everythig to MKS---------
    n=len(depth)
    self._vpC = np.zeros(n)
    self._vpD = np.zeros(n)
    self._vsC = np.zeros(n)
    self._vsD = np.zeros(n)
    self._vsh = np.zeros(n)
    self._muSatC = np.zeros(n)
    self._muSatD = np.zeros(n)
    self._muDryC = np.zeros(n)
    self._muDryD = np.zeros(n)
    self._kDryA = np.zeros(n)

    self._lambdaC = np.zeros(n)
    self._lambdaD = np.zeros(n)
    self._kSatC = np.zeros(n)
    self._kSatD = np.zeros(n)
    self._vsh = np.zeros(n)
    self._rhodry = np.zeros(n)
    self._rho_cal = np.zeros(n)
    self._vp_cal = np.zeros(n)
    self._vs_cal = np.zeros(n)
    self._vpvs = np.zeros(n)
    self._ip = np.zeros(n)
    #self._rhos = rhos
    self._slowp = slowp
    self._slows = slows
    self._dtc = dtc
    self._dts = dts
    self._rhoD = rhoD
    self._so = so
    self._phit = phit
    self._facies = facies
    self._depth = depth 
    self._gr = gr
    self._bvi = bvi
    self._cbw = cbw
    self._cal = cal

    maxGr = 150 #150 # 170 #200 #150#112 # max(gr)
    minGr = 35 #20 #35
    for i in range(n):
      #if self._
      self._vsh[i] = (self._gr[i]-minGr)/(maxGr-minGr) #self._bvi[i] + self._cbw[i] 
      if self._vsh[i]>=1:
        self._vsh[i]=1
      if self._vsh[i]<0:
        self._vsh[i]=0
    #maxvsh = max(self._vsh)
    #print(maxvsh)
    #for i in range(n):
      #self._vsh[i] = self._vsh[i]
      #print(self._vsh[i],self._bvi[i],self._cbw[i])

    for i in range(n):
      if self._depth[i]>float(self._dmin) and self._depth[i]<float(self._dmax):
        self._vpC[i] = 1.0/self._slowp[i]*1000000.0
        #self._rho_cal[i]=0.23*np.power(self._vpC[i],0.25)*1000

        self._vpC[i] *= 0.3048 # *.9
        self._vpD[i] = 1.0/self._dtc[i]*1000000.0
        #self._rho_cal[i]=0.23*np.power(self._vpD[i],0.25)*1000

        self._vpD[i] *= 0.3048 # *.93
        self._rho_cal[i]=-0.0115*self._vpD[i]/1000*self._vpD[i]/1000+0.261*self._vpD[i]/1000+1.515

        self._vsC[i] = 1.0/self._slows[i]*1000000.0
        self._vsC[i] *= 0.3048 #*.9
        self._vsD[i] = 1.0/dts[i]*1000000.0
        self._vsD[i] *= 0.3048 
        self._rhoD[i] = self._rhoD[i]*1000
        #print(rhoD[i],'puuu')
        self._so[i] = self._so[i]/100
        self._phit[i] = self._phit[i]/100  
        self._vs_cal[i] = (3.89-7.07*self._phit[i] - 2.04*self._vsh[i])*1000


        if self._vsD[i] <-200:
          self._vsD[i] = self._vs_cal[i]
       
        self._vpvs[i] = self._vpD[i]/self._vsD[i]
        #self._vp_cal[i] = (5.81-9.42*self._phit[i] - 2.21*.15)*1000
        self._ip[i] = self._vpD[i]*self._rhoD[i]

        #self._rho_cal[i]=0.23*np.power(self._vpC[i],0.25)*1000
        self._muSatC[i] = self._vsC[i]*self._vsC[i]*self._rhoD[i]
        self._muSatD[i] = self._vsD[i]*self._vsD[i]*self._rhoD[i]


        #self._muSatD[i] = self._vs_cal[i]*self._vs_cal[i]*self._rhoD[i]


        self._lambdaC[i] = self._vpC[i]*self._vpC[i]*self._rhoD[i]-2.0*self._muSatC[i]
        self._lambdaD[i] = self._vpD[i]*self._vpD[i]*self._rhoD[i]-2.0*self._muSatD[i]
        self._kSatC[i] = self._lambdaC[i]+2.0*self._muSatC[i]/3.0
        self._kSatD[i] = self._lambdaD[i]+2.0*self._muSatD[i]/3.0
        #print(self._vpD[i],self._rhoD[i])
        #if self._vsh[i]!=self._gr[i]:
          #self._vsh[i] = self._vsh2[i]
          #if self._vsh2[i]>1:
          #  self._vsh[i] = 1
          #if self._vsh2[i]<0:
            #self._vsh[i] = 0
        self._vp_cal[i] = (5.81-9.42*self._phit[i] - 2.21*self._vsh[i])*1000
        #self._muSatD[i] = self._vs_cal[i]*self._vsD[i]*self._rhoD[i] # OJOOOOOOOOOO
        #print(self._vsD[i],self._vs_cal[i],depth[i])
        vo  = 5590.0
        vfl = 1500.0 #1250.0
        #self._vp_cal[i] = (1.0-self._phit[i])*(1.0-self._phit[i])*vo+self._phit[i]*vfl
        #print(self._vp_cal[i],self._vsh[i])

        self._rhodry[i] = self._rhos[0]*(1.0-self._vsh[i])+self._rhos[1]*self._vsh[i]
        self._muDryC[i] = self._vsC[i]*self._vsC[i]*(1.0-self._phit[i])*self._rhodry[i]
        self._muDryD[i] = self._vsD[i]*self._vsD[i]*(1.0-self._phit[i])*self._rhodry[i]
        #print(self._vsh[i],self._vp_cal[i])
        kquartz = 37.0e09
        kclay = 25.0e9 #21.0e09
        gquartz = 44.0e09
        rhoquartz = 2650
        rhoclay = 2550
        gclay = 9.0e09
   
        mv = kquartz*(1.0-self._vsh[i])+kclay*self._vsh[i]
        mr = (1.0-self._vsh[i])/kquartz+self._vsh[i]/kclay
        ks = 0.5*(mv+1.0/mr)
    
        mv = gquartz*(1.0-self._vsh[i])+gclay*self._vsh[i]
        mr = (1.0-self._vsh[i])/gquartz+self._vsh[i]/gclay
        gs = 0.5*(mv+1.0/mr)

        self._kDryA[i] = self._muDryD[i]*ks/gs
        #self._kDryA[i] = self._muSatD[i]*ks/gs

      else:
        self._vpC[i] = -999.25
        self._vpD[i] = -999.25
        self._vsC[i] = -999.25
        self._vsD[i] = -999.25
        self._muSatC[i] = -999.25
        self._muSatD[i] = -999.25
        self._muDryC[i] = -999.25
        self._muDryD[i] = -999.25
        self._lambdaC[i] = -999.25
        self._lambdaD[i] = -999.25
        self._kSatC[i] = -999.25
        self._kSatD[i] = -999.25
        #vsh[i] = -999.25
        self._rhodry[i] = -999.25
        self._kDryA[i] = -999.25


    return (self._kSatC, self._kSatD, self._muDryC, self._muDryD, self._muSatC, \
            self._muSatD, self._phit, self._so, self._vsh, self._facies, self._depth, \
            self._gr,self._vpvs,self._vsD,self._ip,self._rho_cal,self._rhoD,self._vpD,self._vp_cal,self._kDryA)


  def plot(self,title):

    daxmin = 3186 #3125 #14 #3170
    daxmax = 3301 #3440
    dminfor = 3186 #3207
    dmaxfor = 3301 #3440
    #daxmin = 3207 #3170
    #daxmax = 3440
    #dminfor = 3207
    #dmaxfor = 3440

    self._title = title   
    plt.rcParams.update({'font.size': 18,'legend.fontsize': 15})

 
    fig = figure(2, figsize=(18, 8))

    cx = plt.subplot(1,5,1)
    cx.set_xlabel('Caliper (in)')
    cx.set_ylabel('Depth (ft)')
    cx.xaxis.set_major_locator(MaxNLocator(4))
    cx.set_ylim(daxmax,daxmin ) 

    #cx.set_title('Well 159-2')
    p2=plot(self._cal,self._depth,c='k')
    xlim(6,11)

    ax = plt.subplot(1,5,2)
    ax.set_xlabel('Gamma Ray')
    p=ax.plot(self._gr,self._depth,c='k')
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.set_ylim(daxmax,daxmin ) 
    #ax.set_title('Well 159-2')

    colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0]
    n = len(self._depth)
    for i in range(n):
      if self._depth[i]>dminfor and self._depth[i]<dmaxfor:
        for j in range(len(colors)):
          if self._facies[i] == colors[j]:
            if colors[j] >=5 and colors[j]<=9:
              j = 4
              c1='k'
              p=ax.scatter(self._gr[i],self._depth[i],c=c1)
            elif colors[j] >=50 and colors[j]<=90:
              j = 13
              c1='k'
              p=ax.scatter(self._gr[i],self._depth[i],c=c1)
            else:
              c1=cm.hsv((j+1)/float(len(colors)),1)
              p=ax.scatter(self._gr[i],self._depth[i],c=c1,edgecolor='none')
    ax.set_xlim(0,200)
    ax.set_yticklabels([])

    bx = plt.subplot(1,5,3)
    bx.set_xlabel(r'$\rho$ (g/cm$^3$)')
    bx.xaxis.set_major_locator(MaxNLocator(4))
    bx.set_ylim(daxmax,daxmin ) 
    #bx.set_title('Well 159-2')
    p2=plot(self._rhoD/1000,self._depth,c='k',label=r'$\rho$ Obs')
    p1=plot(self._rho_cal,self._depth,c='r',label=r'$\rho$ Cal')
    plt.legend()
    xlim(1.8,2.9)

    bx.set_title(title)


    bx.set_yticklabels([])

    ax = plt.subplot(1,5,4)
    ax.set_xlabel(r'$V_{p}$ (km/s)')
    #ax.set_ylabel('depth')
    p1=plot(self._vpD/1000,self._depth, c ='k',label=r"$V_{p}$ Obs")
    p2=plot(self._vp_cal/1000,self._depth,c='r',label=r"$V_{p}$ Cal")
    ylim(self._dmax,self._dmin) 
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim(1.5,4) 
    ax.set_ylim(daxmax,daxmin ) 
    #ax.set_title('Well 140-1')
    ax.set_yticklabels([])
    plt.legend()
    
    ax = plt.subplot(1,5,5)
    ax.set_xlabel(r'$V_{s}$ (km/s)')
    #ax.set_ylabel('depth')
    p1=plot(self._vsC/1000,self._depth,c='k',label="Vs Obs")
    p2=plot(self._vs_cal/1000,self._depth,c='r',label="Vs Cal")
    ylim(self._dmax,self._dmin) 
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim(.1,2.8) 
    ax.set_ylim(daxmax,daxmin ) 
    ax.set_yticklabels([])
    plt.legend()
    fig.savefig(self._title+'_QC.pdf', bbox_inches='tight')


    fig = figure(3, figsize=(8, 6))
    fig.suptitle(title)

    n = len(self._gr)
    lambd = np.zeros(n) 
    vp1    = np.zeros(n)
    vp    = np.zeros(n)

    kbri=2.81e9
    koil=0.829e9
    rhobri=1030
    rhooil=0.69e3

    ksatR,rhosatR,soR,kfluid=Gassman.kgassman(self._muSatC,self._phit,self._so,self._vsh,kbri,koil,rhobri,rhooil)


    for i in range(n):
      lambd[i]= ksatR[i]-2.0*self._muSatD[i]/3.0
      vp1[i]=np.power((lambd[i]+2*self._muSatC[i])/rhosatR[i],0.5)

    ax = plt.subplot(1,3,1)  
    ax.set_ylabel('Depth (feet)')

    ax.set_xlabel('Gamma Ray')
    p=ax.plot(self._gr,self._depth,c='k')
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.set_ylim(daxmax,daxmin ) 

    colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0]
    n = len(self._depth)
    for i in range(n):
      if self._depth[i]>dminfor and self._depth[i]<dmaxfor:
        for j in range(len(colors)):
          if self._facies[i] == colors[j]:
            if colors[j] >=5 and colors[j]<=9:
              j = 4
              c1='k'
              p=ax.scatter(self._gr[i],self._depth[i],c=c1)
            elif colors[j] >=50 and colors[j]<=90:
              j = 13
              c1='k'
              p=ax.scatter(self._gr[i],self._depth[i],c=c1)
            else:
              c1=cm.hsv((j+1)/float(len(colors)),1)
              p=ax.scatter(self._gr[i],self._depth[i],c=c1,edgecolor='none')
    ax.set_xlim(0,200)

    bx = plt.subplot(1,3,2)  
    bx.set_xlabel('Oil saturation')
    bx.plot(self._so*100,self._depth,c='k')
    bx.set_xlim(0,50)
    bx.xaxis.set_major_locator(MaxNLocator(4))
    bx.set_yticklabels([])
    bx.set_ylim(daxmax,daxmin ) 

    ax = plt.subplot(1,3,3)  

    ax.set_xlabel('Vp (km/s)')
    p1=ax.plot(self._vpD/1000,self._depth,c='k',label='Vp obs')
    p2=ax.plot(vp1/1000,self._depth,c='r',label='Vp cal')
    plt.legend()
    ax.set_yticklabels([])

    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.set_ylim(daxmax,daxmin ) 
    #ax.set_title('Well 159-2')
    ax.set_xlim(2,4.0)
   
    fig.savefig(self._title+'_vp.pdf')




  plt.show()
