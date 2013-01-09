#!/usr/bin/env python
import Gassman as Gassman
import numpy as np
import math
from pylab import *
import numpy as np
import matplotlib.pyplot as plt


def logread140(fileName):
  depth,caliper,rhoC,depthP,bvi,bvw,cbw,phie,phit,gr,slowp,so,slows,vsh1,vsh2,boilD,dmrp, \
  dtc,dts,rhoD,facies,phit = np.loadtxt(fileName,unpack=True)
  return slowp,slows,dtc,dts,rhoC,so,phit,facies,depth,gr,bvi,cbw,vsh2

def logread159(fileName):
  depth,caliper,rhoC,depthP,bvi,bvw,cbw,phie,phit,gr,slowp,so,somril,slows,s1,s2,cal2, \
  dmroil,dmrp,dtc,dts,rhoD,t2k,facies = np.loadtxt(fileName,unpack=True)
  return slowp,slows,dtc,dts,rhoC,so,phit,facies,depth,gr,bvi,cbw,somril

def mergelogs(kSatD2,muDryD2,phit2,so2,gr2,facies2,vsh2,vpvs2,ip2, \
              depth2,kSatD,muDryD,phit,so,gr,facies,vsh,vpvs,ip,depth):
  depth = np.append(depth,depth2,axis=0) 
  ksat = np.append(kSatD,kSatD2,axis=0) 
  mudry = np.append(muDryD,muDryD2,axis=0)
  phit = np.append(phit,phit2,axis=0)
  so = np.append(so,so2,axis=0)
  gr = np.append(gr,gr2,axis=0)
  facies = np.append(facies,facies2,axis=0)
  vsh = np.append(vsh,vsh2,axis=0)
  vpvs = np.append(vpvs,vpvs2,axis=0)
  ip = np.append(ip,ip2,axis=0)

  return depth,ksat,mudry,phit,so,gr,facies,vsh,vpvs,ip

class Logs():

  def __init__(self,fileName,dmin,dmax,rhos):
    self._fileName = fileName
    self._dmin = dmin
    self._dmax = dmax
    self._rhos = rhos

  def parameters(self,slowp,slows,dtc,dts,rhoD,so,phit,facies,depth,gr,bvi,cbw,vsh2):
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
    self._lambdaC = np.zeros(n)
    self._lambdaD = np.zeros(n)
    self._kSatC = np.zeros(n)
    self._kSatD = np.zeros(n)
    self._vsh = np.zeros(n)
    self._rhodry = np.zeros(n)
    self._rho_cal = np.zeros(n)
    self._vp_cal = np.zeros(n)
    self._vpvs = np.zeros(n)
    self._ip = np.zeros(n)
    self._vsh2 = vsh2
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
    

    for i in range(n):

      if self._vsh[i]==self._gr[i]:
        self._vsh[i] = .8 #self._bvi[i] + self._cbw[i] 
      #print(self._vsh[i],self._bvi[i],self._cbw[i])

    #maxvsh = max(self._vsh)
    #print(maxvsh)
    #for i in range(n):
      #self._vsh[i] = self._vsh[i]
      #print(self._vsh[i],self._bvi[i],self._cbw[i])

    for i in range(n):
      if self._depth[i]>float(self._dmin) and self._depth[i]<float(self._dmax):
        self._vpC[i] = 1.0/self._slowp[i]*1000000.0
        self._rho_cal[i]=0.23*np.power(self._vpC[i],0.25)*1000
        self._vpC[i] *= 0.3048 # *.9
        self._vpD[i] = 1.0/self._dtc[i]*1000000.0
        self._vpD[i] *= 0.3048 *.93
        self._vsC[i] = 1.0/self._slows[i]*1000000.0
        self._vsC[i] *= 0.3048 #*.9
        self._vsD[i] = 1.0/dts[i]*1000000.0
        self._vsD[i] *= 0.3048 
        self._rhoD[i] = self._rhoD[i]*1000
        self._so[i] = self._so[i]/100
        self._phit[i] = self._phit[i]/100 *.93
        self._vpvs[i] = self._vpD[i]/self._vsD[i]
        #self._vp_cal[i] = (5.81-9.42*self._phit[i] - 2.21*.15)*1000
        self._ip[i] = self._vpD[i]*self._rhoD[i]

        #self._rho_cal[i]=0.23*np.power(self._vpC[i],0.25)*1000
        self._muSatC[i] = self._vsC[i]*self._vsC[i]*self._rhoD[i]
        self._muSatD[i] = self._vsD[i]*self._vsD[i]*self._rhoD[i]
        self._lambdaC[i] = self._vpC[i]*self._vpC[i]*self._rhoD[i]-2.0*self._muSatC[i]
        self._lambdaD[i] = self._vpD[i]*self._vpD[i]*self._rhoD[i]-2.0*self._muSatD[i]
        self._kSatC[i] = self._lambdaC[i]+2.0/3.0*self._muSatC[i]
        self._kSatD[i] = self._lambdaD[i]+2.0/3.0*self._muSatD[i]

        if self._vsh[i]!=self._gr[i]:
          self._vsh[i] = self._vsh2[i]
          if self._vsh2[i]>1:
            self._vsh[i] = 1
          if self._vsh2[i]<0:
            self._vsh[i] = 0
        self._vp_cal[i] = (5.81-9.42*self._phit[i] - 2.21*self._vsh[i])*1000


        self._rhodry[i] = self._rhos[0]*(1.0-self._vsh[i])+self._rhos[1]*self._vsh[i]
        self._muDryC[i] = self._vsC[i]*self._vsC[i]*(1.0-self._phit[i])*self._rhodry[i]
        self._muDryD[i] = self._vsD[i]*self._vsD[i]*(1.0-self._phit[i])*self._rhodry[i]
        #print(self._vsh[i],self._vp_cal[i])
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
    return (self._kSatC, self._kSatD, self._muDryC, self._muDryD, self._muSatC, \
           self._muSatD, self._phit, self._so, self._vsh, self._facies, self._depth,self._gr,self._vpvs,self._ip)

  def plot(self,title):
    self._title = title    
    #self._somril = somril   
    fig = figure(1, figsize=(15, 8))
    bx = plt.subplot(1,4,1)
    plt.suptitle(self._title)
    bx.set_xlabel('Density')
    bx.set_ylabel('Depth (ft)')
    p1=plot(self._rhoD,self._depth)
    p2=plot(self._rho_cal,self._depth,c='r')
    ylim(self._dmax,self._dmin) 
    bx.xaxis.set_major_locator(MaxNLocator(4))
    xlim(2000,2750) 
      
    plt.rcParams.update({'font.size': 18})
    ax = plt.subplot(1,4,2)
    ax.set_xlabel('P velocity (m/s)')
    #ax.set_ylabel('depth')
    p1=plot(self._vpC,self._depth)
    p2=plot(self._vp_cal,self._depth,c='r')
    ylim(self._dmax,self._dmin) 
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim(1500,5000) 
    #ax.set_title('Well 140-1')
    ax.set_yticklabels([])
    
    ax = plt.subplot(1,4,3)
    ax.set_xlabel('S velocity (m/s)')
    #ax.set_ylabel('depth')
    p1=plot(self._vsC,self._depth)
    #p2=plot(self._vsD,self._depth,c='r')
    ylim(self._dmax,self._dmin) 
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim(800,2000) 
    ax.set_yticklabels([])

    ax = plt.subplot(1,4,4)
    ax.set_xlabel('Gamma Ray')
    #ax.set_ylabel('depth')
    p1=plot(self._gr,self._depth)
    #p1=plot(vsh,depth)
    ylim(self._dmax,self._dmin) 
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim(0,150)     
    ax.set_yticklabels([])

    fig = figure(2, figsize=(11, 8))
    bx = plt.subplot(1,4,1)
    bx.set_xlabel('G')
    bx.set_ylabel('Depth (ft)')
    p1=plot(self._muDryC,self._depth,c='r')
    #p1=plot(self._muDryD,self._depth,c='k')
    p2=plot(self._muSatC,self._depth,c='b')
    #p3=plot(self._muSatD,self._depth,c='g')
    plt.suptitle(self._title)

    ylim(self._dmax,self._dmin) 
    bx.xaxis.set_major_locator(MaxNLocator(4))
    xlim(1e09,1e10) 

    ax = plt.subplot(1,4,2)
    ax.set_xlabel('Ksat')
    p2=plot(self._kSatC,self._depth,c='b')
    #p3=plot(self._kSatD,self._depth,c='r')

    ylim(self._dmax,self._dmin) 
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim(1e09,3e10) 
    ax.set_yticklabels([])

    cx = plt.subplot(1,4,3)
    cx.set_xlabel('Gamma Ray')
    p2=plot(self._gr,self._depth,c='b')
    cx.set_yticklabels([])

    ylim(self._dmax,self._dmin) 
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim(0,150) 


    cx = plt.subplot(1,4,4)
    cx.set_xlabel('SO')
    p2=plot(self._so,self._depth,c='b')
    #p2=plot(self._somril/100,self._depth,c='r')

    cx.set_yticklabels([])

    ylim(self._dmax,self._dmin) 
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim(0,0.5) 


  plt.show()
