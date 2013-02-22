#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt


def kdry(ksatWell,phiWell,soWell,vsh,kbri,koil):
  n1 = len(ksatWell)
  kdry   = np.zeros(n1)
  kquartz = 37.0e9
  kclay = 25.0e09 #15
  gquartz = 44.0e9
  rhoquartz = 2650
  rhoclay = 2550
  gclay = 9.0e9

  #K dry calculation by using Gassman eq in the wells
  for i in range (n1):
    mv = kquartz*(1.0-vsh[i])+kclay*vsh[i]
    mr = (1.0-vsh[i])/kquartz+vsh[i]/kclay
    ks = 0.5*(mv+1.0/mr)
    
    mv = gquartz*(1.0-vsh[i])+gclay*vsh[i]
    mr = (1.0-vsh[i])/gquartz+vsh[i]/gclay
    gs = 0.5*(mv+1.0/mr)
    #soWell[i] =soWell[i]*0.0

    kfluid = (1.0-soWell[i])/kbri+soWell[i]/koil
    kfluid = 1.0/kfluid
    #print(phiWell[i],vsh[i],ksatWell[i],kfluid,soWell[i]) 
    tmp1 = ksatWell[i]*(phiWell[i]*ks/kfluid+1.0-phiWell[i])-ks
    tmp2 = phiWell[i]*ks/kfluid+ksatWell[i]/ks-1.0-phiWell[i] 
    kdry[i] = tmp1/tmp2
  return kdry 

def kgassman(kdry,phi,so,vsh,kbri,koil,rhobri,rhooil):
  n = len(kdry)
  #print(max(kdry),max(so),max(vsh),kbri,koil)
  ksat   = np.zeros(n)
  rhosat   = np.zeros(n)
  soR   = np.zeros(n)
  kfluid   = np.zeros(n)

  kquartz = 37.0e9
  kclay = 25.0e09
  gquartz = 44.0e9
  rhoquartz = 2650
  rhoclay = 2550
  gclay = 9.0e9

  for i in range(n):
    #so[i] =so[i]*0.0
    soR[i] =so[i] #*0.00001

    mv = kquartz*(1.0-vsh[i])+kclay*vsh[i]
    mr = (1.0-vsh[i])/kquartz+vsh[i]/kclay
    ks = 0.5*(mv+1.0/mr)
    #print(so[i]) 
    mv = gquartz*(1.0-vsh[i])+gclay*vsh[i]
    mr = (1.0-vsh[i])/gquartz+vsh[i]/gclay
    gs = 0.5*(mv+1.0/mr)
    rhos = rhoquartz*(1.0-vsh[i])+rhoclay*vsh[i]

    kfluid[i] = (1.0-so[i])/kbri+so[i]/koil
    kfluid[i] = 1.0/kfluid[i]
    #print(kfluid)
    rhofluid = (1.0-so[i])*rhobri+so[i]*rhooil
    tmp3=(phi[i]/kfluid[i]+(1.0-phi[i])/ks-kdry[i]/(ks*ks))
    tmp1=((1.0-kdry[i]/ks)*(1.0-kdry[i]/ks))
    ksat[i]=  kdry[i]+tmp1/tmp3 
    rhosat[i]= rhooil*phi[i]+(1.0-phi[i])*rhos 

  return ksat,rhosat,soR,kfluid
  
class Gassman:
  def __init__(self,ks,rhos,phi):
    self._ks   = ks
    self._rhos = rhos
    self._phi  = phi

  def fluidSub2_p(self,kbri,kco2,koil,rhobri,rhoco2,rhooil,kdry,co2,p,oilInit,brineInit):
    '''K  sat calculation by using Gassman.
    This modules works for combining Bounds with fluid effect
    '''
    self._kbri = kbri
    self._rhobri = rhobri
    self._koil = koil
    self._rhooil = rhooil
    self._rhoco2 = rhoco2
    self._kco2 = kco2
    self._kdry = kdry
    self._co2  = co2
    self._p    = p
    self._rhoco2 =  rhoco2
    n1 = len(self._phi)
    n3 = len(self._p)
    n2 = len(self._co2)
    ksat   = np.zeros(((n3,n2,n1))) 
    rhosat = np.zeros(((n3,n2,n1))) 
    for i3 in range(n3):
      for i2 in range(n2):
        tmp3=0.0
        for i1 in range (n1):
          kfluid = (brineInit-self._co2[i2]/2.0)/self._kbri[i3]+(oilInit-self._co2[i2]/2.0)/self._koil[i3]+self._co2[i2]/self._kco2[i3]
          kfluid = 1.0/kfluid
          rhofluid = (brineInit-self._co2[i2]/2)*self._rhobri[i3]+(oilInit-self._co2[i2]/2)*self._rhooil[i3]+self._co2[i2]*self._rhoco2[i3]
          tmp3=(self._phi[i1]/kfluid+(1.0-self._phi[i1])/self._ks-self._kdry[i3][i1]/(self._ks*self._ks))
          tmp1=((1.0-self._kdry[i3][i1]/self._ks)*(1.0-self._kdry[i3][i1]/self._ks))
          ksat[i3][i2][i1]=  self._kdry[i3][i1]+tmp1/tmp3 
          rhosat[i3][i2][i1]= rhofluid*self._phi[i1]+(1.0-self._phi[i1])*self._rhos 
    return (ksat,rhosat) 

  def fluidSub2_st(self,kbri,kco2,koil,rhobri,rhoco2,rhooil,kdry,co2,p,oilInit,brineInit):
    '''K  sat calculation by using Gassman.
    This modules works for combining Bounds with fluid effect
    '''
    self._kbri = kbri
    self._rhobri = rhobri
    self._koil = koil
    self._rhooil = rhooil
    self._rhoco2 = rhoco2
    self._kco2 = kco2
    self._kdry = kdry
    self._co2  = co2
    self._p    = p
    self._rhoco2 =  rhoco2
    n1 = len(self._phi)
    n3 = len(self._p)
    n2 = len(self._co2)
    ksat   = np.zeros(((n3,n2,n1))) 
    rhosat = np.zeros(((n3,n2,n1))) 
    for i3 in range(n3):
      for i2 in range(n2):
        tmp3=0.0
        for i1 in range (n1):
          kfluid =(1.0-self._co2[i2])/self._kbri[i3]+self._co2[i2]/self._kco2[i3]
          kfluid = 1.0/kfluid
          rhofluid = (1.0-self._co2[i2])*self._rhobri + self._co2[i2]*self._rhoco2
          tmp3=(self._phi[i1]/kfluid+(1.0-self._phi[i1])/self._ks-self._kdry[i3][i1]/(self._ks*self._ks))
          tmp1=((1.0-self._kdry[i3][i1]/self._ks)*(1.0-self._kdry[i3][i1]/self._ks))
          ksat[i3][i2][i1]=  self._kdry[i3][i1]+tmp1/tmp3 
          rhosat[i3][i2][i1]= self._rhoco2[i3]*self._phi[i1]+(1.0-self._phi[i1])*self._rhos 
    return (ksat,rhosat) 


  def fluidSub2(self,kbri,kco2,rhobri,rhoco2,kdry,co2):
    '''K  sat calculation by using Gassman.
    This modules works for combining Bounds with fluid effect
    '''
    self._kbri = kbri
    self._rhobri = rhobri
    self._rhoco2 = rhoco2
    self._kco2 = kco2
    self._kdry = kdry
    self._co2  = co2
    self._rhoco2 =  rhoco2
    n1 = len(self._phi)
    n3 = len(self._p)
    n2 = len(self._co2)
    ksat   = np.zeros((n2,n1)) 
    rhosat = np.zeros((n2,n1)) 
    for i2 in range(n2):
      tmp3=0.0
      for i1 in range (n1):
        kfluid = (1.0-self._co2[i2])/self._kbri+self._co2[i2]/self._kco2
        kfluid = 1.0/kfluid
        rhofluid = (1.0-self._co2[i2])*self._rhobri+self._co2[i2]*self._rhoco2
        tmp3=(self._phi[i1]/kfluid+(1.0-self._phi[i1])/self._ks-self._kdry[i1]/(self._ks*self._ks))
        tmp1=((1.0-self._kdry[i1]/self._ks)*(1.0-self._kdry[i1]/self._ks))
        ksat[i2][i1]=  self._kdry[i1]+tmp1/tmp3 
        rhosat[i2][i1]= self._rhoco2*self._phi[i1]+(1.0-self._phi[i1])*self._rhos 
    return (ksat,rhosat) 


        



