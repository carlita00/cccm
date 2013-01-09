#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt

def hillAverage (k,f):
    n = len(k)    
    mv = 0
    mr =0 
    for i in range(n):
      mv += k[i]*f[i]
      mr += f[i]/k[i]
    mvr = 0.5*(mv+1.0/mr)
    return mvr

class ElasticMedia:
  def __init__(self,ks,gs,phiC,phi,p):
    self._ks   = ks
    self._gs   = gs
    self._phiC = phiC
    self._phi  = phi
    self._p    = p
 
  # Modified Hashin-Shtrikman bound
  
  def kHS_lower(self,khm,ghm):
    self._khm = khm
    self._ghm = ghm
    n2 = len(self._p)
    n1 = len(self._phi)    
    kdry = np.zeros((n2,n1))
    for i2 in range(n2):
      for i1 in range(n1):
        kdry[i2][i1]= (self._phi[i1]/self._phiC)/(self._khm[i2]+4.0/3.0*self._ghm[i2])+ \
                       (1.0-self._phi[i1]/self._phiC)/(self._ks +4.0/3.0*self._ghm[i2])
        kdry[i2][i1]= 1.0/kdry[i2][i1] - 4.0/3.0*self._ghm[i2]
    return kdry

  def gHS_lower(self,ghm,zhm):
    self._zhm = zhm
    self._ghm = ghm
    n2 = len(self._p)
    n1 = len(self._phi)       
    gdry = np.zeros((n2,n1))
    for i2 in range(n2):
      for i1 in range(n1):
        gdry[i2][i1] = (self._phi[i1]/self._phiC)/(self._ghm[i2]+self._zhm[i2])+ \
                       (1.0-self._phi[i1]/self._phiC)/(self._gs+self._zhm[i2])
        gdry[i2][i1] = 1.0/gdry[i2][i1]-self._zhm[i2]
    return gdry 
  
  def kHS_upper(self,khm):
    self._khm = khm
    n2 = len(self._p)
    n1 = len(self._phi)    
    kdry = np.zeros((n2,n1))
    for i2 in range(n2):
      for i1 in range(n1):
        kdry[i2][i1]= (self._phi[i1]/self._phiC)/(self._khm[i2]+4.0/3.0*self._gs)+ \
                       (1.0-self._phi[i1]/self._phiC)/(self._ks +4.0/3.0*self._gs)
        kdry[i2][i1]= 1.0/kdry[i2][i1] - 4.0/3.0*self._gs
    return kdry

  def gHS_upper(self,ghm,z):
    self._z   = z
    self._ghm = ghm
    n2 = len(self._p)
    n1 = len(self._phi)       
    gdry = np.zeros((n2,n1))
    for i2 in range(n2):
      for i1 in range(n1):
        gdry[i2][i1] = (self._phi[i1]/self._phiC)/(self._ghm[i2]+self._z)+ \
                       (1.0-self._phi[i1]/self._phiC)/(self._gs+self._z)
        gdry[i2][i1] = 1.0/gdry[i2][i1]-self._z
    return gdry 
  

   



