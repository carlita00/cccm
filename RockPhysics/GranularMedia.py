#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt



class GranularMedia:

  def __init__(self,n,phiC,gs,ks,prs,p):
    self._n    = n
    self._phiC = phiC
    self._gs   = gs 
    self._ks   = ks
    self._p    = p
    self._prs  = prs

  def applyKhm (self):
    n1  = len(self._p)
    k   = np.zeros(n1)
    po  = 1.0/3.0
    for i in range(n1):
        k[i] = (self._n*self._n*np.power(1.0-self._phiC,2)*self._gs*self._gs)
        k[i] /= (18.0*np.power(np.pi,2)*np.power(1.0-self._prs,2))
        k[i] *= self._p[i] 
        k[i] = np.power(k[i],po)
    return k

  def applyGhm (self):
    n1  = len(self._p)
    g   = np.zeros(n1)
    for i in range(n1):
      g[i] = (3.0*self._n*self._n*np.power(1.0-self._phiC,2.0)*self._gs*self._gs)
      g[i] /= (2.0*np.power(np.pi,2.0)*np.power(1.0-self._prs,2.0))
      g[i] *= self._p[i] 
      g[i] =((5.0-4.0*self._prs)/(10.0-5.0*self._prs)) *np.power(g[i],0.333333) 
    return g


  def applyZhm(self,khm,ghm):
    n1 = len(self._p)
    z  = np.zeros(n1)
    for i in range(n1):
        z[i] = ghm[i]/6.0*(9.0*khm[i]+8.0*ghm[i])
        z[i]/= khm[i]+2.0*ghm[i]  
    return z

  def applyZ(self):
    z = self._gs/6.0*(9.0*self._ks+8.0*self._gs)
    z/= self._ks+2.0*self._gs  
    return z

  def kdryCement(self,c,phi,kc,gc,prc,mc):  
    self._phi = phi   
    n1        = len(self._phi)
    self._kc  = kc  
    self._c   = c
    self._gc  = gc
    self._prc = prc 
    self._mc  = mc
    alpha     = np.zeros(n1)
    k  = np.zeros(n1)
    for i in range(n1):
      n2 = 20.0 - 34.0 *self._phiC +14.0*np.power(self._phiC,2) 
      n2*= self._c

      alpha[i]= (2.0/3.0)*(self._phiC-self._phi[i])/(1.0-self._phiC)
      alpha[i]= np.power(alpha[i], 0.5)
      #alpha[i]=(self._phiC-self._phi[i]) /(3.0*n2*(1.0-self._phiC))
      #alpha[i]= 2*np.power(alpha[i], 0.25)

    vn = 2.0*self._gc*(1.0-self._prs)*(1.0-self._prc)/(math.pi*self._gs*(1.0-2.0*self._prc))

    #vn/= math.pi*self._gs*(1.0-2.0*self._prc)
    print(vn,"vn")        
    an = -.0245153 *np.power(vn,-1.3646)
    bn = 0.20405 *np.power(vn,-0.89008)
    cn = 0.00024649 *np.power(vn,-1.9864)
    for i in range(n1):
      #n2 = 20.0 - 34.0 *phi[i] +14.0*np.power(phi[i],2) 
      n2 = 20.0 - 34.0 *self._phiC +14.0*np.power(self._phiC,2) 
      n2*= self._c

      sn   = an*np.power(alpha[i],2.0)+bn*alpha[i]+cn
      k[i] = n2*(1.0-self._phiC)*self._mc*sn/6.0
    return k

  def gdryCement(self,c,phi,gc,kdryC):
    print(gc)
    self._phi = phi 
    n1    = len(self._phi)
    self._kdryC  = kdryC 
    self._gc  = gc
    alpha = np.zeros(n1)
    g = np.zeros(n1)
    for i in range(n1):
      alpha[i] = (2.0/3.0)*(self._phiC-self._phi[i])/(1.0-self._phiC)
      alpha[i] = np.power(alpha[i], 0.5)  
      n2 = 20.0 - 34.0 *self._phiC +14.0*np.power(self._phiC,2) 
      n2*= self._c
      #alpha[i]=(self._phiC-self._phi[i]) /(3.0*n2*(1.0-self._phiC))
      #alpha[i]= 2*np.power(alpha[i], 0.25)


      #alpha[i]= (2.0/3.0)*(self._phiC-self._phi[i])/(1.0-self._phiC)
      #alpha[i]= np.power(alpha[i], 0.5)
  
    vtao = self._gc/(math.pi*self._gs)
    tmp  = 0.079*np.power(self._prs,2)+0.1754*self._prs-1.342 
    atao = -1.0e-2*(2.26*np.power(self._prs,2)+2.07*self._prs+2.3)*np.power(vtao,tmp)
    tmp  = 0.0274*np.power(self._prs,2)+0.0529*self._prs-0.8765
    btao = (0.0573*np.power(self._prs,2)+0.0937*self._prs+0.202)*np.power(vtao,tmp)
    tmp  = 0.01867*np.power(self._prs,2)+0.4011*self._prs-1.8186
    ctao = 1.0e-4*(9.654*np.power(self._prs,2)+4.945*self._prs+3.1)*np.power(vtao,tmp)
    print(1.0e-4,10e-4)
    #stao = 0.0
    for i in range(n1):
      #n2 = 20.0 - 34.0 *self._phi[i] +14.0*np.power(phi[i],2) 
      n2 = 20.0 - 34.0 *self._phiC +14.0*np.power(self._phiC,2) 
      n2*= self._c
      stao = atao*np.power(alpha[i],2.0)+btao*alpha[i]+ctao
      
      g[i] = (3.0/5.0)*self._kdryC[i]+3.0*n2*(1.0-self._phiC)*self._gc*stao/20.0
    return g

  def ConstantCement(self,kb,gb,phi,khm,ghm,nconstant):
    self._kb = kb 
    self._gb = gb
    self._phi = phi 
    self._khm = khm
    self._ghm = ghm
    self._nconstant = nconstant
    n1 = len(phi)
    n2 = self._nconstant
    k = np.zeros((n2,n1))
    g = np.zeros((n2,n1))
    k_perc = np.zeros(n2)
    g_perc = np.zeros(n2)
    tmp1 = 0
    for i2 in range(n2):
      tmp1 = tmp1 + 1
      k_perc[i2]= ((self._kb[n1-tmp1]-self._khm)*100)/(self._kb[0])
      g_perc[i2]= ((self._gb[n1-tmp1]-self._ghm)*100)/(self._gb[0])
      for i1 in range(n1):
        if i1 <(n1-tmp1):
          k[i2][i1] = (self._phi[i1]/self._phi[n1-tmp1])/(self._kb[n1-tmp1]+4.0/3.0*self._gb[n1-tmp1])+ \
                          (1.0-self._phi[i1]/self._phi[n1-tmp1])/(self._ks+4.0/3.0*self._gb[n1-tmp1])
          k[i2][i1] = 1.0/k[i2][i1]-4.0/3.0*self._gb[n1-tmp1]
          zz = self._gb[n1-tmp1]/6.0*(9.0*self._kb[n1-tmp1]+8.0*self._gb[n1-tmp1])/ \
               (self._kb[n1-tmp1]+2.0*self._gb[n1-tmp1])
          g[i2][i1] = (self._phi[i1]/self._phi[n1-tmp1])/(self._gb[n1-tmp1]+zz )+ \
                      (1.0-self._phi[i1]/self._phi[n1-tmp1])/(zz+self._gs)
          g[i2][i1] = 1.0/g[i2][i1]-zz                            
        else: 
          k[i2][i1] = self._kb[i1]
          g[i2][i1] = self._gb[i1]
           
    return (k,g,k_perc,g_perc)



