#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import elastic_parameters as ep



'''
Elastic parameters
'''

class Parameters:
  def __init__(self,kdry_l,ksat_l,gdry_l,rhosat_l,ksat_loil,rhosat_loil,phi,co2,p):
    self._ksat_l = ksat_l
    self._ksat_loil = ksat_loil
    self._kdry_l = kdry_l
    self._gdry_l = gdry_l
    self._rhosat_l = rhosat_l
    self._rhosat_loil = rhosat_loil

    self._co2 = co2
    self._p = p
    self._phi = phi


  def calculation (self,parameters):
    par = {}
    n = len(parameters)
    n1 = len(self._phi)
    n2 = len(self._co2)
    n3 = len(self._p)
    for i in range(n):
      par[parameters[i]] =  np.zeros(((n3,n2,n1)))
  
    
    vpBrineCo2 = ep.vp(self._ksat_l,self._gdry_l,self._rhosat_l,self._phi,self._co2,self._p)   
    vsBrineCo2 = ep.vs(self._gdry_l,self._rhosat_l,self._phi,self._co2,self._p)
    ipBrineCo2 = ep.Ip(vpBrineCo2,self._rhosat_l,self._phi,self._co2,self._p) 
    isBrineCo2 = ep.Ip(vsBrineCo2,self._rhosat_l,self._phi,self._co2,self._p) 
    rtBrineCo2 = ep.vp_vs(vpBrineCo2,vsBrineCo2,self._phi,self._co2,self._p) 

    vpBrineOil = ep.vp(self._ksat_loil,self._gdry_l,self._rhosat_loil,self._phi,self._co2,self._p)   
    vsBrineOil = ep.vs(self._gdry_l,self._rhosat_loil,self._phi,self._co2,self._p)
    ipBrineOil = ep.Ip(vpBrineOil,self._rhosat_loil,self._phi,self._co2,self._p) 
    isBrineOil = ep.Ip(vsBrineOil,self._rhosat_loil,self._phi,self._co2,self._p) 
    rtBrineOil = ep.vp_vs(vpBrineOil,vsBrineOil,self._phi,self._co2,self._p) 


    par['vpBrineCo2'] = vpBrineCo2        
    par['vsBrineCo2'] = vsBrineCo2
    par['ipBrineCo2'] = ipBrineCo2
    par['isBrineCo2'] = isBrineCo2
    par['rtBrineCo2'] = rtBrineCo2

    par['vpBrineOil'] = vpBrineOil        
    par['vsBrineOil'] = vsBrineOil
    par['ipBrineOil'] = ipBrineOil
    par['isBrineOil'] = isBrineOil
    par['rtBrineOil'] = rtBrineOil

    return par

  
  def differences (self,property,pref,co2ref,type):
    n1 = len(self._phi)
    n2 = len(self._co2)
    n3 = len(self._p)
    property_dif = np.zeros(((n3,n2,n1)))
    #dif = np.zeros(n2)
   
    if type == 'P':
      dif = np.zeros(n3)
      for i3 in range (n3):
        dif[i3]=self._p[i3]-self._p[pref]
        for i2 in range (n2):
          for i1 in range (n1):
            property_dif[i3][i2][i1]=(property[i3][i2][i1]-property[pref][i2][i1])/property[pref][i2][i1]*100
            #print (property[i3][i2][i1],i1,i2,i3)
    elif type == 'S':
      dif = np.zeros(n2)
      for i3 in range (n3):
        for i2 in range (n2):
          dif[i2]=self._co2[i2]-self._co2[co2ref]
          for i1 in range (n1):
            property_dif[i3][i2][i1]=(property[i3][i2][i1]-property[i3][co2ref][i1])/property[i3][co2ref][i1]*100
    elif type == 'B':
      dif = np.zeros(n2)
      for i3 in range (n3):
        for i2 in range (n2):
          dif[i2]=self._co2[i2]-self._co2[co2ref]
          for i1 in range (n1):
            property_dif[i3][i2][i1]=(property[i3][i2][i1]-property[pref][co2ref][i1])/property[pref][co2ref][i1]*100

    else:
      print('Error in the selection, you must type pressure or co2')
    return (property_dif,dif)
         



