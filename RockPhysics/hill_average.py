#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt


'''
Hill's average method is simply the average
between Voigt and Reuss
'''



def Ksget (Kqz, Kclay, Ksid, clayVol,sidVol):
    '''
    Mineral endpoint
    output[n]: combined modulus
    with same dimension as clayVol

    Kqz: scalar with qz B modulus
    Kclay: scalar with clay modulus
    clayVol[n]: vector with volume 
    content
    '''
    
    qzVol = 1.0
    qzVol = qzVol - clayVol - sidVol
    tmp=0.0
    

    tmp  = clayVol * Kclay
    tmp += qzVol   * Kqz    
    tmp += sidVol   * Ksid    
   

    tmp2  = clayVol / Kclay
    tmp2 += qzVol   / Kqz       
    tmp2 += sidVol   / Ksid   
        
    Ks = 0.5*(tmp + 1/tmp2)


    return Ks

def Gsget (Gqz,Gclay, clayVo,sidVoll):
    '''
    Mineral endpoint
    output[n]: combined modulus
    with same dimension as clayVol

    Gqz: scalar with qz B modulus
    Gclay: scalar with clay modulus
    clayVol[n]: vector with volume 
    content
    '''
    Gs= Ksget(Gqz,Gclay,Gsid,clayVol,sidVol)
    return Gs

