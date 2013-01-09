#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt




def K_hertz_mindlin ( n, phi_c, Gs,PRs, P):
    '''
    Hertz-Mindlin:
    output is a vector[n]
    phi_c: is an scalar with critical porisity
    Gs[n]: vector with Shear modulus from min. endpoint
    P    : effective pressure
    PRs  : poisson ratio
    '''
    nvol = P.size

    khm = np.zeros(nvol)
    po=1.0/3.0

    for i in range(nvol):
        khm[i] = (n*n*np.power(1-phi_c,2)*Gs*Gs)

        khm[i] /= (18*np.power(np.pi,2)*np.power(1-PRs,2))
        
        khm[i] *=P[i] 
        khm[i] = np.power(khm[i],po)

    return khm



def G_hertz_mindlin ( n, phi_c, Gs,PRs, P,Ft):
    '''
    Hertz-Mindlin:
    output is a vector[n]
    phi_c: is an scalar with critical porisity
    Gs[n]: vector with Shear modulus from min. endpoint
    P    : effective pressure
    PRs  : poisson ratio
    '''
    nvol = P.size

    ghm = np.zeros(nvol)
    ghm_slip = np.zeros(nvol)

    for i in range(nvol):
      #  Without tanget slip effect
        ghm[i] = (3.0*n*n*np.power(1.0-phi_c,2.0)*Gs*Gs)
        ghm[i] /= (2.0*np.power(np.pi,2.0)*np.power(1.0-PRs,2.0))
        ghm[i] *=P[i] 
        ghm[i] =((5.0-4.0*PRs)/(10-5*PRs)) *np.power(ghm[i],0.33333333333333333333333)
        tmp=n*n*P[i]*1.0/10.0*12.0*np.power(1.0-phi_c,2.0)*Gs*Gs/(3.1416*3.1416*np.power(1-PRs,2))
        ghm_slip[i] = np.power(tmp,0.33333333)
        tmp = 3.0/10.0*12.0*np.power(1.0-phi_c,2.0)*Gs*Gs*(1.0-PRs)/(3.141516*3.141516*np.power(2.0-PRs,3))
        ghm_slip[i]  = ghm_slip[i] + np.power(tmp*n*n*P[i],0.333333)*Ft  #Ft slip factor

    
 
    return (ghm,ghm_slip)


def Z_bound(Khm,Ghm,P):
    nvol = P.size

    Zhm = np.zeros(nvol)
    for i in range(nvol):
        Zhm[i]=Ghm[i]/6.0*(9*Khm[i]+8*Ghm[i])
        Zhm[i] /=Khm[i]+2*Ghm[i]  
    return Zhm

def Z_value(Ks,Gs,P):
    Z=Gs/6.0*(9*Ks+8*Gs)
    Z /=Ks+2*Gs  
    return Z


