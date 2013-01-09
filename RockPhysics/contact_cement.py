#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt





# The contact-cement model
def Kdry_cement(Kc,Gc,Ks,PRc,PRs,phi_c,Mc,phi,Gs,P):
    
    Nk=P.size
    Np=phi.size

    alpha = np.zeros(Np)

    for i in range(Np):
        alpha[i]= (2.0/3.0)*(phi_c-phi[i])/(1-phi_c)
        alpha[i]= np.power(alpha[i], 0.5)
 
    Kdry_c = np.zeros(Np)

    
   
    Vn  = 2*Gc*(1-PRs)*(1-PRc)
    Vn /= math.pi*Gs*(1-2*PRc)
        
    An=-.0245153*np.power(Vn,-1.3646)
    Bn=0.20405*np.power(Vn,-0.89008)
    Cn=0.00024649*np.power(Vn,-1.9864)


    for j in range(Np):
        n2= 20.0 - 34.0 *phi[j] +14.0*np.power(phi[j],2) 
        n2 *= 1.3

        Sn= An*np.power(alpha[j],2.0) + Bn*alpha[j] + Cn
        Kdry_c[j]=n2*(1-phi_c)*Mc*Sn/6

    return Kdry_c


def Gdry_cement(phi,phi_c,Gc,Gs,PRs,Kdry_c,Ks):

    Np=phi.size

    alpha = np.zeros(Np)

    for i in range(Np):
        alpha[i]= (2.0/3.0)*(phi_c-phi[i])/(1-phi_c)
        alpha[i]= np.power(alpha[i], 0.5)
 
    Gdry_c = np.zeros(Np)
    
    
    Vtao  = Gc/(math.pi*Gs)
    tmp=0.079*np.power(PRs,2)+0.1754*PRs-1.342 
    Atao=-1e-2*(2.26*np.power(PRs,2)+2.07*PRs+2.3)*np.power(Vtao,tmp)
    tmp=0.0274*np.power(PRs,2)+0.0529*PRs-0.8765
    Btao=(0.0573*np.power(PRs,2)+0.0937*PRs+0.202)*np.power(Vtao,tmp)
    tmp=0.01867*np.power(PRs,2)+0.4011*PRs-1.8186
    Ctao=1e-4*(9.654*np.power(PRs,2)+4.945*PRs+3.1)*np.power(Vtao,tmp)

    Stao=0.0
    for j in range(Np):
        n2= 20.0 - 34.0 *phi[j] +14.0*np.power(phi[j],2) 
        n2 *=1.3
        Stao= Atao*np.power(alpha[j],2.0) + Btao*alpha[j] + Ctao
        Gdry_c[j]=(3.0/5.0)*Kdry_c[j] + 3.0*n2*(1.0-phi_c)*Gc *Stao/20.0
    return Gdry_c



