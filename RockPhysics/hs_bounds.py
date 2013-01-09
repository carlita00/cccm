#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt


# Modified upper Hashin-Shtrikman bound
def K_sashin_shtrikman_upper(Ks,Khm,Gs,phi_c,phi,P):

    Nk=P.size
    Np=phi.size
    
    Kdry = np.zeros((Nk,Np))

    
    for i in range(Nk):
        for j in range(Np):
            Kdry[i][j]= (phi[j]/phi_c)/(Khm[i] +4.0/3.0*Gs[i])+ \
                       (1-phi[j]/phi_c)/(Ks +4.0/3.0*Gs[i])


            Kdry[i][j]= 1/Kdry[i][j] - 4.0/3.0*Gs[i]

    return Kdry
     
        
def G_sashin_shtrikman_upper(Ghm,Gs,phi_c,phi,Z,Ks,P):
    Nk=P.size
    Np=phi.size
    
    Gdry = np.zeros((Nk,Np))

    
    for i in range(Nk):
        for j in range(Np):
            Gdry[i][j]= (phi[j]/phi_c)/(Ghm[i] +Z[i])+ \
                       (1-phi[j]/phi_c)/(Gs+Z[i])

            Gdry[i][j]= 1/Gdry[i][j] - Z[i]

    return Gdry


# Modified upper Hashin-Shtrikman bound 
def K_sashin_shtrikman_upper2(Ks,Khm,Gs,phi_c,phi,P):

    Nk=P.size
    Np=phi.size
    
    Kdry = np.zeros((Nk,Np))

    
    for i in range(Nk):
        for j in range(Np):
            Kdry[i][j]= (phi[j]/phi_c)/(Khm[i] +4.0/3.0*Gs)+ \
                       (1-phi[j]/phi_c)/(Ks +4.0/3.0*Gs)


            Kdry[i][j]= 1/Kdry[i][j] - 4.0/3.0*Gs

    return Kdry
     
        


def G_sashin_shtrikman_upper2(Ghm,Gs,phi_c,phi,Z,Ks,P):
    Nk=P.size
    Np=phi.size
    
    Gdry = np.zeros((Nk,Np))

    
    for i in range(Nk):
        for j in range(Np):
            Gdry[i][j]= (phi[j]/phi_c)/(Ghm[i] +Z)+ \
                       (1-phi[j]/phi_c)/(Gs+Z)

            Gdry[i][j]= 1/Gdry[i][j] - Z

    return Gdry

            

# Modified lower Hashin-Shtrikman bound
def K_sashin_shtrikman_lower(Ks,Khm,Ghm,phi_c,phi,P):
    Kdry=  K_sashin_shtrikman_upper(Ks,Khm,Ghm,phi_c,phi,P)
    return Kdry

def G_sashin_shtrikman_lower(Ghm,Gs,phi_c,phi,Z,Ks,P):
    Gdry= G_sashin_shtrikman_upper(Ghm,Gs,phi_c,phi,Z,Ks,P)
    return Gdry


