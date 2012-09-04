#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt


def K_constant_cement(Ks,Kb,Gb,Gs,phi,Khm,Ghm):

    Np=phi.size
    Nconstant=50
    tmp1 = 1
    Kconstant = np.zeros((Nconstant,Np))
    Gconstant = np.zeros((Nconstant,Np))
    Kconstant_perc = np.zeros((Nconstant))
    Gconstant_perc = np.zeros((Nconstant))


    for k in range(Nconstant):

        tmp1= tmp1 + 1
        phi_plot = np.zeros(Np-tmp1)
        Kconstant_perc[k]= ((Kb[Np-tmp1]-Khm[0]) *100)/(Kb[0])
        Gconstant_perc[k]= ((Gb[Np-tmp1]-Ghm[0]) *100)/(Gb[0])

        for j in range(Np):

            if j <(Np-tmp1):

               Kconstant[k][j]= (phi[j]/phi[Np-tmp1])/(Kb[Np-tmp1] +4.0/3.0*Gb[Np-tmp1])+ \
                    (1.0-phi[j]/phi[Np-tmp1])/(Ks +4.0/3.0*Gb[Np-tmp1])
         

               Kconstant[k][j]= 1.0/Kconstant[k][j] - 4.0/3.0*Gb[Np-tmp1]
               ZZ=Gb[Np-tmp1]/6.0*(9*Kb[Np-tmp1]+8*Gb[Np-tmp1])/(Kb[Np-tmp1]+2*Gb[Np-tmp1])
               Gconstant[k][j]= (phi[j]/phi[Np-tmp1])/(Gb[Np-tmp1]+ZZ )+ \
                    (1.0-phi[j]/phi[Np-tmp1])/(ZZ+Gs)
               Gconstant[k][j]= 1.0/Gconstant[k][j] - ZZ                
                
            else: 
               Kconstant[k][j]=Kb[j]
               Gconstant[k][j]=Gb[j]
           

    return (Kconstant,Gconstant,Kconstant_perc,Gconstant_perc)



