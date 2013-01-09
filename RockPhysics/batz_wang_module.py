#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt



'''
Elastic models based on North Sea geology,
the modeled rock has two components:
quartz and clay
'''


def Gassman(Kdry,Ks,phi,Kwater,Rhos,Kco2,rhowater,rhoco2):

    '''
    K  sat calculation by using Gassman
    This modules works for combining Bounds
    with fluid effect
    '''
    
    nvol = phi.size
    CO2 =np.arange(0.0,1.0,0.1)
    nCO2=CO2.size
    Ksat = np.zeros((nvol,nCO2)) 
    rhosat = np.zeros((nvol,nCO2)) 
    phi[0]=0.00000000000001
    for i in range(nvol):
        tmp3=0.0
        for j in range (nCO2):
            Kfluid= (1.0-CO2[j])/Kwater  + CO2[j]/Kco2
            Kfluid=1.0/Kfluid
            rhofluid=(1.0-CO2[j])*rhowater + CO2[j]*rhoco2
            
            tmp3=(phi[i]/Kfluid+(1.0-phi[i])/Ks-Kdry[i]/(Ks*Ks))
            tmp1=((1.0-Kdry[i]/Ks)*(1.0-Kdry[i]/Ks))
            Ksat[i][j]=  Kdry[i] +tmp1/tmp3 
            rhosat[i][j]= rhofluid*phi[i] + (1.0-phi[i])*Rhos 
   
    return (Ksat,rhosat) 




