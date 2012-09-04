#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt


def Differences(Property,P,phi,Po):
    Np=phi.size
    Nk=P.size
    CO2 =np.arange(0.0,1.0,0.1)
    nCO2=CO2.size

    for i in range (Nk):

        if P[i]==Po:
           ipo=i
        else:
           i=0

    dif = np.zeros((Nk,nCO2,Np))
    dif_pressure = np.zeros((Nk))
   
    for i in range (Nk):
        dif_pressure[i]=P[ipo]-P[i]
        
        for j in range (nCO2):
            for k in range (Np):

         #       dif[i][j][k]=(Property[i,j,k]-Property[Nk-1,j,k])/Property[Nk-1,j,k]*100
                dif[i][j][k]=(Property[i,j,k]-Property[ipo,j,k])/Property[ipo,j,k]*100

    return (dif,dif_pressure)

def Differences_sat(Property,P,phi,Po,ico2):
    Np=phi.size
    Nk=P.size
    CO2 =np.arange(0.0,1.0,0.1)
    nCO2=CO2.size
    

    dif = np.zeros((Nk,nCO2,Np))
    dif_pressure = np.zeros((nCO2))
   
    for i in range (Nk):


        for j in range (Nk):
            dif_pressure[j]=CO2[j]-CO2[0]

            for k in range (Np):

                dif[i][j][k]=(Property[i,j,k]-Property[i,ico2,k])/Property[i,ico2,k]*100

    return (dif,dif_pressure)


def Differences_both(Property,P,phi,Po,ico2):
    Np=phi.size
    Nk=P.size
    CO2 =np.arange(0.0,1.0,0.1)
    nCO2=CO2.size
   
    for i in range (Nk):
 
        if P[i]==Po:
           ipo=i
        else:
           i=0


    dif_both = np.zeros((Nk,nCO2,Np))
    dif_pressure = np.zeros((nCO2))
   
    for i in range (Nk):

        if P[i]==Po:
           ipo=i
        else:
           i=0
 
 
    for i in range (Nk):
        for j in range (nCO2):
            dif_pressure[j]=CO2[j]-CO2[0]

            for k in range (Np):

           #     dif_both[i][j][k]=(Property[Nk-1-i,j,k]-Property[Nk-1,0,k])/Property[Nk-1,0,k]*100
                dif_both[i][j][k]=(Property[i,j,k]-Property[ipo,ico2,k])/Property[ipo,ico2,k]*100

    return (dif_both)


def Differences_test(Property,P,phi,Po,ico2):
    Np=phi.size
    Nk=P.size
    CO2 =np.arange(0.0,1.0,0.1)
    nCO2=CO2.size
   
    for i in range (Nk):
 
        if P[i]==Po:
           ipo=i
        else:
           i=0


    dif_both = np.zeros((Nk,nCO2,Np))
    dif_pressure = np.zeros((nCO2))
   
    for i in range (Nk):

        if P[i]==Po:
           ipo=i
        else:
           i=0

  
 
 
    for i in range (nCO2):
        for j in range (Nk):
            dif_pressure[j]=CO2[j]-CO2[0]

            for k in range (Np):

           #     dif_both[i][j][k]=(Property[Nk-1-i,j,k]-Property[Nk-1,0,k])/Property[Nk-1,0,k]*100
                dif_both[i][j][k]=(Property[i,j,k]-Property[ipo,ico2,k])/Property[ipo,ico2,k]*100

    return (dif_both)


