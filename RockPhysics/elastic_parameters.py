#!/usr/bin/env python
import numpy as np
import math


from pylab import *
import numpy as np
import matplotlib.pyplot as plt



'''
Elastic parameters
'''


def vp_vs (vp,vs,phi,co2,p):
    '''
    Vp/Vs ratio,
    output is a vector[n]
    Vp[n]: input vector with Vp velocity
    Vs[n]: input vector with Vs velocity
    phi[n]: input vector with porosity
    '''
    n1 = len(phi)
    n2 = len(co2)
    n3 = len(p)
    vp_vs = np.zeros(((n3,n2,n1)))
    for i3 in range (n3):
        for i2 in range(n2):
          for i1 in range(n1):
              vp_vs[i3][i2][i1] = vp[i3][i2][i1]/vs[i3][i2][i1]
    return vp_vs




def Ip (vp,rho,phi,co2,p):
    '''
    Acoustic impedance,
    output is a vector[n]
    Vp[n]: input vector with Velocity
    rho[n]: input vecotr with density
    phi[n]: input vector with porosity

    '''
    n1 = len(phi)
    n2 = len(co2)
    n3 = len(p)
    Ip = np.zeros(((n3,n2,n1)))
    for i3 in range (n3):
        for i2 in range(n2):
          for i1 in range(n1):
              Ip[i3][i2][i1] = vp[i3][i2][i1]*rho[i3][i2][i1]
    return Ip


def poisson_value (K,G):
    '''
    Poisson ratio,
    output is a vector[n]
    K[n]: input vector with Bulk's modulus
    G[n]: input vecotr with Shear modulus
    phi[n]: input vector with porosity

    '''
   
    poisson = (3*K -2*G) /(2*(3*K+G))

    return poisson



def vs (G,rho,phi,co2,p):
    '''
    Shear velocity,
    output is a vector[n]
    rho[n]: input vecotr with density
    G[n]: input vector with Shear modulus
    phi[n]: input vector with porosity

    '''
    n1 = len(phi)
    n2 = len(co2)
    n3 = len(p)
    vs = np.zeros(((n3,n2,n1)))
    for i3 in range (n3):
        for i2 in range(n2):
          for i1 in range(n1):
            vs[i3][i2][i1] = np.power((G[i3][i1]/rho[i3][i2][i1]),0.5)
    return vs

def vp (K,G,rho,phi,co2,p):
    '''
    Compresional velocity,
    output is a vector[n]
    K[n]: input vector with Bulk's modulus
    G[n]: input vecotr with Shear modulus
    rho[n]: input vecotr with density
    phi[n]: input vector with porosity
    '''
    n1 = len(phi)
    n2 = len(co2)
    n3 = len(p)
    vp = np.zeros(((n3,n2,n1)))
    for i3 in range (n3):
        for i2 in range(n2):
            for i1 in range(n1):
                tmp1  =(K[i3][i2][i1]+4.0/3.0*G[i3][i1])/rho[i3][i2][i1]
                vp[i3][i2][i1] = np.power(tmp1,0.5)

    return vp


def M_value (K,G):
    '''
    Shear moduli,
    output is a value
    K: input value with Bulk's modulus
    G: input value with Shear modulus
    '''
    M_value =  K + (4.0/3.0)*G

    return M_value





