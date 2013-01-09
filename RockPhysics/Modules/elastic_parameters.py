#!/usr/bin/env python



'''
Elastic parameters
'''


def vp_vs (vp,vs,phi):
    '''
    Vp/Vs ratio,
    output is a vector[n]
    Vp[n]: input vector with Vp velocity
    Vs[n]: input vector with Vs velocity
    phi[n]: input vector with porosity
    '''
    nvol=phi.size
    vp_vs = np.zeros(nvol)
    for i in range (nvol):
        vp_vs[i] = vp[i]/vs[i]
    return vp_vs




def Ip (vp,rho,phi):
    '''
    Acoustic impedance,
    output is a vector[n]
    Vp[n]: input vector with Velocity
    rho[n]: input vecotr with density
    phi[n]: input vector with porosity

    '''
    nvol=phi.size
    Ip = np.zeros(nvol)
    for i in range (nvol):
        Ip[i] = vp[i]*rho[i]
    return Ip


def poisson (K,G ):
    '''
    Poisson ratio,
    output is a vector[n]
    K[n]: input vector with Bulk's modulus
    G[n]: input vecotr with Shear modulus
    phi[n]: input vector with porosity

    '''
   
    poisson = (3*K -2*G) /(2*(3*K+G))

    return poisson



def vs (G,rho,phi):
    '''
    Shear velocity,
    output is a vector[n]
    rho[n]: input vecotr with density
    G[n]: input vector with Shear modulus
    phi[n]: input vector with porosity

    '''
    nvol=phi.size
    vs = np.zeros(nvol)
    
    tmp = np.zeros(nvol)
    for i in range (nvol):
       # vs=G[i]
        vs[i] = np.power((G[i]/rho[i]),0.5)

    return vs

def vp (K,G,rho,phi):
    '''
    Compresional velocity,
    output is a vector[n]
    K[n]: input vector with Bulk's modulus
    G[n]: input vecotr with Shear modulus
    rho[n]: input vecotr with density
    phi[n]: input vector with porosity
    '''
    nvol=phi.size
    vp = np.zeros(nvol)
    
    tmp = np.zeros(nvol)
    for i in range (nvol):
        tmp1  =(K[i]+4.0/3.0*G[i])/rho[i]
        vp[i] = np.power(tmp1,0.5)

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





