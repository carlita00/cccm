#!/usr/bin/env python
import numpy as np
import elastic_module as em
import elastic_parameters as ep
import hill_average as ha
import gassman_module as gm
import hertz_mindlin as hm
import hs_bounds as mhs
import contact_cement as cc
import constant_cement as ctc
import matplotlib.pyplot as plt

from pylab import *
import numpy as np
import matplotlib.pyplot as plt

'''
Mineral enpoint
'''
Kqz = 37.0e09   #Pa

Gqz = 44.0e9  #44.0e09   #Pa

Kclay = 25.5e09 #25.0e09 #Pa

Gclay = 9.4e09 #9.0e09  #Pa 

Ksid=123e09

Gsid=51e09

Kc = 25.0e09 # 70.2e09   #Pa

Gc = 5.0e09 #29.e09    #Pa

Rhoqz = 2650    #Pa

Rhoclay = 2550    #Pa

clayVol = 0.1 # 0.12 #= Clay proportion Paluxy con mucho clay

sidVol  = 0.0 # 0.1


Ks  = ha.Ksget(Kqz, Kclay, Ksid, clayVol,sidVol)

Gs  = ha.Ksget(Gqz, Gclay, Ksid, clayVol,sidVol)

Rhos  =  2632 # 2632 for Paluxy dirty em.Ksget(Rhoqz, Rhoclay, clayVol)

PRs = ep.poisson(Ks,Gs)

PRc = ep.poisson(Kc,Gc)

Mc  = ep.M_value(Kc,Gc)


'''
Herts-Mindlin endpoint
'''
phi_c = 0.4001 # Critical porosity, hard coded
#phi_c = 0.36010 # Critical porosity, hard coded

'''
Relation between cordination number and porosity 
input parameter in Herts-Mindlin model and
contact dement model ###VERY IMPORTANT#####
'''

n = 20.0 - 34.0 *phi_c +14.0*np.power(phi_c,2) 
n *=1.2

'''
Slip factor in Hertz-Mindlin model, underconstruction
'''

Ft=0  


'''
Gassman substitution parameters for P[3.5 ... 12.5]MPa
'''

rhoco2  = [0.67*1e3, 0.66*1e3,  0.64*1e3, 0.60*1e3, 0.56*1e3, 0.52*1e3, 0.48*1e3,0.44*1e3,0.33*1e3, 0.35*1e3]

rhobrine= [0.998*1e3,0.99*1e3,  0.99*1e3, 0.99*1e3,0.99*1e3, 0.998*1e3,0.99*1e3,0.99*1e3,0.99*1e3, 0.99*1e3] 

rhooil  = [0.796*1e3,0.79*1e3,  0.79*1e3, 0.79*1e3,0.77*1e3, 0.75*1e3, 0.73*1e3,0.70*1e3,0.63*1e3, 0.64*1e3] 

Kco2    = [0.125*1e9,0.112*1e9, 0.107*1e9,0.08*1e9,0.07*1e9, 0.05*1e9, 0.04*1e9,0.03*1e9,0.02*1e9, 0.019*1e9] 

Kbrine  = [2.53*1e9, 2.50*1e9,  2.52*1e9, 2.51*1e9,2.5*1e9,  2.49*1e9, 2.49*1e9,2.48*1e9,2.47*1e9, 2.47*1e9] 

Koil    = [1.13*1e9, 1.1*1e9,   1.09*1e9, 1.07*1e9,1.07*1e9, 1.05*1e9, 1.04*1e9,1.03*1e9,1.02*1e9, 1.00*1e9]


'''
rhoco2  = [0.6*1e3, 0.6*1e3,  0.6*1e3, 0.6*1e3, 0.6*1e3, 0.6*1e3, 0.6*1e3,0.6*1e3,0.6*1e3, 0.6*1e3]

rhobrine= [0.99*1e3,0.99*1e3,  0.99*1e3, 0.99*1e3,0.99*1e3, 0.99*1e3,0.99*1e3,0.99*1e3,0.99*1e3, 0.99*1e3] 

rhooil  = [0.79*1e3,0.79*1e3,  0.79*1e3, 0.79*1e3,0.79*1e3, 0.79*1e3, 0.79*1e3,0.79*1e3,0.79*1e3, 0.79*1e3] 

Kco2    = [0.125*1e9,0.125*1e9, 0.125*1e9,0.125*1e9,0.125*1e9, 0.125*1e9, 0.125*1e9,0.125*1e9,0.125*1e9, 0.125*1e9] 

Kbrine  = [2.47*1e9, 2.47*1e9,  2.47*1e9, 2.47*1e9,2.47*1e9,  2.47*1e9, 2.47*1e9,2.47*1e9,2.47*1e9, 2.47*1e9] 

Koil    = [1.*1e9, 1.*1e9,   1.*1e9, 1.*1e9,1.*1e9, 1.*1e9, 1.*1e9,1.*1e9,1.0*1e9, 1.00*1e9]

'''


CO2 =np.arange(0.0,1.0,0.1)

nCO2=CO2.size

'''
Pressure modeling parameters 
'''
Po = 12.5e6                                 #Reference pressure

#P = np.arange( 3.5*1e6, 13.5*1e6, 1.0*1e6 ) #Effective pressure

P = np.arange( 8.5*1e6, 18.5*1e6, 1.0*1e6 ) 

ipo=4

ico2=0

ibrine=9

print(P)
NP=P.size

phi= np.arange(0.0001,0.4001000050,0.001)

#phi= np.arange(0.0001,0.3601000050,0.001)

Np=phi.size
#------------------Elastic modeling------------------------#
'''Hertz mindlin'''

Khm = hm.K_hertz_mindlin(n, phi_c,Gs, PRs, P)
(Ghm,Ghm_slip) = hm.G_hertz_mindlin(n, phi_c,Gs, PRs, P,Ft)
Zhm = hm.Z_bound(Khm,Ghm,P)
Zhm_slip = hm.Z_bound(Khm,Ghm_slip,P)


'''
Modified lower Hashin-Shtrikman bound
'''

Kdry_l = mhs.K_sashin_shtrikman_lower(Ks,Khm,Ghm,phi_c,phi,P)
Gdry_l = mhs.G_sashin_shtrikman_lower(Ghm,Gs,phi_c,phi,Zhm,Ks,P)

'''
Modified upper Hashin-Shtrikman bound
'''
Z_u = hm.Z_value(Ks,Gs,P)
Kdry_u = mhs.K_sashin_shtrikman_upper2(Ks,Khm,Gs,phi_c,phi,P)
Gdry_u = mhs.G_sashin_shtrikman_upper2(Ghm,Gs,phi_c,phi,Z_u,Ks,P)


'''
Upper bound contact cement
'''

Kdry_cement = cc.Kdry_cement(Kc,Gc,Ks,PRc,PRs,phi_c,Mc,phi,Gs,P)

Gdry_cement = cc.Gdry_cement(phi,phi_c,Gc,Gs,PRs,Kdry_cement,Ks)


#--------------K Sat (Gassman) low HS bound-----------------------------#


Ksat_l1 = np.zeros(((NP,Np,nCO2)))

rhosat_l1 = np.zeros(((NP,Np,nCO2)))

#colors=['k','b','r','m','y','g','c','b']

for i in range(NP):
    (Ksat_l,rhosat_l)= gm.Gassman(Kdry_l[i,:],Ks,phi, Kbrine[i],Rhos,Kco2[i],rhobrine[i],rhoco2[i])
    for j in range(Np):     
        for k in range(nCO2): 
            Ksat_l1[i][j][k]=Ksat_l[j][k]
            rhosat_l1[i][j][k]=rhosat_l[j][k]


(Ksat_c,rhosat)= gm.Gassman(Kdry_cement[:],Ks,phi,Kbrine[0],Rhos,Kco2[0],rhobrine[0],rhoco2[0])


#-----------------Selecting the best trend iterative process--------------#

Nconstant=50 #Number of constant cement lines

Ksat_constant = np.zeros(((NP,Nconstant,Np)))

Gsat_constant = np.zeros(((NP,Nconstant,Np)))

Kconstant_perc = np.zeros(((NP,Nconstant)))

Kdry_constant = np.zeros(((NP,Np)))

Kdry_constant_plot = np.zeros(((NP,Nconstant,Np)))

Gdry_constant_plot = np.zeros(((NP,Nconstant,Np)))

KW = np.zeros(((Np))) #Weight for cement model

GW = np.zeros(((Np))) #Weight for cement model

Gdry_constant = np.zeros(((NP,Np)))

nmodel=4  # 4 6 5% de cement
Gnmodel=4 # 3 6 5%de cement 

Kdry_stiff = np.zeros(((NP,Np)))

Gdry_stiff = np.zeros(((NP,Np)))

'''
Selecting the best trend iterative process
'''

for i in range(NP):
    (Kdry_constant1,Gdry_constant1,Kconstant_perc2,Gconstant_perc2)=ctc.K_constant_cement(Ks,Kdry_cement[:],Gdry_cement[:],Gs,phi,Khm,Ghm)

    for k in range(Np): 
        Kdry_constant[i][k]=Kdry_constant1[nmodel][k]
        Gdry_constant[i][k]=Gdry_constant1[Gnmodel][k]



for i in range(NP):
    (Kdry_constant1,Gdry_constant1,Kconstant_perc2,Gconstant_perc2)=ctc.K_constant_cement(Ks,Kdry_cement[:],Gdry_cement[:],Gs,phi,Khm,Ghm)
    for j in range(Nconstant):
  #      print('Kdry percentaje',Kconstant_perc2)
        for k in range(Np): 
            Kdry_constant_plot[i][j][k]=Kdry_constant1[j][k]
            Gdry_constant_plot[i][j][k]=Gdry_constant1[j][k]





'''
Selecting the model corresponding to cement 10%
'''

nmodel=27

Gnmodel=27 

for i in range(NP):
    (Kdry_constant1,Gdry_constant1,Kconstant_perc2,Gconstant_perc2)=ctc.K_constant_cement(Ks,Kdry_cement[:],Gdry_cement[:],Gs,phi,Khm,Ghm)

    for k in range(Np): 
        Kdry_stiff[i][k]=Kdry_constant1[nmodel][k]
        Gdry_stiff[i][k]=Gdry_constant1[Gnmodel][k]




'''
Pressure dependent modeling 
'''



for i in range(NP):
    if P[i]==Po:
       Kdry_constant=Kdry_constant
       Gdry_constant=Gdry_constant

    else:
       for k in range(Np): 
           KW[k]=(Kdry_constant[ipo][k]- Kdry_l[ipo][k])/(Kdry_stiff[ipo][k]- Kdry_l[ipo][k])    
       #   GW[k]=(1.0-KW[k])*Gdry_l[4][k]/Gs + KW[k]*Gdry_stiff[4][k]/Gs
       #   GW[k]=(1.0-KW[k])*Gdry_l[4][k]/Gdry_constant[4][k] + KW[k]*Gdry_stiff[4][k]/Gdry_constant[4][k]

           GW[k]=(Gdry_constant[ipo][k]- Gdry_l[ipo][k])/(Gdry_stiff[ipo][k]- Gdry_l[ipo][k])   
           Kdry_constant[i][k]=(1.0-KW[k])*Kdry_l[i][k]+KW[k]*Kdry_stiff[i][k]
           Gdry_constant[i][k]=(1.0-GW[k])*Gdry_l[i][k]+GW[k]*Gdry_stiff[i][k]

'''
Gassman substitution, to calculate Kconstant brine saturated lines
just for plotting and visualize the results  
'''

for i in range(NP):
    (Ksat_constant1,Gsat_constant1,Kconstant_perc1,Gconstant_perc1)=ctc.K_constant_cement(Ks,Ksat_c[:,0],Gdry_cement,Gs,phi,Khm,Ghm)
    
    for j in range(Nconstant):     
        Kconstant_perc[i][j]=Kconstant_perc1[j]
        for k in range(Np): 
            Ksat_constant[i][j][k]=Ksat_constant1[j][k]
            Gsat_constant[i][j][k]=Gsat_constant1[j][k]
print('paso')

#-------------K Sat once we find the best model---------------------#

Ksat_l_contact   = np.zeros(((NP,Np,nCO2)))

rhosat_l_contact = np.zeros(((NP,Np,nCO2)))

Ksat_l_contact_oil   = np.zeros(((NP,Np,nCO2)))

rhosat_l_contact_oil = np.zeros(((NP,Np,nCO2)))

for i in range(NP):
    (Ksat_l2,rhosat_l2)= gm.Gassman(Kdry_constant[i,:],Ks,phi, Kbrine[i],Rhos,Kco2[i],rhobrine[i],rhoco2[i])
    (Ksat_l3,rhosat_l3)= gm.Gassman(Kdry_constant[i,:],Ks,phi, Kbrine[i],Rhos,Koil[i],rhobrine[i],rhooil[i])

#    (Ksat_l2,rhosat_l2)= gm.Gassman(Kdry_l[i,:],Ks,phi, Kbrine[i],Rhos,Kco2[i],rhobrine[i],rhoco2[i])
#    (Ksat_l3,rhosat_l3)= gm.Gassman(Kdry_l[i,:],Ks,phi, Kbrine[i],Rhos,Koil[i],rhobrine[i],rhooil[i])





#    (Ksat_l2,rhosat_l2)= gm.Gassman(Kdry_l[i,:],Ks,phi, Kbrine[i],Rhos,Kco2[i],rhobrine[i],rhoco2[i])
#    (Ksat_l3,rhosat_l3)= gm.Gassman(Kdry_l[i,:],Ks,phi, Kbrine[i],Rhos,Koil[i],rhobrine[i],rhooil[i])


    for j in range(Np):     
        for k in range(nCO2): 
            Ksat_l_contact[i][j][k]=Ksat_l2[j][k]
            rhosat_l_contact[i][j][k]=rhosat_l2[j][k]
            Ksat_l_contact_oil[i][j][k]=Ksat_l3[j][k]
            rhosat_l_contact_oil[i][j][k]=rhosat_l3[j][k]
print('paso')

#----------Calculating the elastic parameters----------------------#

Vpsat_l   = np.zeros(((NP,nCO2,Np)))

Vssat_l   = np.zeros(((NP,nCO2,Np)))

Vp_Vs_sat_l   = np.zeros(((NP,nCO2,Np)))

Ip_sat_l   = np.zeros(((NP,nCO2,Np)))

Is_sat_l   = np.zeros(((NP,nCO2,Np)))

Vpsat_l_oil   = np.zeros(((NP,nCO2,Np)))

Vssat_l_oil   = np.zeros(((NP,nCO2,Np)))

Vp_Vs_sat_l_oil   = np.zeros(((NP,nCO2,Np)))

Ip_sat_l_oil   = np.zeros(((NP,nCO2,Np)))

Is_sat_l_oil   = np.zeros(((NP,nCO2,Np)))

rho_sat_l   = np.zeros(((NP,nCO2,Np)))

rho_sat_l_oil   = np.zeros(((NP,nCO2,Np)))


for i in range(NP):
    for k in range(nCO2):    
            
        Vssat_l1= ep.vs(Gdry_constant[i,:],rhosat_l_contact[i,:,k],phi)
        Vpsat_l1= ep.vp(Ksat_l_contact[i,:,k],Gdry_constant[i,:],rhosat_l_contact[i,:,k],phi)
        Vp_Vs_sat_l1= ep.vp_vs(Vpsat_l1,Vssat_l1,phi)
        Ip_sat_l1=ep.Ip(Vpsat_l1,rhosat_l_contact[i,:,k],phi) 
        Is_sat_l1=ep.Ip(Vssat_l1,rhosat_l_contact[i,:,k],phi) 


        Vssat_l3= ep.vs(Gdry_constant[i,:],rhosat_l_contact_oil[i,:,k],phi)
        Vpsat_l3= ep.vp(Ksat_l_contact_oil[i,:,k],Gdry_constant[i,:],rhosat_l_contact_oil[i,:,k],phi)
        Vp_Vs_sat_l3= ep.vp_vs(Vpsat_l3,Vssat_l3,phi)
        Ip_sat_l3=ep.Ip(Vpsat_l3,rhosat_l_contact_oil[i,:,k],phi) 
        Is_sat_l3=ep.Ip(Vssat_l3,rhosat_l_contact_oil[i,:,k],phi) 
        
         

        '''
        Vssat_l1= ep.vs(Gdry_l[i,:],rhosat_l_contact[i,:,k],phi)
        Vpsat_l1= ep.vp(Ksat_l_contact[i,:,k],Gdry_l[i,:],rhosat_l_contact[i,:,k],phi)
        Vp_Vs_sat_l1= ep.vp_vs(Vpsat_l1,Vssat_l1,phi)
        Ip_sat_l1=ep.Ip(Vpsat_l1,rhosat_l_contact[i,:,k],phi) 
        Is_sat_l1=ep.Ip(Vssat_l1,rhosat_l_contact[i,:,k],phi) 


        Vssat_l3= ep.vs(Gdry_l[i,:],rhosat_l_contact_oil[i,:,k],phi)
        Vpsat_l3= ep.vp(Ksat_l_contact_oil[i,:,k],Gdry_l[i,:],rhosat_l_contact_oil[i,:,k],phi)
        Vp_Vs_sat_l3= ep.vp_vs(Vpsat_l3,Vssat_l3,phi)
        Ip_sat_l3=ep.Ip(Vpsat_l3,rhosat_l_contact_oil[i,:,k],phi) 
        Is_sat_l3=ep.Ip(Vssat_l3,rhosat_l_contact_oil[i,:,k],phi) 
        '''

       
        for j in range(Np):
            Vpsat_l[i][k][j]= Vpsat_l1[j]
            Vssat_l[i][k][j]= Vssat_l1[j]
            Vp_Vs_sat_l[i][k][j]= Vp_Vs_sat_l1[j]
            Ip_sat_l[i][k][j]= Ip_sat_l1[j]
            Is_sat_l[i][k][j]= Is_sat_l1[j]

            
            Vpsat_l_oil[i][k][j]= Vpsat_l3[j]
            Vssat_l_oil[i][k][j]= Vssat_l3[j]
            Vp_Vs_sat_l_oil[i][k][j]= Vp_Vs_sat_l3[j]
            Ip_sat_l_oil[i][k][j]= Ip_sat_l3[j]
            Is_sat_l_oil[i][k][j]= Is_sat_l3[j]

            rho_sat_l[i][k][j]= rhosat_l_contact[i,j,k]
            rho_sat_l_oil[i][k][j]=rhosat_l_contact_oil[i,j,k]

#--------Calculating percentual differences------------------------

'''
Pressure sensibility
'''

(Vp_co2_dif,dif_pressure)=em.Differences(Vpsat_l[:,:,:],P[:],phi[:],Po)

(Vs_co2_dif,dif_pressure)=em.Differences(Vssat_l[:,:,:],P[:],phi[:],Po)

(Vs_co2_dif,dif_pressure)=em.Differences(Vssat_l[:,:,:],P[:],phi[:],Po)

(VpVs_co2_dif,dif_pressure)=em.Differences(Vp_Vs_sat_l[:,:,:],P[:],phi[:],Po)

(Vp_oil_dif,dif_pressure)=em.Differences(Vpsat_l_oil[:,:,:],P[:],phi[:],Po)

(Vs_oil_dif,dif_pressure)=em.Differences(Vssat_l_oil[:,:,:],P[:],phi[:],Po)

(Vs_oil_dif,dif_pressure)=em.Differences(Vssat_l_oil[:,:,:],P[:],phi[:],Po)

(VpVs_oil_dif,dif_pressure)=em.Differences(Vp_Vs_sat_l_oil[:,:,:],P[:],phi[:],Po)

(rho_oil_dif,dif_pressure)=em.Differences(rho_sat_l_oil[:,:,:],P[:],phi[:],Po)

(rho_co2_dif,dif_pressure)=em.Differences(rho_sat_l[:,:,:],P[:],phi[:],Po)

(Ip_co2_dif,dif_pressure)=em.Differences(Ip_sat_l[:,:,:],P[:],phi[:],Po)

(Ip_oil_dif,dif_pressure)=em.Differences(Ip_sat_l_oil[:,:,:],P[:],phi[:],Po)


(Is_co2_dif,dif_pressure)=em.Differences(Is_sat_l[:,:,:],P[:],phi[:],Po)

(Is_oil_dif,dif_pressure)=em.Differences(Is_sat_l_oil[:,:,:],P[:],phi[:],Po)





'''
Saturation sensibility
'''
(Vp_co2_difs,dif_sat)=em.Differences_sat(Vpsat_l[:,:,:],P[:],phi[:],Po,ico2)

(Vs_co2_difs,dif_sat)=em.Differences_sat(Vssat_l[:,:,:],P[:],phi[:],Po,ico2)

(Vs_co2_difs,dif_sat)=em.Differences_sat(Vssat_l[:,:,:],P[:],phi[:],Po,ico2)

(VpVs_co2_difs,dif_sat)=em.Differences_sat(Vp_Vs_sat_l[:,:,:],P[:],phi[:],Po,ico2)

(Vp_oil_difs,dif_sat)=em.Differences_sat(Vpsat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(Vs_oil_difs,dif_sat)=em.Differences_sat(Vssat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(Vs_oil_difs,dif_sat)=em.Differences_sat(Vssat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(VpVs_oil_difs,dif_sat)=em.Differences_sat(Vp_Vs_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(rho_oil_difs,dif_sat)=em.Differences_sat(rho_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(rho_co2_difs,dif_sat)=em.Differences_sat(rho_sat_l[:,:,:],P[:],phi[:],Po,ico2)

(Ip_co2_difs,dif_sat)=em.Differences_sat(Ip_sat_l[:,:,:],P[:],phi[:],Po,ico2)

(Ip_oil_difs,dif_sat)=em.Differences_sat(Ip_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(Is_co2_difs,dif_sat)=em.Differences_sat(Is_sat_l[:,:,:],P[:],phi[:],Po,ico2)

(Is_oil_difs,dif_sat)=em.Differences_sat(Is_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(Ip_oil_difst)=em.Differences_test(Ip_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)
(VpVs_oil_difst)=em.Differences_test(Vp_Vs_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)
(Ip_co2_difst)=em.Differences_test(Ip_sat_l[:,:,:],P[:],phi[:],Po,ico2)
(VpVs_co2_difst)=em.Differences_test(Vp_Vs_sat_l[:,:,:],P[:],phi[:],Po,ico2)

'''
Saturation and pressure sensibility
'''

(Vp_co2_difps)=em.Differences_both(Vpsat_l[:,:,:],P[:],phi[:],Po,ico2)

(Vs_co2_difps)=em.Differences_both(Vssat_l[:,:,:],P[:],phi[:],Po,ico2)

(Vs_co2_difps)=em.Differences_both(Vssat_l[:,:,:],P[:],phi[:],Po,ico2)

(VpVs_co2_difps)=em.Differences_both(Vp_Vs_sat_l[:,:,:],P[:],phi[:],Po,ico2)

(Vp_oil_difps)=em.Differences_both(Vpsat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(Vs_oil_difps)=em.Differences_both(Vssat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(Vs_oil_difps)=em.Differences_both(Vssat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(VpVs_oil_difps)=em.Differences_both(Vp_Vs_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(rho_oil_difps)=em.Differences_both(rho_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(rho_co2_difps)=em.Differences_both(rho_sat_l[:,:,:],P[:],phi[:],Po,ico2)

(Ip_co2_difps)=em.Differences_both(Ip_sat_l[:,:,:],P[:],phi[:],Po,ico2)

(Ip_oil_difps)=em.Differences_both(Ip_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)

(Is_co2_difps)=em.Differences_both(Is_sat_l[:,:,:],P[:],phi[:],Po,ico2)

(Is_oil_difps)=em.Differences_both(Is_sat_l_oil[:,:,:],P[:],phi[:],Po,ico2)






(Vp_co2_difpsb)=em.Differences_both(Vpsat_l[:,:,:],P[:],phi[:],Po,ibrine)

(Vs_co2_difpsb)=em.Differences_both(Vssat_l[:,:,:],P[:],phi[:],Po,ibrine)

(Vs_co2_difpsb)=em.Differences_both(Vssat_l[:,:,:],P[:],phi[:],Po,ibrine)

(VpVs_co2_difpsb)=em.Differences_both(Vp_Vs_sat_l[:,:,:],P[:],phi[:],Po,ibrine)

(Vp_oil_difpsb)=em.Differences_both(Vpsat_l_oil[:,:,:],P[:],phi[:],Po,ibrine)

(Vs_oil_difpsb)=em.Differences_both(Vssat_l_oil[:,:,:],P[:],phi[:],Po,ibrine)

(Vs_oil_difpsb)=em.Differences_both(Vssat_l_oil[:,:,:],P[:],phi[:],Po,ibrine)

(VpVs_oil_difpsb)=em.Differences_both(Vp_Vs_sat_l_oil[:,:,:],P[:],phi[:],Po,ibrine)

(rho_oil_difpsb)=em.Differences_both(rho_sat_l_oil[:,:,:],P[:],phi[:],Po,ibrine)

(rho_co2_difpsb)=em.Differences_both(rho_sat_l[:,:,:],P[:],phi[:],Po,ibrine)

(Ip_co2_difpsb)=em.Differences_both(Ip_sat_l[:,:,:],P[:],phi[:],Po,ibrine)

(Ip_oil_difpsb)=em.Differences_both(Ip_sat_l_oil[:,:,:],P[:],phi[:],Po,ibrine)




#---------Plotting-------------------------------------------------------

'''
Input files, well-log data
'''
#phi_d,phi_t,Ksat_well,Vsh_iw,G,depth = np.loadtxt('./xplot1_tmp.txt',unpack=True)

phi_d,phi_t,Ksat_well,Vsh_iw,G,SO,depth,well,facies = np.loadtxt('./well_140_159_holt_bryant.txt',unpack=True)



#---------Calculating Kdry in the wells by using Gassman -----------------

Ndepth=depth.size
dry_ratio = np.zeros(Ndepth)
Ks_well = np.zeros(Ndepth)
Ksid_well = np.zeros(Ndepth)
sidVol_well = np.zeros(Ndepth)

for i in range(Ndepth):
    Ks_well[i]  = ha.Ksget(Kqz, Kclay, Ksid, Vsh_iw[i]/100,sidVol_well[i])

Kdry_well= gm.Kdry(Ksat_well,Ks_well,phi_t/100,SO/100,Kbrine[4],Koil[4])


print()

for i in range(Ndepth):    
    dry_ratio[i]=(Kdry_well[i] + 4.0/3.0*G[i])/G[i] 
    dry_ratio[i]=np.power(dry_ratio[i],0.5)


fig = figure(16, figsize=(10, 8))


bx = plt.subplot(1,1,1)
bx.set_xlabel('Vp/Vs ratio')
bx.set_ylabel('Depth (ft)')
p=bx.scatter(dry_ratio,depth,c=Vsh_iw)
cbar= plt.colorbar(p)
cbar.set_label('Vsh')
ylim(3500, 3200) 



#--------------------------------------------------------------------------



fig = figure(12, figsize=(10, 8))

bx = plt.subplot(1,1,1)
bx.set_xlabel('porosity $\phi$')
bx.set_ylabel('Dry Bulk modulus $K$')
p=bx.scatter(phi_t/100,Kdry_well,c=Vsh_iw)
#p=bx.scatter(phi_t/100,Ksat_well,color='b')
cbar= plt.colorbar(p)
cbar.set_label('Vsh')


axis([0, 0.5, 0, 5*1e10])
        
plot(phi[:],Kdry_constant_plot[4,5,:],'b')
plot(phi[:],Kdry_constant_plot[4,13,:],'b')
plot(phi[:],Kdry_constant_plot[4,30,:],'b')
plot(phi[:],Kdry_cement[:],'r')
plot(phi[:],Kdry_l[4,:],'r')

#-----------Plot tusc from Paluxy -----------------------------------------

fig = figure(17, figsize=(10, 8))

bx = plt.subplot(1,1,1)
bx.set_xlabel('porosity $\phi$')
bx.set_ylabel('Dry Bulk modulus $K$')
axis([0, 0.5, 0, 5*1e10])
        
plot(phi[:],Kdry_constant_plot[4,5,:],'b')
plot(phi[:],Kdry_constant_plot[4,13,:],'b')
plot(phi[:],Kdry_constant_plot[4,30,:],'b')
plot(phi[:],Kdry_cement[:],'r')
plot(phi[:],Kdry_l[4,:],'r')

Ndepth=depth.size

for i in range(Ndepth):
    if Vsh_iw[i] <25:
        if facies[i] <10:
            p=bx.scatter(phi_t[i]/100.0,Kdry_well[i],color='r')
        else:
            p=bx.scatter(phi_t[i]/100.0,Kdry_well[i],color='g')
    else:
        k=k


#-----------Plot Facies -----------------------------------------

fig = figure(18, figsize=(10, 8))

bx = plt.subplot(1,1,1)
bx.set_xlabel('porosity $\phi$')
bx.set_ylabel('Dry Bulk modulus $K$')


axis([0, 0.5, 0, 5*1e10])
        
plot(phi[:],Kdry_constant_plot[4,5,:],'b')
plot(phi[:],Kdry_constant_plot[4,13,:],'b')
plot(phi[:],Kdry_constant_plot[4,30,:],'b')
plot(phi[:],Kdry_cement[:],'r')
plot(phi[:],Kdry_l[4,:],'r')

colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.,10.,20.,30.,40.,50.,60.,90.]
for i in range(Ndepth):
    for j in range(len(colors)):
        if Vsh_iw[i] <25:
            if facies[i] == colors[j]:
                c=cm.spectral((j+1)/float(len(colors)),1)
                p=bx.scatter(phi_t[i]/100.0,Kdry_well[i],color=c)
            else:
                k=k
        else:
            k=k

p1 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((1)/float(len(colors)),1))
p2 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((2)/float(len(colors)),1))
p3 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((3)/float(len(colors)),1))
p5 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((5)/float(len(colors)),1))

p20 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((11)/float(len(colors)),1))
p30 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((12)/float(len(colors)),1))
p40 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((13)/float(len(colors)),1))
p50 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((14)/float(len(colors)),1))


labels = ('Tusc: Beach / Barrier bar' , 'Tusc: Washover', 'Tusc: Transitional', \
          'Tusc: Poor rock', 'Paluxy: Distributary sand','Paluxy: mediocre distributary sand', \
          'Paluxy: mediocre distributary sand','Paluxy: Poor rock')
legend = plt.legend([p1,p2,p3,p5,p20,p30,p40,p50],labels, loc=(0.25, .7), labelspacing=0.1)
ltext = gca().get_legend().get_texts()


#-------------------------------------------------------------------------



fig = figure(1, figsize=(10, 8))

bx = plt.subplot(1,1,1)
bx.set_xlabel('porosity $\phi$')
bx.set_ylabel('Bulk modulus $K$')
p=bx.scatter(phi_t/100,Ksat_well,c=Vsh_iw)
cbar= plt.colorbar(p)
cbar.set_label('Vsh')
axis([0.02, 0.5, 0, 5*1e10])
        
plot(phi[:],Ksat_constant[4,5,:],'b')

plot(phi[:],Ksat_constant[4,18,:],'b')
plot(phi[:],Ksat_c[:,0],'r')

plot(phi[:],Ksat_constant[9,36,:],'b')


plot(phi[:],Ksat_l1[9,:,0],'r')

#--------------------------------

fig = figure(2, figsize=(8, 6))
ax = plt.subplot(1,1,1)

phi_t=0.0
phi_d=0.0

ax.set_ylabel('Vp/Vs ratio ')
ax.set_xlabel('Acoustic impedance kg/(m2*s)')
PR,Ip,Vsh_iw,phi_d,phi_t,SO,vp,vs,mu,depth,well,facies= np.loadtxt('./well_140_159_holt_bryant_2b.txt',unpack=True)
Ndepth=depth.size

p=ax.scatter(Ip_sat_l_oil[9,0,:],Vp_Vs_sat_l_oil[9,0,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,1,:],Vp_Vs_sat_l_oil[9,1,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,9,:],Vp_Vs_sat_l_oil[9,9,:],color='g')

p=ax.scatter(Ip_sat_l[9,0,:],Vp_Vs_sat_l[9,0,:],color='r')
p=ax.scatter(Ip_sat_l[9,1,:],Vp_Vs_sat_l[9,1,:],color='r')
p=ax.scatter(Ip_sat_l[9,9,:],Vp_Vs_sat_l[9,9,:],color='r')

ax.set_xlabel('Impedance')
ax.set_ylabel('Vp/Vs ratio')

axis([4*1e6, 12*1e6, 1.0, 3])
        

colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.,10.,20.,30.,40.,50.,60.,90.]
for i in range(Ndepth):
    for j in range(len(colors)):
        if Vsh_iw[i] <25:
            if facies[i] == colors[j]:
                print(facies[i])
                c=cm.spectral((j+1)/float(len(colors)),1)
                p=ax.scatter(Ip[i],PR[i],color=c)
            else:
                k=k
        else:
            k=k

p1 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((1)/float(len(colors)),1))
p2 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((2)/float(len(colors)),1))
p3 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((3)/float(len(colors)),1))
p5 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((5)/float(len(colors)),1))

p20 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((11)/float(len(colors)),1))
p30 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((12)/float(len(colors)),1))
p40 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((13)/float(len(colors)),1))
p50 = Rectangle((0, 0), 0.2, 0.2, fc=cm.spectral((14)/float(len(colors)),1))


labels = ('Tusc: Beach / Barrier bar' , 'Tusc: Washover', 'Tusc: Transitional', \
          'Tusc: Poor rock', 'Paluxy: Distributary sand','Paluxy: mediocre distributary sand', \
          'Paluxy: mediocre distributary sand','Paluxy: Poor rock')
legend = plt.legend([p1,p2,p3,p5,p20,p30,p40,p50],labels, loc=(0.25, .7), labelspacing=0.1)
ltext = gca().get_legend().get_texts()

axis([4*1e6, 12*1e6, 1.0, 3])


fig = figure(3, figsize=(8, 6))
ax = plt.subplot(1,1,1)

#p=ax.scatter(Ip,PR,c=SO)
#cbar= plt.colorbar(p)
#cbar.set_label('Oil')
axis([4*1e6, 12*1e6, 1.0, 3])

ax.set_xlabel('Impedance')
ax.set_ylabel('Vp/Vs ratio')

'''
p=ax.scatter(Ip_sat_l_oil[9,0,:],Vp_Vs_sat_l_oil[9,0,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,1,:],Vp_Vs_sat_l_oil[9,1,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,9,:],Vp_Vs_sat_l_oil[9,9,:],color='g')


p=ax.scatter(Ip_sat_l[9,0,:],Vp_Vs_sat_l[9,0,:],color='r')
p=ax.scatter(Ip_sat_l[9,1,:],Vp_Vs_sat_l[9,1,:],color='r')
p=ax.scatter(Ip_sat_l[9,9,:],Vp_Vs_sat_l[9,9,:],color='r')

'''

Ip_paluxy=np.zeros(Ndepth)
PR_paluxy=np.zeros(Ndepth)
SO_paluxy=np.zeros(Ndepth)
vp_paluxy=np.zeros(Ndepth)
vs_paluxy=np.zeros(Ndepth)
Vsh_iw_paluxy=np.zeros(Ndepth)
phi_t_paluxy=np.zeros(Ndepth)


j=0
for i in range(Ndepth):
    if depth[i] <=968.0:
	#p=bx.scatter(phi_t[i],Ksat[i],color='r')
       k=10
    elif depth[i] >=968.0 and depth[i]<=972.0:
   #    p=ax.scatter(phi_d[i]/100,Ksat_well[i],color='b')
#       p=ax.scatter(Ip[i],PR[i],c=SO[i])

       k=10
    elif depth[i] >=973.0 and depth[i]<=978.0:
       k=10
#       p=ax.scatter(phi_d[i]/100,Ksat_well[i],color='y')
 #      p=ax.scatter(Ip[i],PR[i],c=SO[i])

    elif depth[i] >=979.0 and depth[i]<=984.0:
       k=10
  #     p=ax.scatter(Ip[i],PR[i],c=SO[i])

#       p=ax.scatter(phi_d[i]/100,Ksat_well[i],color='m')
#    elif depth[i]  >=3388.0 and depth[i]<=3446.0:

#       k=10

#    else:
#       k=10


    elif depth[i]  >=3287.0 and depth[i]<=3446.0:

  #     p=ax.scatter(Ip[i],PR[i],c=SO[i])
       if Vsh_iw[i]<1000.0: 
      # print(Vsh_iw[i])
      # print(phi_t[i])

          Ip_paluxy[j]= Ip[i]
          PR_paluxy[j]= PR[i]
          SO_paluxy[j]= SO[i]
          vp_paluxy[j]= vp[i]
          vs_paluxy[j]= vs[i]
          Vsh_iw_paluxy[j]=Vsh_iw[i]
          phi_t_paluxy[j]= phi_t[i]
          j=j+1

       else:
          k=10

 
    else:
       k=10.0

 #   else:
 #      k=10

#p=ax.scatter(Ip_paluxy,PR_paluxy,c=SO_paluxy)
p=ax.scatter(Ip,PR,c=SO)

cbar= plt.colorbar(p)
cbar.set_label('Oil')
axis([4*1e6, 12*1e6, 1.0, 3])

fig = figure(4, figsize=(8, 6))

ax = plt.subplot(1,1,1)
axis([4*1e6, 12*1e6, 1.0, 3])

'''
p=ax.scatter(Ip_sat_l_oil[9,0,:],Vp_Vs_sat_l_oil[9,0,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,9,:],Vp_Vs_sat_l_oil[9,9,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,1,:],Vp_Vs_sat_l_oil[9,1,:],color='g')


p=ax.scatter(Ip_sat_l[9,0,:],Vp_Vs_sat_l[9,0,:,],color='k')
p=ax.scatter(Ip_sat_l[9,9,:],Vp_Vs_sat_l[9,9,:,],color='k')
p=ax.scatter(Ip_sat_l[9,1,:],Vp_Vs_sat_l[9,1,:,],color='k')
p=ax.scatter(Ip_sat_l[9,9,:],Vp_Vs_sat_l[9,9,:,],color='k')
'''
'''
for i in range(Np):
    if phi[i]==.1051:
       p=ax.plot(Ip_sat_l[9,0,i],Vp_Vs_sat_l[9,0,i],'-o', alpha=0.7,ms=20,lw=2,mfc='b')

    elif  phi[i]==.2001:
       p=ax.plot(Ip_sat_l[9,0,i],Vp_Vs_sat_l[9,0,i],'-o', alpha=0.7,ms=20,lw=2,mfc='g')

    elif  phi[i]==.3051:
       p=ax.plot(Ip_sat_l[9,0,i],Vp_Vs_sat_l[9,0,i],'-o', alpha=0.7,ms=20,lw=2,mfc='r')

    else:
       k=1

'''

#p=ax.scatter(Ip,PR,c=phi_d)

p=ax.scatter(Ip,PR,c=Vsh_iw)


#p=ax.scatter(Ip_paluxy,PR_paluxy,c=Vsh_iw_paluxy)

cbar= plt.colorbar(p)

cbar.set_label('Vsh')


ax.set_xlabel('Impedance')
ax.set_ylabel('Vp/Vs ratio')


fig = figure(5, figsize=(6, 5))

ax = plt.subplot(1,1,1)



plot(phi[:],Vpsat_l_oil[0,0,:],color='g')
plot(phi[:],Vpsat_l_oil[9,0,:],color='b')
plot(phi[:],Vpsat_l_oil[4,0,:],color='r')

p=ax.scatter(phi_t,vp,c=SO)
#p=ax.scatter(phi_t_paluxy,vp_paluxy,c=Vsh_iw_paluxy)
#ax.set_vmin
cbar= plt.colorbar(p)
cbar.set_label('SO')
ax.set_xlabel('Porosity')
ax.set_ylabel('Vp (m/s)')


axis([0, 0.37, 1500, 6000])



fig = figure(6, figsize=(6, 5))
ax = plt.subplot(1,1,1)

#plot(phi[:],Vpsat_l_oil[0,0,:],color='g')
plot(phi[:],Vpsat_l_oil[9,0,:],color='b')
#plot(phi[:],Vpsat_l_oil[9,0,:],color='r')

#p=ax.scatter(phi_d,vp,c=Vsh_iw)

p=ax.scatter(phi_t_paluxy,vp_paluxy,c=Vsh_iw_paluxy)

cbar= plt.colorbar(p)
cbar.set_label('Vsh')
ax.set_xlabel('Porosity')
ax.set_ylabel('Vp (m/s)')

axis([0, 0.37, 1500, 6000])


fig = figure(7, figsize=(6, 5))

ax = plt.subplot(1,1,1)

#plot(phi[:],Vssat_l_oil[0,1,:],color='g')
#plot(phi[:],Vssat_l_oil[2,1,:],color='b')
plot(phi[:],Vssat_l_oil[4,0,:],color='b')
p=ax.scatter(phi_t_paluxy,vs_paluxy,c=Vsh_iw_paluxy)
axis([0, 0.37, 500, 4000])

ax.set_xlabel('Porosity')
ax.set_ylabel('Vs (m/s)')

#p=ax.scatter(phi_t,vs,c=SO)
cbar= plt.colorbar(p)
cbar.set_label('Vsh')
ax.set_xlabel('Porosity')
ax.set_ylabel('Vs (m/s)')


'''

fig = figure(8, figsize=(10, 8))
#print(rho)
#print(phi_t)

ax = plt.subplot(1,1,1)

#p=ax.scatter(phi_t,rho,c=SO)
#cbar= plt.colorbar(p)
#cbar.set_label('so')

#plot(rho_sat_l_oil[9,0,:],phi,color='b')
#plot(rho_sat_l_oil[9,1,:],phi,color='g')
#plot(rho_sat_l_oil[9,9,:],phi,color='g')
'''





fig = figure(8, figsize=(10, 8))

plot(phi[:],Gdry_constant[4,:],'b')
plot(phi[:],Gdry_constant[0,:],'c')
plot(phi[:],Gdry_l[0,:],'r')
plot(phi[:],Gdry_l[4,:],'r')

plot(phi[:],Gdry_cement[:],'r')

phi_d,phi_t,Ksat_well,Vsh_iw,G,depth = np.loadtxt('./xplot1_tmp.txt',unpack=True)


ax = plt.subplot(1,1,1)
ax.set_xlabel('porosity $\phi$')
ax.set_ylabel('Shear modulus $K$')
p=ax.scatter(phi_t/100,G,c=Vsh_iw)
cbar= plt.colorbar(p)
cbar.set_label('Vsh')
axis([0, 0.5, 0, 5*1e10])



#rho_co2_dif=em.Differences(rhosat_l[:,:,:],P[:],phi[:])


fig =  figure(9, figsize=(6, 5))

ax = plt.subplot(1,1,1)
ax.set_xlabel('Differential pressure')
ax.set_ylabel('Difference % ')



for i in range (NP-1):
    ax.scatter(dif_pressure[i],VpVs_oil_dif[i,0,320],c='b')
    ax.scatter(dif_pressure[i],Vp_oil_dif[i,0,320],c= 'r')
    ax.scatter(dif_pressure[i],Vs_oil_dif[i,0,320],c= 'g')
    ax.scatter(dif_pressure[i],rho_oil_dif[i,0,320],c='y')
ax.set_xlabel('Pore pressure change Pa')
ax.set_ylabel('Difference % ')


axis([-5*1e6, 5*1e6, -25.0,5])



fig =  figure(14, figsize=(10, 8))
ax.set_xlabel('Differential pressure Pa')
ax.set_ylabel('Difference % ')


for i in range (NP):
    scatter(P[i],Vp_co2_dif[i,0,320])
    scatter(P[i],Vp_co2_dif[i,1,320],c= 'r')
    scatter(P[i],Vp_co2_dif[i,5,320],c= 'c')
    scatter(P[i],Vp_co2_dif[i,9,320],c='m')
    scatter(P[i],Vs_co2_dif[i,0,320],c= 'b',marker='v',s=50)
    scatter(P[i],Vs_co2_dif[i,1,320],c= 'r',marker='v',s=50)
    scatter(P[i],Vs_co2_dif[i,5,320],c= 'c',marker='v',s=50)
    scatter(P[i],Vs_co2_dif[i,9,320],c='m',marker='v',s=50)



fig =  figure(10, figsize=(6, 5))

ax = plt.subplot(1,1,1)
#ax.set_xlabel('Differential pressure')
ax.set_ylabel('Difference % ')
#ax.set_title('Pore pressure change ')


axis([0, 1, -25.0,5])




for i in range (nCO2):
    ax.scatter(dif_sat[i],VpVs_co2_difs[4,i,320],c='b')
    ax.scatter(dif_sat[i],Vp_co2_difs[4,i,320],c= 'r')
    ax.scatter(dif_sat[i],Vs_co2_difs[4,i,320],c= 'g')
    ax.scatter(dif_sat[i],rho_co2_difs[4,i,320],c='y')
ax.set_xlabel('CO2 saturation')
ax.set_ylabel('Difference % ')
#ax.set_title('Saturation change ')

fig =  figure(11, figsize=(10, 8))

ax = plt.subplot(1,1,1)

'''
for i in range (NP):
    ax.scatter(Ip_co2_dif[i,0,320],VpVs_co2_dif[i,0,320] , c='b')
    ax.scatter(Ip_co2_dif[i,5,320],VpVs_co2_dif[i,5,320], c='b')
    ax.scatter(Ip_co2_dif[i,9,320],VpVs_co2_dif[i,9,320], c='b')
'''

for i in range (nCO2):
    ax.scatter(Ip_oil_dif[0,i,320],VpVs_oil_dif[0,i,320] , c='b')
    ax.scatter(Ip_oil_dif[1,i,320],VpVs_oil_dif[1,i,320] , c='k')
    ax.scatter(Ip_oil_dif[2,i,320],VpVs_oil_dif[2,i,320], c='k')
    ax.scatter(Ip_oil_dif[3,i,320],VpVs_oil_dif[3,i,320], c='k')
    ax.scatter(Ip_oil_dif[4,i,320],VpVs_oil_dif[4,i,320], c='k')
    ax.scatter(Ip_oil_dif[5,i,320],VpVs_oil_dif[5,i,320], c='k')
    ax.scatter(Ip_oil_dif[6,i,320],VpVs_oil_dif[6,i,320], c='k')
    ax.scatter(Ip_oil_dif[7,i,320],VpVs_oil_dif[7,i,320], c='k')
    ax.scatter(Ip_oil_dif[8,i,320],VpVs_oil_dif[8,i,320], c='m')
    ax.scatter(Ip_oil_dif[9,i,320],VpVs_oil_dif[9,i,320], c='k')


'''
for i in range (nCO2):
 #   ax.scatter(Ip_co2_difs[4,i,320],VpVs_co2_difs[4,i,320] ,c= 'c',marker='v',s=50)
#    ax.scatter(Ip_co2_difs[0,i,320],VpVs_co2_difs[0,i,320] ,c= 'c',marker='v',s=50)
   
   
    
    ax.scatter(Ip_co2_difs[1,i,320],Is_co2_difs[1,i,320] , c='r')
    ax.scatter(Ip_co2_difs[2,i,320],Is_co2_difs[2,i,320], c='r')
    ax.scatter(Ip_co2_difs[3,i,320],Is_co2_difs[3,i,320], c='r')
    ax.scatter(Ip_co2_difs[4,i,320],Is_co2_difs[4,i,320], c='r')
    ax.scatter(Ip_co2_difs[5,i,320],Is_co2_difs[5,i,320], c='r')
    ax.scatter(Ip_co2_difs[6,i,320],Is_co2_difs[6,i,320], c='r')
    ax.scatter(Ip_co2_difs[7,i,320],Is_co2_difs[7,i,320], c='r')
    ax.scatter(Ip_co2_difs[8,i,320],Is_co2_difs[8,i,320], c='r')
    ax.scatter(Ip_co2_difs[9,i,320],Is_co2_difs[9,i,320], c='r')
    
'''
'''

for i in range (nCO2):
    ax.scatter(Ip_co2_difps[0,i,320],VpVs_co2_difps[0,i,320] , c='b')
    ax.scatter(Ip_co2_difps[1,i,320],VpVs_co2_difps[1,i,320] , c='b')
    ax.scatter(Ip_co2_difps[2,i,320],VpVs_co2_difps[2,i,320], c='b')
    ax.scatter(Ip_co2_difps[3,i,320],VpVs_co2_difps[3,i,320], c='b')
    ax.scatter(Ip_co2_difps[4,i,320],VpVs_co2_difps[4,i,320], c='b')
    ax.scatter(Ip_co2_difps[5,i,320],VpVs_co2_difps[5,i,320], c='b')
    ax.scatter(Ip_co2_difps[6,i,320],VpVs_co2_difps[6,i,320], c='b')
    ax.scatter(Ip_co2_difps[7,i,320],VpVs_co2_difps[7,i,320], c='b')
    ax.scatter(Ip_co2_difps[8,i,320],VpVs_co2_difps[8,i,320], c='b')
    ax.scatter(Ip_co2_difps[9,i,320],VpVs_co2_difps[9,i,320], c='b')
'''


'''
for i in range (nCO2):
    ax.scatter(Ip_co2_difpsb[0,i,320],VpVs_co2_difpsb[0,i,320] , c='g')
    ax.scatter(Ip_co2_difpsb[1,i,320],VpVs_co2_difpsb[1,i,320] , c='g')
    ax.scatter(Ip_co2_difpsb[2,i,320],VpVs_co2_difpsb[2,i,320], c='g')
    ax.scatter(Ip_co2_difpsb[3,i,320],VpVs_co2_difpsb[3,i,320], c='g')
    ax.scatter(Ip_co2_difpsb[4,i,320],VpVs_co2_difpsb[4,i,320], c='g')
    ax.scatter(Ip_co2_difpsb[5,i,320],VpVs_co2_difpsb[5,i,320], c='g')
    ax.scatter(Ip_co2_difpsb[6,i,320],VpVs_co2_difpsb[6,i,320], c='g')
    ax.scatter(Ip_co2_difpsb[7,i,320],VpVs_co2_difpsb[7,i,320], c='g')
    ax.scatter(Ip_co2_difpsb[8,i,320],VpVs_co2_difpsb[8,i,320], c='g')
    ax.scatter(Ip_co2_difpsb[9,i,320],VpVs_co2_difpsb[9,i,320], c='g')

'''
'''
for i in range (nCO2):
    ax.scatter(Ip_oil_difps[0,i,320],Is_oil_difps[0,i,320] , c='m')
    ax.scatter(Ip_oil_difps[1,i,320],Is_oil_difps[1,i,320] , c='m')
    ax.scatter(Ip_oil_difps[2,i,320],Is_oil_difps[2,i,320], c='m')
    ax.scatter(Ip_oil_difps[3,i,320],Is_oil_difps[3,i,320], c='m')
    ax.scatter(Ip_oil_difps[4,i,320],Is_oil_difps[4,i,320], c='m')
    ax.scatter(Ip_oil_difps[5,i,320],Is_oil_difps[5,i,320], c='m')
    ax.scatter(Ip_oil_difps[6,i,320],Is_oil_difps[6,i,320], c='m')
    ax.scatter(Ip_oil_difps[7,i,320],Is_oil_difps[7,i,320], c='m')
    ax.scatter(Ip_oil_difps[8,i,320],Is_oil_difps[8,i,320], c='m')
    ax.scatter(Ip_oil_difps[9,i,320],Is_oil_difps[9,i,320], c='m')

#ax.scatter(Ip_co2_difst[0,i,320],VpVs_co2_difst[0,i,320], c='c')
#ax.scatter(Ip_co2_difst[9,i,320],VpVs_co2_difst[9,i,320], c='c')

'''
'''





'''

'''
for i in range (nCO2):
    ax.scatter(Ip_co2_difs[4,i,320],VpVs_co2_difs[4,i,320] ,c= 'r',marker='v',s=50)
    ax.scatter(Ip_oil_difs[4,i,320],VpVs_oil_difs[4,i,320] ,c= 'r',marker='v',s=50)
    ax.scatter(Ip_co2_dif[4,i,320],VpVs_co2_dif[4,i,320] ,c= 'k',marker='v',s=50)
   

   
''' 



ax.set_xlabel('Impedance %')
ax.set_ylabel('Vp/Vs ration %')
ax.set_title('Saturation ')

print(CO2)


plt.show()


