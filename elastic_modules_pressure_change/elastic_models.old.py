#!/usr/bin/env python
import numpy as np
import elastic_module as em
import elastic_module as ep
import gassman_module as gm

import matplotlib.pyplot as plt

from pylab import *
import numpy as np
import matplotlib.pyplot as plt

'''
Mineral enpoint
'''
Kqz = 37.0e09   #Pa

Gqz = 44.0e9  #44.0e09   #Pa

Kclay = 25e09 #25.0e09 #Pa

Gclay = 9e09 #9.0e09  #Pa 

Kc = 94.0e09 # 70.2e09   #Pa

Gc = 45.0e09 #29.e09    #Pa

Rhoqz = 2650    #Pa

Rhoclay = 2550    #Pa

clayVol =  0.1 #= Clay proportion Paluxy con mucho clay

Ks  = em.Ksget(Kqz, Kclay, clayVol)

Gs  = em.Ksget(Gqz, Gclay, clayVol)

Rhos  =  2632 # 2632 for Paluxy dirty em.Ksget(Rhoqz, Rhoclay, clayVol)

PRs = em.poisson(Ks,Gs)

PRc = em.poisson(Kc,Gc)

Mc  = em.M_value(Kc,Gc)


'''
Herts-Mindlin endpoint
'''
#phi_c = 0.4001 # Critical porosity, hard coded
phi_c = 0.36010 # Critical porosity, hard coded



rhoco2  = [0.67*1e3, 0.66*1e3, 0.64*1e3,0.6*1e3,0.56*1e3, 0.52*1e3, 0.48*1e3,0.44*1e3,0.33*1e3, 0.35*1e3  ] #0.0636*1e9
rhobrine= [0.998*1e3,0.99*1e3,0.99*1e3,0.99*1e3,0.99*1e3, 0.998*1e3,0.99*1e3,0.99*1e3,0.99*1e3,0.99*1e3] #2.5016*1e9
rhooil  = [0.796*1e3,0.79*1e3,0.79*1e3,0.79*1e3,0.77*1e3,0.75*1e3, 0.73*1e3,0.70*1e3,0.63*1e3,0.64*1e3] 
Kco2    = [0.125*1e9,0.1125*1e9,0.107*1e9,0.08*1e9,0.07*1e9,0.05*1e9,0.04*1e9,0.03*1e9,0.02*1e9,0.019*1e9] #0.3721*1e3
Kbrine  = [2.53*1e9,2.50*1e9,2.52*1e9,2.51*1e9,2.5*1e9,2.49*1e9,2.49*1e9,2.48*1e9,2.47*1e9,2.47*1e9] #0.9972*1e3
Koil    = [1.13*1e9,1.1*1e9,1.09*1e9,1.07*1e9,1.07*1e9,1.05*1e9,1.04*1e9,1.03*1e9,1.02*1e9,1.00*1e9]

'''
 
rhoco2  = [0.3511*1e3, 0.3388*1e3, 0.4447*1e3,0.4875*1e3,0.5285*1e3 ] #0.0636*1e9
rhobrine= [0.9954*1e3,0.9953*1e3,0.9962*1e3,0.9967*1e3,0.9971*1e3] #2.5016*1e9
rhooil  =[0.6422*1e3,0.6328*1e3,0.702*1e3,0.7339*1e3,0.753*1e3] 
Kco2    = [0.01996*1e9,0.0214*1e9,0.0312*1e9,0.0442*1e9,0.0576*1e9] #0.3721*1e3
Kbrine  = [2.47*1e9, 2.48*1e9, 2.49*1e9, 2.4911*1e9, 2.49*1e9] #0.9972*1e3
Koil    =[0.0365*1e9,0.0338*1e9,0.0716*1e9,0.1065*1e9,0.1725*1e9]
''' 





n = 20.0 - 34.0 *phi_c +14.0*np.power(phi_c,2) #relation between cordination number and porosity
n *=1.

Po = 12.5e6
P = np.arange( 3.5*1e6, 13.5*1e6, 1.0*1e6 ) #Effective pressure
NP=P.size
Pstiff = np.arange( 600.0*1e6,601.0*1e6, 1.0*1e6   )

##################Pressure changes##########################
'''Hertz mindlin'''


#phi= np.arange(0.0001,0.40010000050,0.001)

phi= np.arange(0.0001,0.3601000050,0.001)

Ft=0

Np=phi.size
Khm = em.K_hertz_mindlin(n, phi_c,Gs, PRs, P)
(Ghm,Ghm_slip) = em.G_hertz_mindlin(n, phi_c,Gs, PRs, P,Ft)
Zhm = em.Z_bound(Khm,Ghm,P)
Zhm_slip = em.Z_bound(Khm,Ghm_slip,P)

Khm_stiff = em.K_hertz_mindlin(n, phi_c,Gs, PRs, Pstiff)
Ghm_stiff = em.G_hertz_mindlin(n, phi_c,Gs, PRs, Pstiff,Ft)
Zhm_stiff = em.Z_bound(Khm,Ghm,P)



#############################################################

'''
Modified lower Hashin-Shtrikman bound
'''

Kdry_l = em.K_sashin_shtrikman_lower(Ks,Khm,Ghm,phi_c,phi,P)
Gdry_l = em.G_sashin_shtrikman_lower(Ghm,Gs,phi_c,phi,Zhm,Ks,P)


'''
Modified upper Hashin-Shtrikman bound
'''
Z_u = em.Z_value(Ks,Gs,P)
Kdry_u = em.K_sashin_shtrikman_upper2(Ks,Khm,Gs,phi_c,phi,P)
Gdry_u = em.G_sashin_shtrikman_upper2(Ghm,Gs,phi_c,phi,Z_u,Ks,P)



'''
Upper bound contact cement
'''


Kdry_cement = em.Kdry_cement(Kc,Gc,Ks,PRc,PRs,phi_c,Mc,phi,Gs,P)

Gdry_cement = em.Gdry_cement(phi,phi_c,Gc,Gs,PRs,Kdry_cement,Ks)



################K Sat##########################################


Np=phi.size
CO2 =np.arange(0.0,1.0,0.1)
nCO2=CO2.size
Ksat_l1 = np.zeros(((NP,Np,nCO2)))
rhosat_l1 = np.zeros(((NP,Np,nCO2)))

colors=['k','b','r','m','y','g','c','b']
for i in range(NP):
    (Ksat_l,rhosat_l)= gm.Gassman(Kdry_l[i,:],Ks,phi, Kbrine[i],Rhos,Kco2[i],rhobrine[i],rhoco2[i])
    for j in range(Np):     
        for k in range(nCO2): 
            Ksat_l1[i][j][k]=Ksat_l[j][k]
            rhosat_l1[i][j][k]=rhosat_l[j][k]


(Ksat_c,rhosat)= gm.Gassman(Kdry_cement[:],Ks,phi,Kbrine[0],Rhos,Kco2[0],rhobrine[0],rhoco2[0])




############################################################


fig = figure(3, figsize=(10, 8))
Nconstant=50

Ksat_constant = np.zeros(((NP,Nconstant,Np)))
Gsat_constant = np.zeros(((NP,Nconstant,Np)))

Kconstant_perc = np.zeros(((NP,Nconstant)))



Kdry_constant = np.zeros(((NP,Np)))

KW = np.zeros(((Np)))
GW = np.zeros(((Np)))


Gdry_constant = np.zeros(((NP,Np)))
nmodel=9  #   3 = model que setee para Paluxy, con cementacion y mas contenido de clay
Gnmodel=7 
Kdry_stiff = np.zeros(((NP,Np)))


Gdry_stiff = np.zeros(((NP,Np)))


for i in range(NP):
    (Kdry_constant1,Gdry_constant1,Kconstant_perc2,Gconstant_perc2)=em.K_constant_cement(Ks,Kdry_cement[:],Gdry_cement[:],Gs,phi,Khm,Ghm)

    for k in range(Np): 
        Kdry_constant[i][k]=Kdry_constant1[nmodel][k]
        Gdry_constant[i][k]=Gdry_constant1[Gnmodel][k]



print(Kconstant_perc2)  
print(Gconstant_perc2)  
print(Ghm)  
print(Ghm_slip) 
nmodel=22
Gnmodel=17 

for i in range(NP):
    (Kdry_constant1,Gdry_constant1,Kconstant_perc2,Gconstant_perc2)=em.K_constant_cement(Ks,Kdry_cement[:],Gdry_cement[:],Gs,phi,Khm,Ghm)

    for k in range(Np): 
        Kdry_stiff[i][k]=Kdry_constant1[nmodel][k]
        Gdry_stiff[i][k]=Gdry_constant1[Gnmodel][k]


print(P)

for i in range(NP):
    if P[i]==Po:
       Kdry_constant=Kdry_constant
       Gdry_constant=Gdry_constant
    #   GW[k]=(Gdry_constant[4][k]- Gdry_l[4][k])/(Gdry_stiff[4][k]- Gdry_l[4][k])   
    #   GW[k]=(1.0-KW[k])*Gdry_l[4][k]/Gs + KW[k]*Gdry_stiff[4][k]/Gs
    #   GW[k]=(1.0-KW[k])*Gdry_l[4][k]/Gdry_constant[4][k] + KW[k]*Gdry_stiff[4][k]/Gdry_constant[4][k]


   #    Gdry_constant[i][k]=(1.0-GW[k])*Gdry_l[i][k]+GW[k]*Gdry_stiff[i][k]


    else:
       for k in range(Np): 
           KW[k]=(Kdry_constant[9][k]- Kdry_l[9][k])/(Kdry_stiff[9][k]- Kdry_l[9][k])    
       #    GW[k]=(1.0-KW[k])*Gdry_l[4][k]/Gs + KW[k]*Gdry_stiff[4][k]/Gs
        #   GW[k]=(1.0-KW[k])*Gdry_l[4][k]/Gdry_constant[4][k] + KW[k]*Gdry_stiff[4][k]/Gdry_constant[4][k]

           GW[k]=(Gdry_constant[9][k]- Gdry_l[9][k])/(Gdry_stiff[9][k]- Gdry_l[9][k])   
           Kdry_constant[i][k]=(1.0-KW[k])*Kdry_l[i][k]+KW[k]*Kdry_stiff[i][k]
           Gdry_constant[i][k]=(1.0-GW[k])*Gdry_l[i][k]+GW[k]*Gdry_stiff[i][k]



for i in range(NP):
    (Ksat_constant1,Gsat_constant1,Kconstant_perc1,Gconstant_perc1)=em.K_constant_cement(Ks,Ksat_c[:,0],Gdry_cement,Gs,phi,Khm,Ghm)
    
    for j in range(Nconstant):     
        Kconstant_perc[i][j]=Kconstant_perc1[j]
        for k in range(Np): 
            Ksat_constant[i][j][k]=Ksat_constant1[j][k]
            Gsat_constant[i][j][k]=Gsat_constant1[j][k]
  
        
plot(phi[:],Ksat_constant[4,9,:],'b')
plot(phi[:],Ksat_constant[4,18,:],'b')

plot(phi[:],Gdry_cement[:],'k')
#plot(phi[:],Gdry_cement[:],'c')

plot(phi[:],Gdry_constant[4,:],'b')

#plot(phi[:],Gdry_constant[4,:],'m')

#plot(phi[:],Gdry_constant[0,:],'m')

#plot(phi[:],Kdry_constant[0,:],'c')

#plot(phi[:],Gdry_l[0,:],'c')
#plot(phi[:],Gdry_l[4,:],'c')
plot(phi[:],Gdry_l[4,:],'r')
plot(phi[:],Gdry_stiff[4,:],'c')


plot(phi[:],Ksat_l1[4,:,0],'r')
plot(phi[:],Ksat_l1[0,:,0],'r')
plot(phi[:],Ksat_c[:,0],'r')



phi_d,phi_t,Ksat_well,Vsh_iw,G,depth = np.loadtxt('./xplot1_tmp.txt',unpack=True)

ax = plt.subplot(1,1,1)
ax.set_xlabel('porosity $\phi$')
ax.set_ylabel('Bulk modulus $K$')
p=ax.scatter(phi_t/100,Ksat_well,c=Vsh_iw)
cbar= plt.colorbar(p)
cbar.set_label('Vsh')
axis([0, 0.5, 0, 5*1e10])

Ndepth=depth.size

for i in range(Ndepth):
    if depth[i] <=968.0:
	#p=bx.scatter(phi_t[i],Ksat[i],color='r')
       k=10
    elif depth[i] >=968.0 and depth[i]<=972.0:
       p=ax.scatter(phi_d[i]/100,Ksat_well[i],color='b')
    elif depth[i] >=973.0 and depth[i]<=978.0:
       p=ax.scatter(phi_d[i]/100,Ksat_well[i],color='y')
    elif depth[i] >=979.0 and depth[i]<=984.0:
       p=ax.scatter(phi_d[i]/100,Ksat_well[i],color='m')
   # elif depth[i] >=989.0 and depth[i]<=1001:
    elif depth[i] >=3388.0 and depth[i]<=3446:
       k=10
       p=ax.scatter(phi_t[i],Ksat_well[i],color='k')
    else:
       k=10

p0 = Rectangle((0, 0), 0.5, 0.5, fc="b")
p1 = Rectangle((0, 0), 0.5, 0.5, fc="g")
p2 = Rectangle((0, 0), 0.5, 0.5, fc="r")
p3 = Rectangle((0, 0), 0.5, 0.5, fc="c")
p4 = Rectangle((0, 0), 0.5, 0.5, fc="m")
p5 = Rectangle((0, 0), 0.5, 0.5, fc="y")
p6 = Rectangle((0, 0), 0.5, 0.5, fc="k")

labels = ('Paluxy' , 'Tusc3', 'Tusc5', \
          'Tusc7')
legend = plt.legend([p6,p4,p5,p0],labels, loc=(0.5, .7), labelspacing=0.1)
ltext = gca().get_legend().get_texts()






################K Sat once we find the best model#############################

Ksat_l_contact   = np.zeros(((NP,Np,nCO2)))
rhosat_l_contact = np.zeros(((NP,Np,nCO2)))
Ksat_l_contact_oil   = np.zeros(((NP,Np,nCO2)))
rhosat_l_contact_oil = np.zeros(((NP,Np,nCO2)))



for i in range(NP):
    (Ksat_l2,rhosat_l2)= gm.Gassman(Kdry_constant[i,:],Ks,phi, Kbrine[i],Rhos,Kco2[i],rhobrine[i],rhoco2[i])
    (Ksat_l3,rhosat_l3)= gm.Gassman(Kdry_constant[i,:],Ks,phi, Kbrine[i],Rhos,Koil[i],rhobrine[i],rhooil[i])
#    (Ksat_l2,rhosat_l2)= gm.Gassman(Kdry_l[i,:],Ks,phi, Kbrine[i],Rhos,Kco2[i],rhobrine[i],rhoco2[i])
#    (Ksat_l3,rhosat_l3)= gm.Gassman(Kdry_l[i,:],Ks,phi, Kbrine[i],Rhos,Koil[i],rhobrine[i],rhooil[i])


    for j in range(Np):     
        for k in range(nCO2): 
            Ksat_l_contact[i][j][k]=Ksat_l2[j][k]
            rhosat_l_contact[i][j][k]=rhosat_l2[j][k]
            Ksat_l_contact_oil[i][j][k]=Ksat_l3[j][k]
            rhosat_l_contact_oil[i][j][k]=rhosat_l3[j][k]


print(Ksat_l_contact[0][0][:])

print(Ksat_l_contact_oil[0][0][:])
print(rhosat_l_contact[0][0][:])

print(rhosat_l_contact_oil[0][0][:])


fig = figure(3, figsize=(10, 8))


######################VS#####################################




nvol=phi.size
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

       
        Vssat_l1= em.vs(Gdry_constant[i,:],rhosat_l_contact[i,:,k],phi)
        Vpsat_l1= em.vp(Ksat_l_contact[i,:,k],Gdry_constant[i,:],rhosat_l_contact[i,:,k],phi)
        Vp_Vs_sat_l1= em.vp_vs(Vpsat_l1,Vssat_l1,phi)
        Ip_sat_l1=em.Ip(Vpsat_l1,rhosat_l_contact[i,:,k],phi) 
        Is_sat_l1=em.Ip(Vssat_l1,rhosat_l_contact[i,:,k],phi) 


        Vssat_l3= em.vs(Gdry_constant[i,:],rhosat_l_contact_oil[i,:,k],phi)
        Vpsat_l3= em.vp(Ksat_l_contact_oil[i,:,k],Gdry_constant[i,:],rhosat_l_contact_oil[i,:,k],phi)
        Vp_Vs_sat_l3= em.vp_vs(Vpsat_l3,Vssat_l3,phi)
        Ip_sat_l3=em.Ip(Vpsat_l3,rhosat_l_contact_oil[i,:,k],phi) 
        Is_sat_l3=em.Ip(Vssat_l3,rhosat_l_contact_oil[i,:,k],phi) 


       
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


 


fig = figure(5, figsize=(10, 8))

ax = plt.subplot(1,1,1)

PR,Ip,Vsh_iw,phi_d,phi_t,SO,vp,vs,mu,tt = np.loadtxt('./xplot21.txt',unpack=True)
#PR,Ip,Vsh_iw,phi_d,SO,vp,vs,tt = np.loadtxt('./xplot2_159.tmp.txt',unpack=True)

#PR,Ip,Vsh_iw,phi_d,SO,vp,vs,tt= np.loadtxt('./xplot2_cat.txt',unpack=True)

p=ax.scatter(Ip_sat_l_oil[9,0,:],Vp_Vs_sat_l_oil[9,0,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,1,:],Vp_Vs_sat_l_oil[9,1,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,9,:],Vp_Vs_sat_l_oil[9,9,:],color='g')
p=ax.scatter(Ip_sat_l_oil[0,1,:],Vp_Vs_sat_l_oil[0,1,:],color='g')




p=ax.plot(Ip_sat_l[9,0,:],Vp_Vs_sat_l[9,0,:],color='r')
p=ax.plot(Ip_sat_l[9,1,:],Vp_Vs_sat_l[9,1,:],color='r')
p=ax.plot(Ip_sat_l[9,9,:],Vp_Vs_sat_l[9,9,:],color='r')
p=ax.plot(Ip_sat_l[0,1,:],Vp_Vs_sat_l[0,9,:],color='r')



p=ax.scatter(Ip,PR,c=Vsh_iw)
cbar= plt.colorbar(p)
'''
p=ax.scatter(Ip_sat_l_oil[4,0,:],Vp_Vs_sat_l_oil[4,0,:],color='g')
p=ax.scatter(Ip_sat_l_oil[4,9,:],Vp_Vs_sat_l_oil[4,9,:],color='g')
p=ax.scatter(Ip_sat_l_oil[4,1,:],Vp_Vs_sat_l_oil[4,1,:],color='g')


p=ax.scatter(Ip_sat_l[4,0,:],Vp_Vs_sat_l[4,0,:],color='k')
p=ax.scatter(Ip_sat_l[4,9,:],Vp_Vs_sat_l[4,9,:],color='k')
p=ax.scatter(Ip_sat_l[4,1,:],Vp_Vs_sat_l[4,1,:],color='k')
p=ax.scatter(Ip_sat_l[4,9,:],Vp_Vs_sat_l[4,9,:],color='k')
'''
axis([4*1e6, 12*1e6, 1.0,3])
ax.set_ylabel('Vp/Vs ratio ')
ax.set_xlabel('Acoustic impedance kg/(m2*s)')

fig = figure(6, figsize=(10, 8))

ax = plt.subplot(1,1,1)

axis([4*1e6, 12*1e6, 1.0,3])

p=ax.scatter(Ip_sat_l_oil[9,0,:],Vp_Vs_sat_l_oil[9,0,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,1,:],Vp_Vs_sat_l_oil[9,1,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,9,:],Vp_Vs_sat_l_oil[9,9,:],color='g')


p=ax.scatter(Ip_sat_l[9,0,:],Vp_Vs_sat_l[9,0,:],color='r')
p=ax.scatter(Ip_sat_l[9,1,:],Vp_Vs_sat_l[9,1,:],color='r')
p=ax.scatter(Ip_sat_l[9,9,:],Vp_Vs_sat_l[9,9,:],color='r')


#p=ax.plot(Ip_sat_l_oil[4,0,:],Vp_Vs_sat_l_oil[4,0,:],color='k')
#p=ax.plot(Ip_sat_l_oil[0,0,:],Vp_Vs_sat_l_oil[0,0,:],color='k')

#p=ax.plot(Ip_sat_l_oil[4,9,:],Vp_Vs_sat_l_oil[4,9,:],color='c')
#p=ax.plot(Ip_sat_l_oil[0,9,:],Vp_Vs_sat_l_oil[0,9,:],color='r')


#p=ax.plot(Ip_sat_l_oil[4,1,:],Vp_Vs_sat_l_oil[4,1,:],color='c')
#p=ax.plot(Ip_sat_l_oil[0,1,:],Vp_Vs_sat_l_oil[0,1,:],color='r')


#p=ax.plot(Ip_sat_l_oil[3,1,:],Vp_Vs_sat_l_oil[3,1,:],color='g')
#p=ax.plot(Ip_sat_l_oil[3,9,:],Vp_Vs_sat_l_oil[3,9,:],color='g')

#p=ax.plot(Ip_sat_l[3,0,:],Vp_Vs_sat_l[3,0,:],color='k')
#p=ax.plot(Ip_sat_l[3,1,:],Vp_Vs_sat_l[3,1,:],color='r')
#p=ax.plot(Ip_sat_l[3,9,:],Vp_Vs_sat_l[3,9,:],color='r')



#p=ax.plot(Ip_sat_l_oil[0,0,:],Vp_Vs_sat_l_oil[0,0,:],color='k')
#p=ax.plot(Ip_sat_l_oil[0,2,:],Vp_Vs_sat_l_oil[0,2,:],color='k')
#p=ax.plot(Ip_sat_l_oil[0,9,:],Vp_Vs_sat_l_oil[0,9,:],color='k')

'''
Ip_paluxy=np.zeros(13)
PR_paluxy=np.zeros(13)
SO_paluxy=np.zeros(13)
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
    elif depth[i] >=989.0 and depth[i]<=1001:
  #     p=ax.scatter(Ip[i],PR[i],c=SO[i])
       
       Ip_paluxy[j]= Ip[i]
       PR_paluxy[j]= PR[i]
       SO_paluxy[j]= SO[i]
       j=j+1

    else:
       k=10
'''
p=ax.scatter(Ip,PR,c=SO)
cbar= plt.colorbar(p)
cbar.set_label('Oil')
axis([4*1e6, 12*1e6, 1.0, 3])

ax.set_xlabel('Impedance')
ax.set_ylabel('Vp/Vs ratio')



fig = figure(7, figsize=(10, 8))

ax = plt.subplot(1,1,1)

p=ax.scatter(Ip_sat_l_oil[9,0,:],Vp_Vs_sat_l_oil[9,0,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,9,:],Vp_Vs_sat_l_oil[9,9,:],color='g')
p=ax.scatter(Ip_sat_l_oil[9,1,:],Vp_Vs_sat_l_oil[9,1,:],color='g')


p=ax.scatter(Ip_sat_l[9,0,:],Vp_Vs_sat_l[9,0,:,],color='k')
p=ax.scatter(Ip_sat_l[9,9,:],Vp_Vs_sat_l[9,9,:,],color='k')
p=ax.scatter(Ip_sat_l[9,1,:],Vp_Vs_sat_l[9,1,:,],color='k')
p=ax.scatter(Ip_sat_l[9,9,:],Vp_Vs_sat_l[9,9,:,],color='k')


for i in range(Np):
    if phi[i]==.1051:
       p=ax.plot(Ip_sat_l[4,0,i],Vp_Vs_sat_l[4,0,i],'-o', alpha=0.7,ms=20,lw=2,mfc='b')

    elif  phi[i]==.2001:
       p=ax.plot(Ip_sat_l[4,0,i],Vp_Vs_sat_l[4,0,i],'-o', alpha=0.7,ms=20,lw=2,mfc='g')

    elif  phi[i]==.3051:
       p=ax.plot(Ip_sat_l[4,0,i],Vp_Vs_sat_l[4,0,i],'-o', alpha=0.7,ms=20,lw=2,mfc='r')

    else:
       k=1



p=ax.scatter(Ip,PR,c=phi_d)
cbar= plt.colorbar(p)
cbar.set_label('Porosity')
axis([4*1e6, 18*1e6, 1, 3])

ax.set_xlabel('Impedance')
ax.set_ylabel('Vp/Vs ratio')


fig = figure(8, figsize=(10, 8))

ax = plt.subplot(1,1,1)

plot(phi[:],Vpsat_l_oil[0,0,:],color='g')
plot(phi[:],Vpsat_l_oil[2,0,:],color='b')
plot(phi[:],Vpsat_l_oil[9,0,:],color='r')

scatter(phi_d,vp,c=Vsh_iw)

cbar= plt.colorbar(p)
cbar.set_label('Vsh')
ax.set_xlabel('Porosity')
ax.set_ylabel('Vp (m/s)')




fig = figure(13, figsize=(10, 8))

ax = plt.subplot(1,1,1)

plot(phi[:],Vpsat_l_oil[0,0,:],color='g')
plot(phi[:],Vpsat_l_oil[2,0,:],color='b')
plot(phi[:],Vpsat_l_oil[9,0,:],color='r')

scatter(phi_d,vp,c=SO)

cbar= plt.colorbar(p)
cbar.set_label('SO')


fig = figure(9, figsize=(10, 8))

ax = plt.subplot(1,1,1)

plot(phi[:],Vssat_l_oil[0,0,:],color='g')
plot(phi[:],Vssat_l_oil[2,0,:],color='b')
plot(phi[:],Vssat_l_oil[9,0,:],color='r')

plot(phi[:],Vssat_l_oil[0,1,:],color='g')
plot(phi[:],Vssat_l_oil[2,1,:],color='b')
plot(phi[:],Vssat_l_oil[9,1,:],color='r')
ax.set_xlabel('Porosity')
ax.set_ylabel('Vs (m/s)')

scatter(phi_d,vs,c=Vsh_iw)
cbar= plt.colorbar(p)
cbar.set_label('Vsh')
ax.set_xlabel('Porosity')
ax.set_ylabel('Vs (m/s)')




fig = figure(12, figsize=(10, 8))

plot(phi[:],Vssat_l_oil[0,0,:],color='g')
plot(phi[:],Vssat_l_oil[2,0,:],color='b')
plot(phi[:],Vssat_l_oil[4,0,:],color='r')

plot(phi[:],Vssat_l_oil[0,1,:],color='g')
plot(phi[:],Vssat_l_oil[2,1,:],color='b')
plot(phi[:],Vssat_l_oil[4,1,:],color='r')

scatter(phi_d,vs,c=SO)
cbar= plt.colorbar(p)
cbar.set_label('so')


fig = figure(10, figsize=(10, 8))

plot(phi[:],Vpsat_l[9,1,:],color='c')
plot(phi[:],Vpsat_l[0,1,:],color='r')
plot(phi[:],Vpsat_l[9,9,:],color='c')
plot(phi[:],Vpsat_l[0,9,:],color='r')

plot(phi[:],Vssat_l[9,1,:],color='c')
plot(phi[:],Vssat_l[0,1,:],color='r')

plot(phi[:],Vssat_l[9,9,:],color='c')
plot(phi[:],Vssat_l[0,9,:],color='r')


fig = figure(11, figsize=(10, 8))

ax = plt.subplot(1,1,1)
ax.set_title('Pore pressure diff 500 PSI')

ax.set_xlabel('Norm diff Impedance')
ax.set_ylabel('Norm diff Vp/Vs')

labels = ('100% brine' , '90% CO2', '50% CO2', \
          '10% CO2')
legend = plt.legend(labels, loc=(0.5, .7), labelspacing=0.1)
ltext = gca().get_legend().get_texts()





fig = figure(14, figsize=(10, 8))

p=scatter(phi,KW[:],c=KW[:])
cbar= plt.colorbar(p)
cbar.set_label('W')

#plot(phi,GW[:],'m')


#axis([0, 5, 0, 1])


fig = figure(15, figsize=(10, 8))

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

#print(Vpsat_l_oil[0,0,0])

(Vp_co2_dif,dif_pressure)=em.Differences(Vpsat_l[:,:,:],P[:],phi[:])
(Vs_co2_dif,dif_pressure)=em.Differences(Vssat_l[:,:,:],P[:],phi[:])
(Vs_co2_dif,dif_pressure)=em.Differences(Vssat_l[:,:,:],P[:],phi[:])
(VpVs_co2_dif,dif_pressure)=em.Differences(Vp_Vs_sat_l[:,:,:],P[:],phi[:])







(Vp_oil_dif,dif_pressure)=em.Differences(Vpsat_l_oil[:,:,:],P[:],phi[:])

(Vs_oil_dif,dif_pressure)=em.Differences(Vssat_l_oil[:,:,:],P[:],phi[:])
(Vs_oil_dif,dif_pressure)=em.Differences(Vssat_l_oil[:,:,:],P[:],phi[:])
(VpVs_oil_dif,dif_pressure)=em.Differences(Vp_Vs_sat_l_oil[:,:,:],P[:],phi[:])
(rho_oil_dif,dif_pressure)=em.Differences(rho_sat_l_oil[:,:,:],P[:],phi[:])
(rho_co2_dif,dif_pressure)=em.Differences(rho_sat_l[:,:,:],P[:],phi[:])






#rho_co2_dif=em.Differences(rhosat_l[:,:,:],P[:],phi[:])


fig =  figure(16, figsize=(10, 8))

ax = plt.subplot(1,1,1)
ax.set_xlabel('Differential pressure')
ax.set_ylabel('Difference % ')



for i in range (NP-1):
    ax.scatter(dif_pressure[i],VpVs_oil_dif[i,9,320],c='b')
    ax.scatter(dif_pressure[i],Vp_oil_dif[i,9,320],c= 'r')
    ax.scatter(dif_pressure[i],Vs_oil_dif[i,9,320],c= 'g')
    ax.scatter(dif_pressure[i],rho_oil_dif[i,9,320],c='y')
ax.set_xlabel('Pore pressure change Pa')
ax.set_ylabel('Difference % ')


axis([0*1e6, 10*1e6, -15.0,10])



fig =  figure(17, figsize=(10, 8))
ax.set_xlabel('Differential pressure Pa')
ax.set_ylabel('Difference % ')


for i in range (NP-1):
    scatter(P[i],Vp_co2_dif[i,0,320])
    scatter(P[i],Vp_co2_dif[i,1,320],c= 'r')
    scatter(P[i],Vp_co2_dif[i,5,320],c= 'c')
    scatter(P[i],Vp_co2_dif[i,9,320],c='m')
    scatter(P[i],Vs_co2_dif[i,0,320],c= 'b',marker='v',s=50)
    scatter(P[i],Vs_co2_dif[i,1,320],c= 'r',marker='v',s=50)
    scatter(P[i],Vs_co2_dif[i,5,320],c= 'c',marker='v',s=50)
    scatter(P[i],Vs_co2_dif[i,9,320],c='m',marker='v',s=50)






fig =  figure(18, figsize=(10, 8))
CO2 =np.arange(0.0,1.0,0.1)
nCO2=CO2.size

#for i in range (nCO2):

scatter(phi,Vpsat_l[8,0,:],c='b')

scatter(phi,Vpsat_l_oil[8,0,:],c='m')
scatter(phi,Vssat_l[8,0,:],c='b')

scatter(phi,Vssat_l_oil[8,0,:],c='m')





#scatter(P,Vpsat_l[:,3,320],c='y')
#scatter(P,Vpsat_l[:,5,320],c='r')
#scatter(P,Vpsat_l[:,7,320],c='g')

#scatter(P,Vpsat_l[:,9,320],c='m')



#vs(Gdry_constant[i,:],rhosat_l_contact[i,:,k],phi)
 #       Vpsat_l1= em.vp(Ksat_l_contact[i,:,k]
plt.show()


