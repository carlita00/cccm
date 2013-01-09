#!/usr/bin/env python

import ElasticMedia as ElasticMedia
import GranularMedia as GranularMedia
import Gassman as Gassman
import GranularMedia as GranularMedia
import Parameters as Parameters

import Plot as Plot
import Init_logs as Logs
import elastic_parameters as ep
import constant_cement as ctc
import matplotlib.pylab as pylab
import matplotlib.mlab as ml

from pylab import *
import numpy as np
import matplotlib.pyplot as plt


#-----Global variables--------
  
k = [37.0e09,25.0e09,123.0e09]    #Pa  
g = [44.e09,9.4e09,51.0e09]       #Pa
f = [.9,.09,.01] # Mineral fraccion
#f = [.65,.15,.2] # Mineral fraccion
rho = [2650.0, 2550.0,3960.0]
phiC = 0.4001 # critical porosity
n = 20.0 - 34.0 *phiC +14.0*np.power(phiC,2) 
c = 1.3
n = n*c
print(n)

#-------Elastic Media------------
  
ks  = ElasticMedia.hillAverage(k,f)
gs  = ElasticMedia.hillAverage(g,f)
prs = ep.poisson_value(ks,gs)
kc = k[1]
gc = g[1]
prc = ep.poisson_value(kc,gc)
mc  = ep.M_value(kc,gc)
print(prc,"poisson",prs)

#--------Elastic Modeling---------

p = np.arange(3.5e06,17.5e06,1.0e06) 
phi = np.arange(0.0001,0.4001,0.001)

'''Hertz mindlin''' 
  
gm = GranularMedia.GranularMedia(n,phiC,gs,ks,prs,p)
khm = gm.applyKhm() 
ghm = gm.applyGhm() 
zhm = gm.applyZhm(khm,ghm) 
z = gm.applyZ() 

'''Modified lower Hashin-Shtrikman bound'''

em = ElasticMedia.ElasticMedia(ks,gs,phiC,phi,p)
kdry_l = em.kHS_lower(khm,ghm) 
gdry_l = em.gHS_lower(ghm,zhm) 

'''Modified upper Hashin-Shtrikman bound'''
kdry_u = em.kHS_upper(khm)
gdry_u = em.gHS_upper(ghm,z)

'''Upper bound contact cement'''
kdry_c = gm.kdryCement(c,phi,kc,gc,prc,mc)
#print(gc)
gdry_c = gm.gdryCement(c,phi,gc,kdry_c)


#---------K Sat-- Gassmman -------
rhos  =  2632 #make sure of this value!!!!!! and check the pressure !!!
pref  = 9 # CHECKKK
gassman = Gassman.Gassman(ks,rhos,phi)
co2 = np.arange(0.0,1,0.1)
n1 = len(phi)
n3 = len(p)
n2 = len(co2)
rhobri,rhooil,rhoco2,kbri,koil,kco2 = np.zeros(n3),np.zeros(n3),np.zeros(n3),\
                                      np.zeros(n3),np.zeros(n3),np.zeros(n3)
for i in range(n3):
  rhoco2[i] = .082e3 #0.6*1e3
  rhobri[i] = 0.99*1e3
  rhooil[i] = 0.666e3 #0.79*1e3
  kco2[i] = 0.018e9 #0.125*1e9
  kbri[i] = 2.47*1e9
  koil[i] = .5e9 #1.*1e9

(ksat_l,rhosat_l)= gassman.fluidSub2_p(kbri,kco2,rhobri,rhoco2,kdry_l,co2,p)
(ksat_c,rhosat_c)= gassman.fluidSub2(kbri[pref],kco2[pref],rhobri[pref],rhooil[pref],kdry_c,co2)

(ksat_loil,rhosat_loil)= gassman.fluidSub2_p(kbri,koil,rhobri,rhooil,kdry_l,co2,p)
(ksat_c,rhosat_c)= gassman.fluidSub2(kbri[pref],kco2[pref],rhobri[pref],rhooil[pref],kdry_c,co2)



#-------Selecting the best trend----
nconstant = 50 # line number of constant cementation
#print(khm)
(kdryConstant,gdryConstant,k_per,g_per) = gm.ConstantCement(kdry_c,gdry_c,phi, \
                                          khm[pref],ghm[pref],nconstant)
print(p)
for i in range(nconstant):
  print(i,k_per[i],g_per[i])
#------------Plot-----------------
dmin = 3207
dmax = 3440
dmin2 = 3186
dmax2 = 3301

fileName1 = './Logs/140.txt'
fileName2 = './Logs/159.txt'

#--- Reading the logs files-------
slowp1,slows1,dtc1,dts1,rhoD1,so1,phit1,facies1,depth1,gr1,bvi1,cbw1,vsh1 = Logs.logread140(fileName1)
slowp2,slows2,dtc2,dts2,rhoD2,so2,phit2,facies2,depth2,gr2,bvi2,cbw2,somril = Logs.logread159(fileName2)


#--- Init a class each log will have one ------
log1 = Logs.Logs(fileName1,dmin,dmax,rho)
log2 = Logs.Logs(fileName2,dmin2,dmax2,rho)

#---- Calculating parameters--------
kSatC1,kSatD1,muDryC1,muDryD1,muSatC1,muSatD1,phit1,so1,vsh1,facies1,depth1,gr1,vpvs1,ip1 = \
log1.parameters(slowp1,slows1,dtc1,dts1,rhoD1,so1,phit1,facies1,depth1,gr1,bvi1,cbw1,vsh1)

kSatC2,kSatD2,muDryC2,muDryD2,muSatC2,muSatD2,phit2,so2,vsh2,facies2,depth2,gr2,vpvs2,ip2 = \
log2.parameters(slowp2,slows2,dtc2,dts2,rhoD2,so2,phit2,facies2,depth2,gr2,bvi2,cbw2,gr2)


#print len(a)
#---- Plotting the logs ------------
log1.plot("140-1")
#log2.plot("159-2")

#------Merge the logs----------------


depth,ksat,mudry,phit,so,gr,facies,vsh,vpvs,ip = Logs.mergelogs(kSatD2,muDryD2,phit2, \
                                                       so2,gr2,facies2,vsh2,vpvs2,ip2,depth2,kSatD1, \
                                                       muDryD1,phit1,so1,gr1,facies1,vsh1,vpvs1,ip1,depth1)
#print(depth)  #n=len(kSatD2)+len(kSatD)
print(len(depth1),len(depth2),len(depth)) 


plots = Plot.Plot(kdry_l,gdry_l,kdry_c,gdry_c,kdryConstant,gdryConstant,phi,ks)
#plots.plotDry(kbri[pref],koil[pref],pref,kSatD2,muDryD2,phit2,so2,gr2,facies2,vsh2)
plots.plotDry(kbri[pref],koil[pref],pref,ksat,mudry,phit,so,gr,facies,vsh)

#kdry_well = plots.plotDryratio(kbri[pref],koil[pref],kSatD2,muDryD2,phit2,so2,gr2,facies2,depth2,dmin2,dmax2,vsh2)
#plots.plotfacies(kSatD2,muDryD2,phit2,vsh2,facies2,depth2,kdry_well,pref)
#kdry_well = plots.plotDryratio(kbri[pref],koil[pref],ksat,mudry,phit,so,gr,facies,depth,dmin2,dmax2,vsh)
#plots.plotfacies(kSatD2,muDryD2,phit2,vsh2,facies2,depth2,kdry_well,pref)
#plots.plotfacies(ksat,mudry,phit,vsh,facies,depth,kdry_well,pref)

#----- Calculating elastic parameters ---------

parameters1 = ['vpBrineCo2','vsBrineCo2','ipBrineCo2','isBrineCo2','vpvsBrineCo2', \
              'vpBrineOil','vsBrineOil','ipBrineOil','isBrineOil','vpvsBrineOil', ]

differences = ['P','S','B']

co2ref=0 # 100% brine 
par = Parameters.Parameters(kdry_l,ksat_l,gdry_l,rhosat_l,ksat_loil,rhosat_loil,phi,co2,p)

par1 = par.calculation(parameters1)

m2 = len(parameters1)
m1 = len(differences)
par2 = {}
par3 = {}
for i1 in range(m1):

  string1 = "%s"%(differences[i1])
  par3[string1] =  np.zeros(((n3,n2,n1)))
  for i2 in range(m2):
    print("paso")

    string2 = parameters1[i2]+string1
    par2[string2] =  np.zeros(((n3,n2,n1)))
    property_dif,dif = par.differences(par1[parameters1[i2]],pref,co2ref,string1)
    par2[string2]= property_dif
    par3[string1] =  dif




#---Plotting static template---------------
plots.plottemplate(par1,pref,vpvs,ip,so)


#---Plotting dynamic template--------------
print("paso222")

g = par2['vpvsBrineCo2P']
print("paso3")

h = par3['P']
fig =  figure(8, figsize=(6, 5))
ax = plt.subplot(1,1,1)
#for i3 in range (n3):
#  for i2 in range(n2):
#    print("paso2")
#    ax.scatter(e[i3][i2][320],d[i3][i2][320] , c='r')
print("paso222")

for i3 in range(n3):
  for i2 in range(n2):
    print("paso2")
    ax.scatter(h[i3],g[i3][i2][320] , c='r')
      
plt.show()
