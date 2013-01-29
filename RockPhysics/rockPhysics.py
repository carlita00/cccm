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
f = [.97,.03,.0] # Mineral fraccion
#f = [.85,.15,.00] # Mineral Paluxy
#f = [1,.0,0.0] # Mineral fraccion
#f = [.7,.20,.10] # Mineral fraccion Tusc Fluvia;

rho = [2650.0, 2550.0,3960.0]
phiC = 0.4001 # critical porosity
n = 20.0 - 34.0 *phiC +14.0*np.power(phiC,2) 
c = 1.3 #TusF 1 #Tusc Beach 1.5!! # Paluxy 1.2
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
pref  = 10 # 9 # CHECKKK
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
#dmin = 3170 #3207
#dmax = 3440 #3440 #3440
#dmin2 = 3114 #3186
#dmax2 = 3301 #3440 #3301

dmin = 3207
dmax = 3440 #3440 #3440
dmin2 = 3186
dmax2 = 3301 #3440 #3301

fileName1 = './Logs/140.txt'
fileName2 = './Logs/159.txt'
#input140 = ['slowp','slows','dtc','dts','rhoD','so','phit','facies','depth','gr','bvi','cbw','vsh1']
#input159 = ['slowp','slows','dtc','dts','rhoD','so','phit','facies','depth','gr','bvi','cbw','somril']


#log140inp = {}
#log140out = {}
#log159inp = {}
#log159out = {}
#--- Reading the logs files-------

slowp1,slows1,dtc1,dts1,rhoD1,so1,phit1,facies1,depth1,gr1,bvi1,cbw1,vsh1,cal1 = Logs.logread140(fileName1)
slowp2,slows2,dtc2,dts2,rhoD2,so2,phit2,facies2,depth2,gr2,bvi2,cbw2,somril,cal2 = Logs.logread159(fileName2)


#--- Init a class each log will have one ------
log1 = Logs.Logs(fileName1,dmin,dmax,rho)
log2 = Logs.Logs(fileName2,dmin2,dmax2,rho)

#---- Calculating parameters--------
kSatC1,kSatD1,muDryC1,muDryD1,muSatC1,muSatD1,phit1,\
so1,vsh1,facies1,depth1,gr1,vpvs1,ip1,rhoCal1,rhoDD1,vp1,vpCal1,kDry1 = \
log1.parameters(slowp1,slows1,dtc1,dts1,rhoD1,so1,phit1,facies1,depth1,gr1,bvi1,cbw1,cal1)

kSatC2,kSatD2,muDryC2,muDryD2,muSatC2,muSatD2,phit2,\
so2,vsh2,facies2,depth2,gr2,vpvs2,ip2,rhoCal2,rhoDD2,vp2,vpCal2,kDry2 = \
log2.parameters(slowp2,slows2,dtc2,dts2,rhoD2,so2,phit2,facies2,depth2,gr2,bvi2,cbw2,cal2)


#print len(a)
#---- Plotting the logs ------------
#log1.plot("140-1")
#log2.plot("159-2")

#------Merge the logs----------------


depth,ksat,kdry,mudry,phit,so,gr,facies,vsh,vpvs,ip = Logs.mergelogs(kSatD2,muDryD2,phit2, \
                                                       so2,gr2,facies2,vsh2,vpvs2,ip2,depth2,cal2,rhoDD2,kDry2,kSatD1, \
                                                       muDryD1,phit1,so1,gr1,facies1,vsh1,vpvs1,ip1,depth1,cal1,rhoDD1,kDry1)


#Logs.plotqc(kSatD2,muDryD2,phit2,so2,gr2,facies2,vsh2,vpvs2,ip2,depth2,cal2,rhoDD2,rhoCal2,vp2,vpCal2, \
#            kSatD1,muDryD1,phit1,so1,gr1,facies1,vsh1,vpvs1,ip1,depth1,cal1,rhoDD1,rhoCal1,vp1,vpCal1)


#print(depth)  #n=len(kSatD2)+len(kSatD)
print(len(depth1),len(depth2),len(depth)) 


plots = Plot.Plot(kdry_l,gdry_l,kdry_c,gdry_c,kdryConstant,gdryConstant,phi,ks)
#plots.plotDry(kbri[pref],koil[pref],pref,kSatD2,mudryD2,phit2,so2,gr2,facies2,vsh2,kDry2)
#plots.plotDry(kbri[pref],koil[pref],pref,kSatD1,muDryD1,phit1,so1,gr1,facies1,vsh1,kDry1)
plots.plotDry(kbri[pref],koil[pref],pref,ksat,mudry,phit,so,gr,facies,vsh,kdry)

plots.plotDryratio('140-1',kbri[pref],koil[pref],kSatD1,muDryD1,kDry1,phit1,so1,gr1,facies1,depth1,dmin,dmax,vsh1)
#plots.plotDryratio('159-2',kbri[pref],koil[pref],kSatD2,muDryD2,kDry2,phit2,somril/100,gr2,facies2,depth2,dmin2,dmax2,vsh2)


#kdry_well = plots.plotDryratio(kbri[pref],koil[pref],kSatD2,muDryD2,phit2,so2,gr2,facies2,depth2,dmin2,dmax2,vsh2)

#plots.plotfacies(kSatD2,muDryD2,phit2,vsh2,facies2,depth2,kDry2,pref)
#kdry_well = plots.plotDryratio(kbri[pref],koil[pref],ksat,mudry,phit,so,gr,facies,depth,dmin2,dmax2,vsh)
#plots.plotfacies(kSatD2,muDryD2,phit2,vsh2,facies2,depth2,kdry_well,pref)
#plots.plotfacies(kSatD1,muSatD1,phit1,vsh1,facies1,depth1,kDry1,pref)

plots.plotfacies(ksat,mudry,phit,vsh,facies,depth,kdry,pref)

#----- Calculating elastic parameters ---------

#parameters1 = ['vpBrineCo2','vsBrineCo2','ipBrineCo2','isBrineCo2','vpvsBrineCo2', \
#              'vpBrineOil','vsBrineOil','ipBrineOil','isBrineOil','vpvsBrineOil', ]

'''
parameters = ['vpBrine','vsBrine','ipBrine','isBrine','rtBrine']
fluid = ['Co2','Oil'] 
differences = ['P','S','B']

co2ref=0 # 100% brine 
par = Parameters.Parameters(kdry_l,ksat_l,gdry_l,rhosat_l,ksat_loil,rhosat_loil,phi,co2,p)
m1 = len(fluid)
m2 = len(parameters)
m3 = len(differences)
i=0

static_par = ['']*(m1*m2)
for i1 in range(m1):
  string1 = "%s"%(fluid[i1])
  for i2 in range(m2):
    string2 = parameters[i2]
    static_par[i] = string2+string1
    i += 1

static = par.calculation(static_par) #dictionary static parameters
'''
#---Plotting static template---------------
'''
plots.plottemplate(static,pref,vpvs,ip,so)
plots.plot3dtemplate(static,pref,vpvs,ip,so,p,phi)

dynamic = {} #dictonary dynamic parameters

diff = {}    #dictionary differences

m2 = len(static_par)
i=0

dynamic_par = ['']*(m3*m2)

for i3 in range(m3):
  string3 = "%s"%(differences[i3])
  for i2 in range(m2):
    string2 = "%s"%(static_par[i2])
    string = string2+string3
    property_dif,dif = par.differences(static[string2],pref,co2ref,string3)
    dynamic[string2+string3]= property_dif
    diff[string3] =  dif
    dynamic_par[i] = string2+string3
    i += 1
print(dynamic_par)
'''

#---Plotting dynamic template--------------

#plots.plotparschange(dynamic,diff,dynamic_par,p,co2,phi) 
'''
g1 = dynamic['vsBrineCo2P']
g2 = dynamic['vpBrineCo2P']
g3 = dynamic['rtBrineCo2B']
g4 = dynamic['ipBrineCo2B']
g5 = dynamic['rtBrineOilB']
g6 = dynamic['ipBrineOilB']

#print("paso3")
'''
'''
h = diff['P']
fig =  figure(9, figsize=(6, 5))
ax = plt.subplot(1,1,1)
for i3 in range (n3):
  for i2 in range(n2):
#    print("paso2")
    ax.scatter(g4[i3][i2][320],g3[i3][i2][320] , c='r')
    ax.scatter(g6[i3][i2][320],g5[i3][i2][320] , c='g')

#print("paso222")

#for i3 in range(n3):
#  for i2 in range(n2):
    #print("paso2")
   #ax.scatter(h[i3],g1[i3][i2][320] , c='r')
  #ax.scatter(h[i3],g2[i3][0][320] , c='b')
  #ax.scatter(h[i3],g2[i3][1][320] , c='b')
  #ax.scatter(h[i3],g2[i3][2][320] , c='b')
  #ax.scatter(h[i3],g2[i3][3][320] , c='b')
  #ax.scatter(h[i3],g2[i3][4][320] , c='b')

'''

      
plt.show()
