import ElasticMedia as ElasticMedia
import GranularMedia as GranularMedia
import PlotTom as Plot
import elastic_parameters as ep
import matplotlib.pylab as pylab
import matplotlib.mlab as mll

from pylab import *
import numpy as np
import matplotlib.pyplot as plt

#Tom, here is a very useful. I hope it helps you#

#-----Global variables--------
  

k = [37.0e09,25.0e09,123.0e09]    # Bulk modulus for Quartz, Clay, Siderite (Pa)  
g = [44.e09,9.0e09,51.0e09]       # Shear modulus for Quartz, Clay, Siderite (Pa)  

f = [.7,.20,.10]                  # Mineral fraccion Tuscaloosa Fluvial;

rho = [2650.0, 2550.0,3960.0]     # Density for Quartz, Clay, Siderite (Pa)  


phiC = 0.4001                     # critical porosity
n = 20.0 - 34.0 *phiC +14.0*np.power(phiC,2) # coordination number caculation
pref = 9                          # Effective pressure reference, corresponds to vector position
#-------Elastic Media------------
  
ks  = ElasticMedia.hillAverage(k,f)
gs  = ElasticMedia.hillAverage(g,f)
prs = ep.poisson_value(ks,gs)
kc = k[0]
gc = g[0]
prc = ep.poisson_value(kc,gc)
mc  = ep.M_value(kc,gc)

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
gdry_c = gm.gdryCement(c,phi,gc,kdry_c)

plots = Plot.Plot(kdry_l,gdry_l,kdry_c,gdry_c,phi)
