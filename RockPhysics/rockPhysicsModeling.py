import sys
from java.awt import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from Modules import *
import elastic_parameters as ep

def main(args):

  #-----Global variables--------
  
  k = [37.0e09,25.0e09,123.0e09]    #Pa
  
  g = [44.e09,9.4e09,51.0e09]       #Pa
  
  f = [.8,.19,.01] # Mineral fraccion
  
  rho = [2650.0, 2550.0,3960.0]
  
  phiC = 0.4001 # critical porosity
  
  n = 20.0 - 34.0 *phiC +14.0*pow(phiC,2) 
  
  n *=1.3
  
  print(n)
  
  
  
  #-------Elastic Media------------
  
  em = ElasticMedia()
  
  ks  = em.hillAverage(k,f)
  
  gs  = em.hillAverage(g,f)
  
  prs = ep.poisson(ks,gs)
  
  prc = ep.poisson(k[1],g[1])
  
  mc  = ep.M_value(k[1],g[1])
  
  #--------Elastic Modeling---------

  p = arange(3.5e06,17.5e06,1.0e06)
  
  phi = arange(0.0001,0.4001,0.001)

  '''Hertz mindlin''' 
  
  gm = GranularMedia(n,phiC,gs,prs,p)
  
  khm = gm.applyKhm() 

  ghm = gm.applyGhm() 

  zhm = gm.applyZhm() 

  z = gm.applyZ() 

  '''Modified lower Hashin-Shtrikman bound'''
  
  '''Modified upper Hashin-Shtrikman bound'''

  CO2 = arange(0.0,1.0,.1)

def arange(init,fin,sample):
  n = (fin-init)/sample + 1
  x = zerofloat(int(n))
  x[0] = init
  for i in range(1,n):
    x[i] = init + i*sample
  return x

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
