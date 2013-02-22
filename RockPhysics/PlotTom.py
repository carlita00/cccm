#!/usr/bin/env python
import Gassman as Gassman
import math
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mll

#---------Plotting----------------
class Plot:       # Here I am defining the Plot class. The variables inside init  
                  # are global meaning that you can use them in every function
                  # that you create inside this class. That is the purpose of the
                  # self. 

  def __init__(self,kdry_l,gdry_l,kdry_c,gdry_c):
    self._kdry_l = kdry_l
    self._kdry_c = kdry_c
    self._gdry_l = gdry_l
    self._gdry_c = gdry_c
    

  def plot(self,phi,pref):

    plt.rcParams.update({'font.size': 14,'legend.fontsize': 14}) # Letter size of the plots 

    fig = figure(1, figsize=(13, 4))
    ax = plt.subplot(1,2,1)
    ax.set_ylabel('K dry (Pa)')
    plot(self._phi[:],self._kdry_l[pref,:],'b')
    plot(self._phi[:],self._kdry_c[:],'b')
    ax.set_xlabel('porosity $\phi$')
    ylim(1.e9,1.5e10) 
    xlim(0,0.5) 


    bx = plt.subplot(1,2,2)
    bx.set_xlabel('porosity $\phi$')
    bx.set_ylabel('G dry (Pa)')
    plot(self._phi[:],self._gdry_l[pref,:],'b')
    plot(self._phi[:],self._gdry_c[:],'b')
    ylim(1e9,1.5e10) 
    xlim(0,0.5) 

    fig.savefig('example.png', bbox_inches='tight')


plt.show()

