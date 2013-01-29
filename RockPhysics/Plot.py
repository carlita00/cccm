#!/usr/bin/env python
import Gassman as Gassman
import numpy as np
import math
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#---------Plotting----------------
class Plot:
  def __init__(self,kdry_l,gdry_l,kdry_c,gdry_c,kdryConstant,gdryConstant,phi,ks):
    self._kdry_l = kdry_l
    self._kdry_c = kdry_c
    self._gdry_l = gdry_l
    self._gdry_c = gdry_c
    self._kdryConstant = kdryConstant
    self._gdryConstant = gdryConstant
    self._ks = ks
    self._phi = phi
    

  def plotDry(self,kbri,koil,pref,ksat_well,gdry,phit,SO,gr,facies,vsh,kdry_well):


    plt.rcParams.update({'font.size': 16,'legend.fontsize': 14})

    #phid,phit,ksat_well,vsh_iw,g,SO,depth,well,facies,sid = np.loadtxt(fileName,unpack=True)
    #kdry_well = Gassman.kdry(ksat_well,phit,SO,vsh,kbri,koil)

    n1 = len(ksat_well)
    dry_ratio = np.zeros(n1)
    #for i in range(n1):  
      #print(vsh[i])
    fig = figure(4, figsize=(7, 11))
    ax = plt.subplot(2,1,1)
    #ax.set_xlabel('porosity $\phi$')
    ax.set_ylabel('K dry (Pa)')
    #plot(self._phi[:],self._kdry_l[1,:],'b')
    #plot(self._phi[:],self._kdry_l[2,:],'b')
    plot(self._phi[:],self._kdry_l[pref,:],'b')
    plot(self._phi[:],self._kdry_c[:],'b')
    #plot(self._phi[:],self._kdryConstant[3,:],'b')
    plot(self._phi[:],self._kdryConstant[5,:],'b')
    plot(self._phi[:],self._kdryConstant[7,:],'b')
    plot(self._phi[:],self._kdryConstant[9,:],'b')
    ax.set_xticklabels([])

    #plot(self._phi[:],self._kdryConstant[18,:],'b')
    ylim(1.e9,1.5e10) 
    xlim(0,0.5) 
    p=ax.scatter(phit,kdry_well,c=vsh)
    #p=ax.scatter(phit,ksat_well,c='r')
    p.set_clim([0,.5])
    cbar= plt.colorbar(p)
    cbar.set_label('Vsh')
    #plt.rcParams.update({'font.size': 12,'legend.fontsize': 12})

    bx = plt.subplot(2,1,2)
    bx.set_xlabel('porosity $\phi$')
    bx.set_ylabel('G dry (Pa)')
    plot(self._phi[:],self._gdry_l[pref,:],'b')
    plot(self._phi[:],self._gdry_c[:],'b')
    plot(self._phi[:],self._gdryConstant[5,:],'b')
    plot(self._phi[:],self._gdryConstant[7,:],'b')
    plot(self._phi[:],self._gdryConstant[9,:],'b')
    ylim(1e9,1.5e10) 
    xlim(0,0.5) 

    p=bx.scatter(phit,gdry,c=vsh)
    cbar= plt.colorbar(p)
    cbar.set_label('Vsh')
    p.set_clim([0,.5])
    fig.savefig('model_calibration.pdf')


    fig = figure(5, figsize=(6, 6))
    ax = plt.subplot(1,1,1)
    ax.set_xlabel('porosity $\phi$')
    ax.set_ylabel('K dry (Pa)')
    #plot(self._phi[:],self._kdry_l[1,:],'b')
    #plot(self._phi[:],self._kdry_l[2,:],'b')
    plot(self._phi[:],self._kdry_l[pref,:],'b')
    plot(self._phi[:],self._kdry_c[:],'b')
    #plot(self._phi[:],self._kdryConstant[2,:],'b')
    plot(self._phi[:],self._kdryConstant[5,:],'b')
    plot(self._phi[:],self._kdryConstant[7,:],'b')
    plot(self._phi[:],self._kdryConstant[9,:],'b')
    ylim(1.e9,4e10) 
    xlim(0,0.5) 
    fig.savefig('teory_model.pdf')


    #plt.rcParams.update({'font.size': 20,'lengend.fontsize':12})

  def plotDryratio(self,title,kbri,koil,ksat_well,gdry,kdry_well,phit,SO,gr,facies,depth,dmin,dmax,vsh):
    #phid,phit,ksat_well,vsh_iw,g,SO,depth,well,facies,sid = np.loadtxt(fileName,unpack=True)
    #kdry_well= Gassman.kdry(ksat_well,phit,SO,vsh,kbri,koil)
    n1 = len(ksat_well)
    dry_ratio = np.zeros(n1)
    for i in range(n1):  
      #print(vsh[i],depth[i])  
      dry_ratio[i]=(kdry_well[i]+4.0/3.0*gdry[i])/gdry[i] 
      dry_ratio[i]=np.power(dry_ratio[i],0.5)




    fig = figure(7, figsize=(16, 8))
    fig.suptitle(title)

    bx = plt.subplot(1,5,1)
    bx.set_xlabel('Gamma Ray')
    bx.set_ylabel('Depth (feet)')
    bx.xaxis.set_major_locator(MaxNLocator(4))
    p=bx.plot(gr,depth,c='k')
    bx.set_ylim(dmax,dmin) 
    #bx.set_title('')
    #bx.set_yticklabels([])


    colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0]
    n = len(depth)
    for i in range(n):
      for j in range(len(colors)):
        if facies[i] == colors[j]:
          if colors[j] >=5 and colors[j]<=9:
            j = 4
            c1='k'
            p=bx.scatter(gr[i],depth[i],c=c1)
          elif colors[j] >=50 and colors[j]<=90:
            j = 13
            c1='k'
            p=bx.scatter(gr[i],depth[i],c=c1)
          else:
            c1=cm.hsv((j+1)/float(len(colors)),1)
            p=bx.scatter(gr[i],depth[i],c=c1,edgecolor='none')
    bx.set_xlim(0,160)


    bx = plt.subplot(1,5,2)
    bx.set_xlabel('Vp/Vs dry ratio')
    #bx.set_ylabel('Depth (ft)')
    #p=bx.scatter(dry_ratio,depth,c=SO,edgecolor='none')
    p1=bx.plot(dry_ratio,depth,c='k')
    bx.xaxis.set_major_locator(MaxNLocator(4))
    #p.set_clim([0,0.5])
    #cbar= plt.colorbar(p)
    #cbar.set_label('Oil saturation')
    ylim(dmax,dmin) 
    bx.set_yticklabels([])

    xlim(1.4,1.8) 
    '''
    bx = plt.subplot(1,5,3)
    bx.set_xlabel('Oil saturation')
    #bx.set_ylabel('Depth (ft)')
    #p=bx.scatter(dry_ratio,depth,c=SO,edgecolor='none')
    p1=bx.plot(SO,depth,c='k')
    bx.xaxis.set_major_locator(MaxNLocator(4))
    #p.set_clim([0,0.5])
    #cbar= plt.colorbar(p)
    #cbar.set_label('Oil saturation')
    bx.set_yticklabels([])

    ylim(dmax,dmin) 
    xlim(0,0.5) 
    '''
    bx = plt.subplot(1,5,4)
    bx.set_xlabel('Kdry (Pa)')
    #bx.set_ylabel('Depth (ft)')
    #p=bx.scatter(dry_ratio,depth,c=SO,edgecolor='none')
    p1=bx.plot(kdry_well,depth,c='k')
    bx.set_yticklabels([])

    bx.xaxis.set_major_locator(MaxNLocator(4))
    #p.set_clim([0,0.5])
    #cbar= plt.colorbar(p)
    #cbar.set_label('Oil saturation')
    ylim(dmax,dmin) 
    bx.set_xlim(1e9 ,2.e10) 

    kbri=2.44e9
    koil=0.04e9
    rhobri=990
    rhooil=0.65e3
    #kgassman(kdry_well,phit,SO,vsh,kbri,koil,rhobri,rhooil)
    ksatR,rhosatR,soR=Gassman.kgassman(kdry_well,phit,SO,vsh,kbri,koil,rhobri,rhooil)


    bx = plt.subplot(1,5,3)
    bx.set_xlabel('Oil saturation')
    #bx.set_ylabel('Depth (ft)')
    #p=bx.scatter(dry_ratio,depth,c=SO,edgecolor='none')
    #p1=bx.plot(SO,depth,c='k')
    p2=bx.plot(soR,depth,c='k')
    bx.xaxis.set_major_locator(MaxNLocator(4))
    #p.set_clim([0,0.5])
    #cbar= plt.colorbar(p)
    #cbar.set_label('Oil saturation')
    bx.set_yticklabels([])

    ylim(dmax,dmin) 
    xlim(0,0.5) 

    bx = plt.subplot(1,5,5)
    bx.set_xlabel('Ksat (Pa)')
    #bx.set_ylabel('Depth (ft)')
    p=bx.plot(ksat_well,depth,c='k',label="Ksat Obs")
    p1=bx.plot(ksatR,depth,c='r',label="Ksat Cal")
    xlim(3.e9,2.5e10) 
    bx.xaxis.set_major_locator(MaxNLocator(4))
    ylim(dmax,dmin) 
    bx.set_yticklabels([])
    plt.legend()
    fig.savefig(title+'QC2.pdf')


    fig = figure(8, figsize=(8, 8))
    fig.suptitle(title)
    bx = plt.subplot(1,2,1)
    bx.set_xlabel('Oil saturation')
    bx.set_ylabel('Depth (ft)')
    #p=bx.scatter(dry_ratio,depth,c=SO,edgecolor='none')
    #p1=bx.plot(SO,depth,c='k')
    p2=bx.plot(soR,depth,c='k')
    bx.xaxis.set_major_locator(MaxNLocator(4))
    #p.set_clim([0,0.5])
    #cbar= plt.colorbar(p)
    #cbar.set_label('Oil saturation')
    #bx.set_yticklabels([])

    ylim(dmax,dmin) 
    xlim(-0.1,0.5) 


    bx = plt.subplot(1,2,2)
    bx.set_xlabel('Ksat (Pa)')
    #bx.set_ylabel('Depth (ft)')
    p=bx.plot(ksat_well,depth,c='k',label="Ksat Obs")
    p1=bx.plot(ksatR,depth,c='r',label="Ksat Cal")
    xlim(3.e9,2.5e10) 
    bx.xaxis.set_major_locator(MaxNLocator(4))
    ylim(dmax,dmin) 
    bx.set_yticklabels([])
    plt.legend()
    fig.savefig(title+'fixed.pdf')


  def plotfacies(self,ksat_well,gdry,phit,vsh_iw,facies,depth,kdry_well,pref):
    plt.rcParams.update({'font.size': 16,'legend.fontsize': 14})

    fig = figure(9, figsize=(7, 11))
    bx = plt.subplot(2,1,1)
    #bx.set_xlabel('porosity $\phi$')
    bx.set_ylabel('K Dry (Pa)')
    #axis([0.0, 0.4, -10, 4*1e10])
    #ylim(0,4*1e10)
    ylim(1e9,1.5e10) 
    xlim(0,0.5) 
    bx.set_xticklabels([])

    colors=[1.0,2.,3.,4.,5.,6.,7.,8.,9.,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0]
    n = len(depth)
    for i in range(n):
      for j in range(len(colors)):
        if vsh_iw[i] <50:
          if facies[i] == colors[j]:
            if colors[j] >=5 and colors[j]<=9:
              j = 4
              c='k'
              p=bx.scatter(phit[i],kdry_well[i],color=c)
            elif colors[j] >=50 and colors[j]<=90:
              j = 13
              c = 'k'
              p=bx.scatter(phit[i],kdry_well[i],color=c)
            else:
              c=cm.hsv((j+1)/float(len(colors)),1)
              p=bx.scatter(phit[i],kdry_well[i],color=c)

            #print(facies[i])
    plot(self._phi[:],self._kdry_l[pref,:],c='b')
    plot(self._phi[:],self._kdry_c[:],c='b')
    #plot(self._phi[:],self._kdryConstant[3,:],c='b')
    plot(self._phi[:],self._kdryConstant[5,:],c='b')
    plot(self._phi[:],self._kdryConstant[7,:],c='b')
    plot(self._phi[:],self._kdryConstant[9,:],c='b')


    bx = plt.subplot(2,1,2)
    bx.set_xlabel('porosity $\phi$')
    bx.set_ylabel('G dry (Pa)')
    plot(self._phi[:],self._gdry_l[pref,:],'b')
    plot(self._phi[:],self._gdry_c[:],'b')
    plot(self._phi[:],self._gdryConstant[5,:],'b')
    plot(self._phi[:],self._gdryConstant[7,:],'b')
    plot(self._phi[:],self._gdryConstant[9,:],'b')
    ylim(1e9,1.5e10) 
    xlim(0,0.5) 

    for i in range(n):
      for j in range(len(colors)):
        if vsh_iw[i] <50:
          if facies[i] == colors[j]:
            if colors[j] >=5 and colors[j]<=9:
              j = 4
              c='k'
              p=bx.scatter(phit[i],gdry[i],color=c)
            elif colors[j] >=50 and colors[j]<=90:
              j = 13
              c = 'k'
              p=bx.scatter(phit[i],gdry[i],color=c)
            else:
              c=cm.hsv((j+1)/float(len(colors)),1)
              p=bx.scatter(phit[i],gdry[i],color=c)




    #p=bx.scatter(d[pref][9][:],c[pref][9][:],c=so)
    p1 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((1)/float(len(colors)),1))
    p2 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((2)/float(len(colors)),1))
    p3 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((3)/float(len(colors)),1))
    p4 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((4)/float(len(colors)),1))
    p5 = Rectangle((0, 0), 0.2, 0.2, fc='k')

    p20 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((11)/float(len(colors)),1))
    p30 = Rectangle((0, 0), 0.2, 0.2, fc=cm.hsv((12)/float(len(colors)),1))
    p40 = Rectangle((0, 0), 0.2, 0.2, fc='k')



    labels = ('Tusc: Beach / Barrier bar' , 'Tusc: Washover', 'Tusc: Transitional','Tusc: Fluvial', \
          'Tusc: Poor rock', 'Paluxy: Distributary sand','Paluxy: mediocre distributary sand', \
          'Paluxy: Poor rock')
    legend = plt.legend([p1,p2,p3,p4,p5,p20,p30,p40],labels, loc='center left', bbox_to_anchor=(1, 1),ncol=1)
    ltext = gca().get_legend().get_texts()
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig('model_facies.pdf',bbox_extra_artists=(legend,), bbox_inches='tight')

    #---- Tuscaloosa fluvial modeling---------
    '''
    fig = figure(10, figsize=(7, 6))
    bx = plt.subplot(1,1,1)
    bx.set_title('Quartz:70%, Clay:20%, Sid:10% ')
    bx.set_xlabel('porosity $\phi$')
    bx.set_ylabel('K Dry (Pa)')
    #axis([0.0, 0.4, -10, 4*1e10])
    #ylim(0,4*1e10)
    ylim(1e9,1.5e10) 
    xlim(0,0.5) 
    #bx.set_xticklabels([])
    plot(self._phi[:],self._kdry_l[pref,:],c='r')
    plot(self._phi[:],self._kdry_c[:],c='b')
    #plot(self._phi[:],self._kdryConstant[3,:],c='b')
    plot(self._phi[:],self._kdryConstant[5,:],c='b')
    plot(self._phi[:],self._kdryConstant[7,:],c='b')
    plot(self._phi[:],self._kdryConstant[9,:],c='b')

    n = len(depth)
    for i in range(n):
      for j in range(len(colors)):
        if facies[i]==4 :
          if facies[i] == colors[j]:
            c=cm.hsv((j+1)/float(len(colors)),1)
            p=bx.scatter(phit[i],kdry_well[i],color=c)

    fig.savefig('facies_modelingTF.pdf')
    '''
    

    #---- Tuscaloosa beach modeling---------
    
    fig = figure(10, figsize=(7, 6))
    bx = plt.subplot(1,1,1)
    bx.set_title('Quartz:97%, Clay:3%')
    bx.set_xlabel('porosity $\phi$')
    bx.set_ylabel('K Dry (Pa)')
    #axis([0.0, 0.4, -10, 4*1e10])
    #ylim(0,4*1e10)
    ylim(1e9,1.5e10) 
    xlim(0,0.5) 
    #bx.set_xticklabels([])
    plot(self._phi[:],self._kdry_l[pref,:],c='b')
    plot(self._phi[:],self._kdry_c[:],c='b')
    #plot(self._phi[:],self._kdryConstant[3,:],c='b')
    plot(self._phi[:],self._kdryConstant[5,:],c='r')
    plot(self._phi[:],self._kdryConstant[7,:],c='b')
    plot(self._phi[:],self._kdryConstant[9,:],c='b')

    n = len(depth)
    for i in range(n):
      for j in range(len(colors)):
        if facies[i]==1 :
          if facies[i] == colors[j]:
            c=cm.hsv((j+1)/float(len(colors)),1)
            p=bx.scatter(phit[i],kdry_well[i],color=c)

    fig.savefig('facies_modelingTB.pdf')
    

    
    #---- Paluxy-Tusc modeling---------
    '''
    fig = figure(10, figsize=(7, 6))
    bx = plt.subplot(1,1,1)
    bx.set_title('Quartz:85%, Clay:15%')
    bx.set_xlabel('porosity $\phi$')
    bx.set_ylabel('K Dry (Pa)')
    #axis([0.0, 0.4, -10, 4*1e10])
    #ylim(0,4*1e10)
    ylim(1e9,1.5e10) 
    xlim(0,0.5) 
    #bx.set_xticklabels([])
    plot(self._phi[:],self._kdry_l[pref,:],c='b')
    plot(self._phi[:],self._kdry_c[:],c='b')
    #plot(self._phi[:],self._kdryConstant[3,:],c='b')
    plot(self._phi[:],self._kdryConstant[10,:],c='b')
    plot(self._phi[:],self._kdryConstant[6,:],c='r')
    plot(self._phi[:],self._kdryConstant[8,:],c='b')

    n = len(depth)
    for i in range(n):
      for j in range(len(colors)):
        if facies[i]>=20 and facies[i]<=40 :
          if facies[i] == colors[j]:
            c=cm.hsv((j+1)/float(len(colors)),1)
            p=bx.scatter(phit[i],kdry_well[i],color=c)
        if facies[i]>=2 and facies[i]<=3 :
          if facies[i] == colors[j]:
            c=cm.hsv((j+1)/float(len(colors)),1)
            p=bx.scatter(phit[i],kdry_well[i],color=c)

    fig.savefig('facies_modelingP.pdf')
    '''

  def plottemplate(self,par,pref,vpvs,ip,so):
    fig = figure(10, figsize=(10, 8))
    bx = plt.subplot(1,1,1)
    bx.set_xlabel('Ip')
    bx.set_ylabel('vp/vs')
    a= par['rtBrineCo2']
    b= par['ipBrineCo2']
    c= par['rtBrineOil']
    d= par['ipBrineOil']


    p=bx.scatter(b[pref][0][:],a[pref][0][:],c='b')
    p=bx.scatter(b[pref][9][:],a[pref][9][:],c='r')
    p=bx.scatter(d[pref][3][:],c[pref][3][:],c='g')
    #p=bx.scatter(d[pref][9][:],c[pref][9][:],c=so)
    p=bx.scatter(ip,vpvs,c=so)


    cbar= plt.colorbar(p)
    cbar.set_label('Oil Saturation')
    #plot(self._phi[:],self._gdry_c[:],'b')
    ylim(1,3) 
    #xlim(0,0.5) 
    p.set_clim([0,0.3])

  def plot3dtemplate(self,par,pref,vpvs,ip,so,p,phi):
    fig = figure(11, figsize=(10, 8))
    #bx = plt.subplot(1,1,1)
    bx = fig.add_subplot(111, projection='3d')
    bx.set_xlabel('Ip')
    bx.set_ylabel('vp/vs')
    bx.set_zlabel('Pe')

    a= par['rtBrineCo2']
    b= par['ipBrineCo2']
    c= par['rtBrineOil']
    d= par['ipBrineOil']
    n2 = len(p)
    n1 = len(phi)
    n = len(ip)
    plotPw = np.zeros(n)
    for i1 in range(n):
      plotPw[i1] = p[pref] 
   
    plotP = np.zeros((n2,n1))
    for i2 in range(n2):
      for i1 in range(n1):
        plotP[i2][i1] = p[i2] 
    for i in range(n2):
      #print(len(b[pref][0][:]),len(p[i]))
      bx.scatter(b[i][0][:],a[i][0][:],plotP[i][:],c='b')
      bx.scatter(b[i][9][:],a[i][9][:],plotP[i][:],c='r')

    ##p=bx.scatter(b[pref][9][:],a[pref][9][:],c='r')
    #p=bx.scatter(d[pref][3][:],c[pref][3][:],c='g')
    #p=bx.scatter(d[pref][9][:],c[pref][9][:],c=so)
    #p1=bx.scatter(ip,vpvs,plotPw,c=so)
    p1=bx.scatter(ip,vpvs,plotPw,c=so)


    cbar= plt.colorbar(p1)
    cbar.set_label('Oil Saturation')
    #plot(self._phi[:],self._gdry_c[:],'b')
    ylim(1,3) 
    #xlim(0,0.5) 
    p1.set_clim([0,0.3])
    plt.show()

  def plotparschange(self,dynamic,diff,dynamic_par,p,co2,phi):
    m1 = len(dynamic)
    n1 = len(phi)
    n2 = len(co2)
    n3 = len(p)
    fig = figure(12, figsize=(10, 8))
    bx = plt.subplot(1,1,1)
    bx.set_xlabel('Difference in P (Pa)')
    bx.set_ylabel('Parameters')
    colors=['y','g','r','b','k','m' ]
    i=0
    for j1 in range(m1):
      a = dynamic_par[j1]
      if(a[10] == 'P' ):
        if(a[7]== 'C'):
          c = dynamic[dynamic_par[j1]]
          d = diff['P']
          i=i+1

          #print(dynamic_par[j1], i)

          for i1 in range(n3):
            #print(c[i1][0][320],dynamic_par[j1], i)
            bx.scatter(d[i1],c[i1][0][320] , c=colors[i])
            #i=i+1

  


plt.show()

