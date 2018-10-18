#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
###     Program which calculates the branching ratios of the charged              
###     Higgs to the channels cs, tn, cb, AW*. Valid for the 2HDM and              
###     MHDM. It also calculate BR(t -> H+b) x BR(H+ -> cb) 
###

### INITIAL SET UP 
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
import vegas
### INVARIANT VALUE (PHYSICS VALUE OF SM PARAMATERS)
mt = 171.2    # mass of top quark
mz = 91.18    # mass of z boson
nf = 5.0      # origin 5
PI = math.pi  # pi number 
mw = 80.33    # mass of w boson
mtau = 1.7771 # mass of tau neutrino
gf = 0.0000116639 # fermi constant
#mhch = 130.00  # mass of charged higgs
massrange = [85.00,130.00]# charged Higgs ranged values
for mhch in massrange :
   ma = 200.00 # mass of Pseudoscalar 
# CKM elements
   vcs = 0.97
   vcb = 0.04
############
#QCD running coupling constant (alp) at energy scale MH. Relevant                                                                                                                                                                        
#for QCD corrections to Higgs decay widths.
############
   coeffmc = 12.0/25.0
   coeffmb = 12.0/23.0
   alpmz = 0.1185 #coupling of z 
   alpmhch = alpmz /(1.0 + (alpmz/(12.0 * PI)) * (33.0 - 2 * nf) * \
                  math.log((mhch**2)/(mz**2))) 

   alpmb= alpmz /(1.0 + (alpmz/(12.0 * PI))*(33.0 - 2 * nf) * \
               math.log(4.89**2/(mz**2)))
   print('alpmb :',alpmb, ' alpmhch:',alpmhch)

# alpha_s at sale 2 GeV.Can read this off from Rev.Part.Properties
#alp2gev = 0.117 /(1.0 + (0.117 /(12.0 * PI)) * (33.0 - 2 * nf) *\
#                  math.log(2.0**2 / mz**2))
   alp2gev = 0.310
# for running of ms at 2 GeV to ms at 4.88 GeV (mb pole mass)                                                                                                                                                                                                                          
   fc2gev=1.0 + 1.101 * alp2gev/PI + 1.39 * (alp2gev/PI)**2 + 1.09 * \
     (alp2gev/PI)**3
   fcmb=1.0 + 1.101 * alpmb/PI + 1.39 * (alpmb/PI)**2 + 1.09 * \
     (alp2gev/PI)**3


#########################################################################
# Running quark masses at energy scale MHch as a function of running                                                                                                                                                   
# alpha strong. Pole mass of b quark is 4.88 GeV                                                                                                                                                                                          
# pole mass of c quark is 1.64 GeV. Lattice gives ms (2 GeV)=93 MeV.
   fb = 1.0 + 1.17 * alpmhch/PI + 1.50 * (alpmhch/PI)**2 + 0.1725 *(alpmhch/PI)**3
   fbmb = 1.0 + 1.17 * alpmb/PI + 1.50 * (alpmb/PI)**2 + 0.1725 *(alpmb/PI)**3

######### mass of charm quark and strange quark 
   mc = 0.9377 * (alpmhch/alpmb)**coeffmb * fb/fbmb
   mb = 4.18 * (alpmhch/alpmb)**coeffmb * fb/fbmb
##########################################################################
# mass of strange at scale of b quark pole mass(4.88 GeV)
   msatmb = 0.093 * (alpmb/alp2gev)**coeffmc * fc2gev/fcmb

######### mass of s quark at scale of charged higgs mass
   ms = msatmb * (alpmhch/alpmb)**coeffmb * fb/fbmb

#mb = 2.95# mass bottom quark
#ms =0.055
   print(mc , mb, ms)
##########################################################################
#4*pi*SQRT(2); factor appears in partial width of 2 fermions
   fac = 4.0 * math.sqrt(2.0) * PI  
#######################################################################
##################    Partial width  and Branching ratio for charged Higgs 
##################    based on different Models
#Higgs-fermion-fermion couplings in MHDM:
#couptn = Z**2
#coupc = Y**2
#coups = X**2
#coupb = X**2

# charged H - tau, neutrino(tau-type)
   def lamtn (Z):
      couptn = Z**2
      return gf * mhch * (mtau**2) * couptn / fac
# Charged H - charm, strange quark
   def lamcs (X,Y):   
      coupc = Y**2
      coups = X**2 
      return (3 * gf * mhch * (vcs**2)) * (mc**2 * coupc + ms**2 * coups) * \
           (1.0 + 17.0 * alpmz/(3.0 * PI)) / fac

# Charged H - charm, bottom quark
   def lamcb (X,Y):
      coupc = Y**2
      coupb = X**2
    
      return (3 * gf * mhch * (vcb**2)) * (mc**2 * coupc + mb**2 * coupb) * \
           (1.0 + 17.0 * alpmz/(3.0 * PI)) / fac
# BR : charged Higgs - A ,W
#braw=LamAW/(LamTN+LamCS+LamCB+LamAW) 
###################################################################################                                                                                                                                                                                                  
   if mhch > ma:
      lamaw=9.0 * gf**2 * mw**4 * mhch * gaw /(16.0 * PI**3)
   else:
      lamaw=0.00

# BR : charged Higgs - tau , neutrino(tau)
   def brtn (X,Y,Z):
    return lamtn(Z)/(lamtn(Z)+lamcs(X,Y)+lamcb(X,Y)+lamaw)
# BR : charged Higgs - charm , strange quark
   def brcs (X,Y,Z):
    return lamcs(X,Y)/(lamtn(Z)+lamcs(X,Y)+lamcb(X,Y)+lamaw)
# BR : charged Higgs - charm, bottom quark
   def brcb (X,Y,Z):
     return lamcb(X,Y)/(lamtn(Z)+lamcs(X,Y)+lamcb(X,Y)+lamaw)
###########################################################
### LOOP CALCULATION and  X, Y, Z dependence
   x = np.arange(0.0,40.2,0.2)
   y = np.arange(0.0,0.62,0.02)
   z = np.arange(0.0,5.02,0.02)
   z_x = np.array([0.1] * len(x)) # create same length of z array for x array
   z_y = np.array([0.1] * len(y)) # create same length of z array for y array
   for Y in y:
       for X in x:
#for Z in np.arange(0.0,5.1,0.5):

    
                 Z = 0.1
#          for Z in z:          


#                 print(X,Y,Z,brtn(X,Y,Z),brcs(X,Y,Z),brcb(X,Y,Z),\
#                       brtn(X,Y,Z) + brcs(X,Y,Z) + brcb(X,Y,Z))

#############################################################
###  Here evaluate the Partial  widths and Branching Ratio of t-> H+,b and t-> w,b
# Partial width for t-> H+,b

   def lamtHb(X,Y):
    return (0.5/fac) * gf * mt * (mt**2 * Y**2 + mb**2 * X**2) * \
           (1.0 - mhch**2/mt**2)**2
# Partial width for t-> w,b
   lamtWb=(0.5/fac) * gf * mt * (mt**2 + 2.0 * mw**2) * \
       (1.0 - mw**2/mt**2)**2

# Branching ratio for t-> H+,b
   def brtHb(X,Y):
    return lamtHb(X,Y) /(lamtHb(X,Y) +lamtWb)
# Branching ratio for t-> w,b
   def brtWb(X,Y):
    return lamtWb/(lamtHb(X,Y)+lamtWb)

###############################################################
### Calculate the c,b events ;tau,nv events; c,s events created from H+ in t-> H+,b 
#c,b events
   def cbevents(X,Y,Z):
    return brtHb(X,Y) * brcb(X,Y,Z)
#c,s events
   def csevents(X,Y,Z):
    return brtHb(X,Y) * brcs(X,Y,Z)
#tau,nv events
   def tnevents(X,Y,Z):
    return brtHb(X,Y) * brtn(X,Y,Z)
#c,b events + c,s events
   def totalcbcs(X,Y,Z): 
    return cbevents(X,Y,Z) + csevents(X,Y,Z)
############################################################
######   PLOTS               
##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
#now x= np.arange(0.0,20.0,0.05) build x array  (0,0.05,0.1,...,20.95)   len(y) times
#Put y= np.arange(0.0,0.5,0.01)) y array  (0.01 X 19),(0.02 X 19),(2 X 19)....(len(x) X 19)
   xarray,yarray = np.meshgrid(x,y) 
##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
# Preparing conotour plots for X,Y,Z 
###################################################

###################################################
   BRCS = brcs(xarray,yarray,0.1)
   BRCB = brcb(xarray,yarray,0.1)
   BRTN = brtn(xarray,yarray,0.1)
# set contour line value

   linecs = np.arange(0.0,1.2,0.2) # percentage limit 
   linecs1 = (0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01)   
   linecs2 = (0.001,0.005,0.01,0.02,0.05,0.10) 
   linecs3 = np.arange(-1.1,1.1,0.4) # REALXY Limit
   linecs4 = np.arange(-0.1,0.11,0.01)#IMXY Limit
   linecs5 = (0.0001,0.0005,0.001 ,0.002 ,0.005 ,0.01 )
   print('--------------------------------------------------')                                                                                                                                                                      
# Contour of Branching ratio for first H+ > c,s x-axis: |X|, y-axis = |Y|; |Z|= 0.1
#plt.figure()
#ax = plt.subplot(1, 1, 1)#axes for minor ticks plot
#ContourBRCS = plt.contour(xarray,yarray,BRCS,levels = linecs)
#plt.title('BR($H^+ \longrightarrow $ cs)')#plot title
#ax.axes.xaxis.set_minor_locator(MultipleLocator(2))# add minor ticks in x-axis
#ax.axes.yaxis.set_minor_locator(MultipleLocator(0.02))# add minor ticks in y-axis
#plt.xlabel('|X|')
#plt.ylabel('|Y|')
#plt.text(12.0, 0.2, r'|Z|= '+ str(Z) +',$M_{H^+}$= '+ str(mhch) +' GeV')# text|Z|= fixed MH+= GeV 
#plt.axis([0.0, 40.0, 0.0, 0.5])# plot x,y-axis limit : ([xmin,xmax,ymin,ymax])
#plt.clabel(ContourBRCS, inline= 0.05, fontsize=10)#plot
#plt.colorbar(ContourBRCS)
#plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#plt.savefig('hcs.png')
#plt.close()
#print('--------------------------------------------------')      
# Contour of Branching ratio for first H+ > c,b x-axis: |X|, y-axis = |Y|; |Z|= 0.1
   plt.figure()
#   ax = plt.subplot(1, 1, 1)#axes for minor ticks plot
   ContourBRCB = plt.contour(xarray,yarray,BRCB,levels = linecs)
   plt.title('BR($H^+ \longrightarrow $ cb) ') #plot title
#   ax.axes.xaxis.set_minor_locator(MultipleLocator(2))# add minor ticks in x-axis
#   ax.axes.yaxis.set_minor_locator(MultipleLocator(0.02))# add minor ticks in y-axis
   plt.xlabel('|X|')
   plt.ylabel('|Y|')
   plt.text(12.0, 0.2, r'|Z|= '+ str(Z) +',$M_{H^+}$= '+ str(mhch) +' GeV')# text|Z|= fixed MH+= GeV 
   plt.axis([0.0, 40.0, 0.0, 0.5])# plot x,y-axis limit : ([xmin,xmax,ymin,ymax])
   plt.clabel(ContourBRCB, inline=1, fontsize=10)#plot
   plt.colorbar(ContourBRCB)
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('hcb.png')
   plt.close()
   print('--------------------------------------------------')
# Contour of Branching ratio for H+ > tau,tau-neutrino x-axis: |X|, y-axis = |Y|; |Z|= 0.1
   plt.figure()
   ContourBRTN = plt.contour(xarray,yarray,BRTN,levels = linecs)
   plt.title('BR($H^+ \longrightarrow $ $\\tau \\nu_\\tau $)')#plot title
   plt.xlabel('|X|')
   plt.ylabel('|Y|')
   plt.text(12.0, 0.2, r'|Z|= '+ str(Z) +',$M_{H^+}$= '+ str(mhch) +' GeV')# text|Z|= fixed MH+= GeV 
   plt.axis([0.0, 40.0, 0.0, 0.6])# plot x,y-axis limit : ([xmin,xmax,ymin,ymax])
   plt.clabel(ContourBRTN, inline=0.05, fontsize=10)#plot 
   plt.colorbar(ContourBRTN)
   plt.grid(axis='y', linestyle='-', color='0.8') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.8') # show x-axis grid line
   plt.savefig('htaunv.png')
   plt.close()
   print('---------------------------------------------------') 
##############################################################################
##############################################################################
# Contour of Branching ratio for t >H+b - H+ > cb
   BRTHBBRCB = cbevents(xarray,yarray,0.1)
   plt.figure()
   plt.axis([0.0, 40.0, 0.0, 0.6])
   ContourBRTHBBRCB = plt.contour(xarray,yarray,BRTHBBRCB,levels = linecs2,\
     colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
   plt.plot(x, 1.1 / x,color='red',label = '|XY$^*$| < 1.1', linestyle='dashed',linewidth = 1.0)
   plt.plot(x, 0.7 / x,color='blue',label = '|XY$^*$| < 0.7', linestyle='dotted',linewidth = 1.0)
   plt.legend(loc= 'center right',frameon=False,fontsize = 'small')
   plt.title('BR($t \longrightarrow $ $H^{\pm}$b) X BR($H^{\pm} \longrightarrow $ cb)')#plot title
   plt.xlabel('|X|')
   plt.ylabel('|Y|')
   plt.text(25.0,0.4, r'|Z|='+ str(Z) +',\n $M_{H^+}$= '+ str(mhch) +' GeV')
   plt.colorbar(ContourBRTHBBRCB)
   plt.clabel(ContourBRTHBBRCB, inline=0.05, fontsize=10)#plot 
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^+}= '+ str(mhch) +' GeV,thbhcbXY.png')
   plt.show()
   print('--------------------------------------------------')
# Contour of Branching ratio for t >H+b - H+ > cs in X,Y plane
   BRTHBBRCS = csevents(xarray,yarray,0.1)
   plt.figure()
   plt.axis([0.0, 40.0, 0.0, 0.5]) 

   plt.title('BR($t \longrightarrow $ $H^{\pm}$b) X BR($H^{\pm} \longrightarrow $ cs)')#plot title
   plt.xlabel('|X|')
   plt.ylabel('|Y|')
   plt.text(25.0,0.35, r'|Z|='+ str(Z) +',\n $M_{H^+}$= '+ str(mhch) +' GeV')
   ContourBRTHBBRCS = plt.contour(xarray,yarray,BRTHBBRCS,levels = linecs2,\
        colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
   plt.plot(x, 1.1 / x,color='red',label = '|XY$^*$| < 1.1', linestyle='dashed',linewidth = 1.0)
   plt.plot(x, 0.7 / x,color='blue',label = '|XY$^*$| < 0.7', linestyle='dotted',linewidth = 1.0)
   plt.legend(loc='center right',frameon=False,fontsize = 'small')
   plt.clabel(ContourBRTHBBRCS, inline= 0.05, fontsize=10)#plot contour line
   plt.colorbar(ContourBRTHBBRCS)
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^+}= '+ str(mhch) +' GeV,thbhcsXY.png')
   plt.show()
   print('--------------------------------------------------')
# Contour of Branching ratio for t >H+b - H+ > tau,nv in Y,Z plane
   for Y in y:
          for Z in z: 
              X = 5.0
   ynum,znum = np.meshgrid(y,z)
   BRTHBBRTN = tnevents(X,ynum,znum)
   plt.figure()
   plt.axis([0.0, 0.3, 0.0, 5.0]) 
   ContourBRTHBBRTN = plt.contour(ynum,znum,BRTHBBRTN,\
                    levels = linecs5,colors = ['black','royalblue','purple','darkgreen','brown','red','black']) 
   plt.title('BR($t \longrightarrow $ $H^{\pm}$b) X BR($H^{\pm} \longrightarrow $ $\\tau \\nu_\\tau $)')#plot title
   plt.plot(np.array([0.22] * len(z)),z,color='red',label = '|XY$^*$| < 1.1', linestyle='dashed',linewidth = 1.0)
   plt.plot(np.array([0.14] * len(z)),z,color='blue',label = '|XY$^*$| < 0.7', linestyle='dotted',linewidth = 1.0)
   plt.legend(loc= 'center right',frameon=False,fontsize = 'small')
   plt.clabel(ContourBRTHBBRTN, fontsize= 8.0,fmt = '%.4f')#plot contour line
   plt.xlabel('|Y|')
   plt.ylabel('|Z|')
   plt.text(0.22, 4.0, r'|X|='+ str(X))
   plt.text(0.22, 3.5, r'$M_{H^{\pm}}$=')
   plt.text(0.22, 3.0, str(mhch) +'GeV')
   plt.colorbar(ContourBRTHBBRTN)
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,thbhtaunvYZ.png')
   plt.show()
##############################################################################
# |X| value for theta, tangentbeta, tangamma, delta (only care about first H+ this time)
# Condition d = 1 u = 2 e =3
#the = random.uniform(- math.pi/2, math.pi/2)#theta (i in loop)
#tbe = random.uniform(1.0,61.0) #tanbeta (j in loop)
#tga = random.uniform(1.0,61.0) #tangamma (k in loop)
#delta = random.uniform(0.0,2 * math.pi)#delta (l in loop)  #delta fixed#
   the = np.arange(- PI / 2, PI/50 ,PI / 50)#theta (i in loop)
   tbe = np.arange(1.0,61.0,1.0) #tanbeta (j in loop)
   tga = np.arange(1.0,61.0,1.0) #tangamma (k in loop)
   delta = np.arange(0.0,2.1 * PI , PI /6)#delta (l in loop)  #delta fixed#
   A = []
   B = []
   read1 = str('')
   read2 = str('')
   def start():
    global read1,A,B,read2
    while True:
      read1 = input('x-axis:number 0 to 3 (0 for theta, 1 for tanbeta,2 for tangamma,3 for delta):')
      if read1 == str(0) :
                 A = the 
                 print("{}".format("theta"))
                 print('A',A)
                 break
      elif read1 == str(1):
                 A = tbe
                 print("{}".format("tanbeta"))
                 print('A',A)
                 break
      elif read1 == str(2):
                 A = tga                 
                 print("{}".format("tangamma"))
                 print('A',A)
                 break
      elif read1 == str(3):
                 A = delta                 
                 print("{}".format("delta")) 
                 print('A',A)
                 break                 
      else :
                 print("Not correct variable")
    
    while True:
      read2 = input('y-axis:number 0 to 3 (0 for theta, 1 for tanbeta,2 for tangamma,3 for delta):')
      if read2 != read1:
              if read2 == str(0):
                  B = the
                  print("{}".format("theta"))
                  print('B',B)
                  break
              elif read2 == str(1):
                  B = tbe
                  print("{}".format("tanbeta"))
                  print('B',B)
                  break
              elif read2 == str(2):
                  B = tga
                  print("{}".format("tangamma"))
                  print('B',B)
                  break
              elif read2 == str(3):
                  B = delta
                  print("{}".format("delta")) 
                  print('B',B)
                  break
              else:
                  print("Not correct variable")
      else :
              print("this is already used")      
    return
   start() #choose what 2 parameters in total 4 parameters (theta,tanbeta,tangamma,delta)
#############################################################################
# create empty list to fill in all results of |X|,|Y|,|Z|
# from theta, tangentbeta, tangamma, delta
   absolutexlist =[]
   absoluteylist =[]
   absolutezlist =[]
   xyfun =[]
#############################################################################
   def U(i,j,k,l):
    sga = k / math.sqrt(k**2 + 1.0) # sinegamma
    cga = 1.0 / math.sqrt(1.0 + k**2)# cosgamma
    sthe = math.sin(i)#sinetheta
    cthe = math.cos(i)#costheta
    cde = math.cos(l)#cosdelta
    sde = math.sin(l)#sindelta
    cbe = 1.0 / math.sqrt(1.0 + j**2)# cosbeta
    sbe = j / math.sqrt(j**2 + 1.0) # sinebeta
    ud1 = sga * cbe
    ud2 = complex(- cthe * sbe * cde  - sthe * cga * cbe, cthe * sbe *sde)
    ud3 = complex(sthe * sbe * cde - cthe * cga * cbe, sthe * sbe * sde)
    uu1 = sga * sbe
    uu2 = complex(cthe * cbe * cde - sthe * cga * sbe, cthe * cbe * sde)
    uu3 = complex(- sthe * cbe * cde - cthe * cga * sbe, - sthe * cbe * sde)
    ul1 = cga
    ul2 = sthe * sga
    ul3 = cthe * sga
    return [[ud1,ud2,ud3],[uu1,uu2,uu3],[ul1,ul2,ul3]] 
##############################################################
################### d = 1 ,u = 2, l = 3 Democratic
#abs(U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0]) #X2
#abs(- U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0]) #Y2
#abs(U(i,j,k,l)[2][1] / U(i,j,k,l)[2][0]) #Z2
#abs(U(i,j,k,l)[0][2] / U(i,j,k,l)[0][0]) #X3
#abs(- U(i,j,k,l)[1][2] / U(i,j,k,l)[1][0])#Y3
#abs(U(i,j,k,l)[2][2] / U(i,j,k,l)[2][0]) #Z3
################### d = 1 ,u = 2, l = 2 Flipped
#abs(U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0]) #X2
#abs(- U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0]) #Y2
#abs(U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0]) #Z2
#abs(U(i,j,k,l)[0][2] / U(i,j,k,l)[0][0])#X3
#abs(- U(i,j,k,l)[1][2] / U(i,j,k,l)[1][0]) #Y3
#abs(U(i,j,k,l)[1][2] / U(i,j,k,l)[1][0]) #Z3
################### d = 2 ,u = 2, l = 1 Type III (Lepton-specific)
#abs(U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0]) #X2
#abs(- U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0]) #Y2
#abs(U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0]) #Z2
#abs(U(i,j,k,l)[1][2] / U(i,j,k,l)[1][0]) #X3
#abs(- U(i,j,k,l)[1][2] / U(i,j,k,l)[1][0]) #Y3
#abs(U(i,j,k,l)[0][2] / U(i,j,k,l)[0][0]) #Z3
################### d = 2 ,u = 2, l = 2 Type I
#abs(U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0]) #X2
#abs(- U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0])#Y2
#abs(U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0]) #Z2
#abs(U(i,j,k,l)[1][2] / U(i,j,k,l)[1][0]) #X3
#abs(- U(i,j,k,l)[1][2] / U(i,j,k,l)[1][0]) #Y3
#abs(U(i,j,k,l)[1][2] / U(i,j,k,l)[1][0]) #Z3
################### d = 1 ,u = 2 , l = 1 Type II
#abs(U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0]) #X2
#abs(- U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0]) #Y2
#abs(U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0]) #Z2
#abs(U(i,j,k,l)[0][2] / U(i,j,k,l)[0][0]) #X3
#abs(- U(i,j,k,l)[1][2] / U(i,j,k,l)[1][0]) #Y3
#abs(U(i,j,k,l)[0][2] / U(i,j,k,l)[0][0]) #Z3
##############################################################
   def start1():
    global read0,X2,Y2,Z2
    while True:
      read0 = input('Choose type of3HDM (1 for I, 2 for II ,\
3 for Leptonic-specific, 4 for flipped,5 for Democratic):')
      if read0 == str(1) :
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Z2
                 print('Model:',read0)
                 break
      elif read0 == str(2):
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[1][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #Z2
                 print('Model:',read0)
                 break
      elif read0 == str(3):
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #Z2
                 print('Model:',read0)
                 break
      elif read0 == str(4):
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Z2
                 print('Model:',read0)
                 break   
      elif read0 == str(5):
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[2][1] / U(i,j,k,l)[2][0] #Z2
                 print('Model:',read0)
                 break 
              
      else :
                 print("Not correct type 3HDM")
      print('Model:',read0)
    return
   start1()
#def X2(i,j,k,l):
#    return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #X2
#def Y2(i,j,k,l):
#    return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
#def Z2(i,j,k,l) :
#    return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Z2
############################################################
# |X| value for theta, tangentbeta, tangamma, delta (only care about first H+ this time)
# Condition d = 1 u = 2 e =3
   def X_function(i,j,k,l):
    sga = k / math.sqrt(k**2 + 1.0) # sinegamma
    cga = 1.0 / math.sqrt(1.0 + k**2)# cosgamma
    sthe = math.sin(i)#sinetheta
    cthe = math.cos(i)#costheta
    cde = math.cos(l)#cosdelta
    sde = math.sin(l)#sindelta
    return complex((- cthe * j * cde - sthe * cga)/sga , (- cthe * j * sde/sga))
#############################################################################
# |Y| value for theta, tangentbeta, tangamma, delta (only care about first H+ this time)
# Condition d = 1 u = 2 e =3
   def Y_function(i,j,k,l):
    sga = k / math.sqrt(k**2 + 1.0) # sinegamma
    cga = 1.0 / math.sqrt(1.0 + k**2)# cosgamma
    sthe = math.sin(i)#sinetheta
    cthe = math.cos(i)#costheta
    cde = math.cos(l)#cosdelta
    sde = math.sin(l)#sindelta
    return complex((- cthe * cde / j + sthe * cga) / sga , - cthe * sde /(j * sga))
#############################################################################   
# |Z| value for theta, tangentbeta, tangamma, delta (only care about first H+ this time)
# Condition d = 1 u = 2 e =3
   def Z_function(i,k):
       sthe = math.sin(i)#sinetheta
       return sthe * k
#############################################################################
   def complexyfunction(i,j,k,l):# (XY^*) function 
       return (X2(i,j,k,l)*(np.conj(Y2(i,j,k,l))))
#############################################################################
   print('___________________________________________')
#print(xyfun)
#BRCS1 = brcs(xarray,yarray,0.1) #|X|,|Y| for BRCS result
# (4 parameters),for BRCS result
   BRCSfinal = []
#BRCB1 = brcb(xarray,yarray,0.1)#|X|,|Y| for BRCB result
# (4 parameters),for BRCB result
#BRCBfinal = brcb(np.array(absolutexlist),np.array(absoluteylist),np.array(absolutezlist))
   BRCBfinal = []
############################################################################# 
   i = - PI / 3 #theta
   j = 40.0    #tangentbeta
   k = 10.0 #tangamma
   l = 0.0 # set delta (phase shift) to 0
   if read1 == '0':
       if read2 == '1':
          for j in B: 
              for i in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
                  
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
       elif read2 == '2':
          for k in B:
              for i in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
       else:
          for l in B:
              for i in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
   elif read1 == '1':
       if read2 == '0':
          for i in B:
              for j in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
       elif read2 == '2':
          for k in B:
              for j in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
       else:
          for l in B:
              for j in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
   elif read1 == '2':
       if read2 == '0':
          for i in B:
              for k in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
       elif read2 == '1':
          for j in B:
              for k in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
       else:
          for l in B:
              for k in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
   else:
        if read2 == '0':
          for i in B:
              for l in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
        elif read2 == '1':
          for j in B:
              for l in A:
                  
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),\
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))
        else:
          for k in B:
              for l in A:
#                  print(':',i,j,k,l,abs(X_function(i,j,k,l)))
                  absolutexlist.append(abs(X2(i,j,k,l)))
                  absoluteylist.append(abs(Y2(i,j,k,l)))
                  absolutezlist.append(abs(Z2(i,j,k,l)))
                  xyfun.append(complexyfunction(i,j,k,l))
#                  print('___________________________________________')
                  BRCBfinal.append(brcb(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
                  BRCSfinal.append(brcs(abs(X2(i,j,k,l)),abs(Y2(i,j,k,l)),abs(Z2(i,j,k,l))))
#                  print(i,j,k,l,'|',abs(X_function(i,j,k,l)),brcb(abs(X_function(i,j,k,l)),abs(Y_function(i,j,k,l)),abs(Z_function(i,k))))
#                        abs(Y_function(i,j,k,l)),abs(Z_function(i,k))\
#                        ,complexyfunction(i,j,k,l))  
   print('_____________')
#print(absolutexlist)
#print(read1,read2,A,B,np.array(xyfun))
#####################################################
#plot labels 
   readlist = ("$\\theta$","tan$\\beta$","tan$\\gamma$","$\\delta$") 
# format to 4 possible parameters
# readlist[int(read1)] = x-axis label
# readlist[int(read2)] = y-axis label
##################
#i = - PI /2#theta
#j = 2.0    #tangentbeta
#k = 40.0 #tangamma
#l = 0 # set delta (phase shift) to 0
#print('||||',brcb(abs(X_function(- PI /4,40.0,2.0,0.0)),abs(Y_function(- PI /4,40.0,2.0,0.0)),abs(Z_function(- PI /4,40.0))))
#print(i,j,k,l,BRCBfinal)
####################
   BRTN1 = brtn(xarray,yarray,0.1)#|X|,|Y| for BRTAUNV result
#4 parameters for BRTAUNV result
   BRTNfinal = brtn(np.array(absolutexlist),np.array(absoluteylist),np.array(absolutezlist))

#4 parameters for t>H+b to H+ >cb result
   BRTHBBRCB1 = cbevents(np.array(absolutexlist),np.array(absoluteylist),np.array(absolutezlist))
#4 parameters for t>H+b to H+ >cs result
   BRTHBBRCS1 = csevents(np.array(absolutexlist),np.array(absoluteylist),np.array(absolutezlist))
#4 parameters for t>H+b to H+ >cs + t>H+b to H+ >cb result
   BRTHQUARK = totalcbcs(np.array(absolutexlist),np.array(absoluteylist),np.array(absolutezlist))
#################################
#building same array length for the,tbe,delta,tga as final solution for contour plot
#print('|||',len(the),len(tbe),len(tga),len(delta),len(BRCS1),len(absolutexlist),\
#      len(np.resize(BRCSfinal,len(BRCSfinal)).reshape(len(B),len(A))))
#then,tben =np.meshgrid(the,tbe)# build same length array
#################################

   print(mhch)
#PLOT OF CONTOUR WITH RANGE OS 4 PARAMETERS INTO |X|,|Y|,|Z = 0.1| WITH CS,CB,TAUNV
   plt.figure()
#ax =plt.subplot(1,1,1)
   Contourbrcs = plt.contour(A,B, \
          np.resize(BRCSfinal,len(BRCSfinal)).reshape(len(B),len(A)),cmap = 'brg')
   plt.axis([min(A), max(A),min(B), max(B)])
#print('||',len(A),len(B),len(BRCSfinal),np.ndim(A),np.ndim(BRCSfinal))
# (4 parameters):A,B, BRCS contour plot [in the :reshape(y,x) not reshape(x,y)]
#plt.clabel(Contourbrcs, inline= 0.01, fontsize=10)# contour level show
   plt.colorbar(Contourbrcs)
   plt.title('BR($H^{\pm} \longrightarrow $ cs),$M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^+}= '+ str(mhch) +' GeV,cs.png')
   print('--------------------------------------------------')
   plt.figure()
   Contourbrcb = plt.contour(A,B, \
           np.resize(BRCBfinal,len(BRCBfinal)).reshape(len(B),len(A)),cmap = 'brg',levels = linecs)
# (4 parameters):A,B, BRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
   plt.clabel(Contourbrcb, inline= 0.02, fontsize= 9)# contour level show
   plt.colorbar(Contourbrcb)
   plt.title('BR($H^{\pm} \longrightarrow $ cb),$M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,cb.png')
   print('--------------------------------------------------')
   plt.figure()
#ax =plt.subplot(1,1,1)
   Contourbrtn = plt.contour(A,B, \
           np.resize(BRTNfinal,len(BRTNfinal)).reshape(len(B),len(A)),cmap = 'brg',levels = linecs)
#(4 parameters): A,B, BRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
#plt.clabel(Contourbrtn)# contour level show
   plt.colorbar(Contourbrtn)
   plt.title('BR($H^{\pm} \longrightarrow $ $\\tau \\nu_\\tau $),$M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,TN.png')
   print('--------------------------------------------------')
#(4 parameters): A,B, BRTHBBRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
   plt.figure()
#ax =plt.subplot(1,1,1)
#plt.axis([1.0, 4.0, 1.0, 60.0])

   ContourBRTHBBRCB1 = plt.contour(A,B, \
           np.resize(BRTHBBRCB1,len(BRTHBBRCB1)).reshape(len(B),len(A)),\
           levels = linecs2,colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
   plt.title('BR($t \longrightarrow $ $H^{\pm}$b) X BR($H^{\pm} \longrightarrow $ cb),\n $M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
   plt.clabel(ContourBRTHBBRCB1, inline= 0.01, fontsize=10)# contour level show
   plt.colorbar(ContourBRTHBBRCB1)
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,thbhcb.png')
   print('--------------------------------------------------')
#(4 parameters): A,B, BRTHBBRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
   plt.figure()
#plt.axis([1.0, 4.0, 1.0, 60.0])
   plt.title('BR($t \longrightarrow $ $H^{\pm}$b) X BR($H^{\pm} \longrightarrow $ cs),\n $M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
   ContourBRTHBBRCS1 = plt.contour(A,B, \
           np.resize(BRTHBBRCS1,len(BRTHBBRCS1)).reshape(len(B),len(A)),\
           levels = linecs2,colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
   plt.clabel(ContourBRTHBBRCS1, inline= 0.01, fontsize=10)# contour level show
   plt.colorbar(ContourBRTHBBRCS1)
#   ax.axes.xaxis.set_minor_locator(MultipleLocator(2))# add minor ticks in x-axis
#   ax.axes.yaxis.set_minor_locator(MultipleLocator(2))# add minor ticks in y-axis
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,thbhcs.png')
   print('--------------------------------------------------')
#(4 parameters): A,B, BRTHBBRCB + BRTHBBRCS contour plot
# [in the :reshape(y,x) not reshape(x,y)]
   plt.figure()
#ax =plt.subplot(1,1,1)
#plt.axis([0.0, 40.0, 0.0, 0.3])
   plt.title('BR($t \longrightarrow $ $H^{\pm}$b) X [BR($H^{\pm} \longrightarrow $ cb)\
+ BR($H^{\pm} \longrightarrow $ cs)],\n $M_{H^{\pm}}}$= '+ str(mhch) +' GeV')#plot title
   ContourTOTAL = plt.contour(A,B, \
           np.resize(BRTHQUARK,len(BRTHQUARK)).reshape(len(B),len(A)),\
           levels = linecs2,colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.clabel(ContourTOTAL, inline= 0.01, fontsize=10)# contour level show
   plt.colorbar(ContourTOTAL)
#ax.axes.xaxis.set_minor_locator(MultipleLocator(2))# add minor ticks in x-axis
#ax.axes.yaxis.set_minor_locator(MultipleLocator(2))# add minor ticks in y-axis
#ContourTOTAL = plt.contour(A,B, \
#           np.resize(BRTHQUARK,len(BRTHQUARK)).reshape(len(B),len(A)), levels = linecs2)
#plt.xlabel(readlist[int(read1)])# x-axis label
#plt.ylabel(readlist[int(read2)])# y-axis label
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,thbtotal.png')
   print('--------------------------------------------------')
#(4 parameters): A,B, |X| plot
# [in the :reshape(y,x) not reshape(x,y)]
#xarray3,yarray3 = np.meshgrid(A,B)
   plt.figure()
   plt.title('|X| in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane')#plot title
   plt.axis([0.0, max(A),min(B), max(B)])
#plt.axis([1.0, 4.0, 1.0, 60.0])
   contourx = plt.contour(A,B,\
            np.resize(np.array(absolutexlist),len(np.array(absolutexlist))).reshape(len(B),len(A)))
   plt.clabel(contourx, inline= 0.01, fontsize=10)# contour level show
   plt.colorbar(contourx)
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('x.png')
   plt.show()
   print('--------------------------------------------------')
#(4 parameters): A,B, REAL(XY^{*}) plot
# [in the :reshape(y,x) not reshape(x,y)]
#xarray3,yarray3 = np.meshgrid(A,B)
   plt.figure()
   plt.title('REAL(XY^*) in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane,$M_{H^+}$= '+ str(mhch) +' GeV')#plot title
#plt.axis([1.0, 4.0, 1.0, 60.0])
   contourxy = plt.contour(A,B,\
            np.resize(np.array(np.real(xyfun)),len(np.array(np.real(xyfun)))).reshape(len(B),len(A)),levels = np.arange(-1.1,0.9,0.6),cmap = 'brg')#,levels = linecs3)#
   plt.clabel(contourxy, inline= 0.01, fontsize=10)# contour level show
   plt.colorbar(contourxy)
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^+}= '+ str(mhch) +' GeV,realxy.png')
   print('--------------------------------------------------')
#(4 parameters): A,B, IM(XY^{*}) plot
# [in the :reshape(y,x) not reshape(x,y)]
#xarray3,yarray3 = np.meshgrid(A,B)
   plt.figure()
   plt.title('IM(XY^*) in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane,$M_{H^+}$= '+ str(mhch) +' GeV')#plot title
#plt.axis([1.0, 4.0, 1.0, 60.0])
   contourxy2 = plt.contour(A,B,\
            np.resize(np.array(np.imag(xyfun)),len(np.array(np.imag(xyfun)))).reshape(len(B),len(A)),levels = linecs4,cmap = 'brg')
#   plt.clabel(contourxy2, inline= 0.01, fontsize=10)# contour level show
#   plt.colorbar(contourxy2)
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('M{H^+}= '+ str(mhch) +' GeV,imxy.png')
   print('--------------------------------------------------')
#(4 parameters): A,B, |Y| plot
# [in the :reshape(y,x) not reshape(x,y)]
#xarray3,yarray3 = np.meshgrid(A,B)
   plt.figure()
   plt.title('|Y| in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane')#plot title
#plt.axis([1.0, 4.0, 1.0, 60.0])
   contoury = plt.contour(A,B,\
            np.resize(np.array(absoluteylist),len(np.array(absoluteylist))).reshape(len(B),len(A)))
   plt.clabel(contoury, inline= 0.01, fontsize=10)# contour level show
   plt.colorbar(contoury)
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('y.png')
   plt.show()
   print('--------------------------------------------------')
#(4 parameters): A,B, |Z| plot
# [in the :reshape(y,x) not reshape(x,y)]
#xarray3,yarray3 = np.meshgrid(A,B)
   plt.figure()
   plt.title('|Z| in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane')#plot title
#plt.axis([1.0, 4.0, 1.0, 60.0])
   contourz = plt.contour(A,B,\
            np.resize(np.array(absolutezlist),len(np.array(absolutezlist))).reshape(len(B),len(A)))
   plt.clabel(contourz, inline= 0.01, fontsize=10)# contour level show
   plt.colorbar(contourz)
   plt.xlabel(readlist[int(read1)])# x-axis label
   plt.ylabel(readlist[int(read2)])# y-axis label
   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
   plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
   plt.savefig('z.png')
   plt.show()
   print('--------------------------------------------------')

##################################################################
##################################################################
#unitary matrix (diagonalization matrix)
#### theta, tangentbeta, tangamma, delta for U^dagger
#def U(i,j,k,l):
#    sga = k / math.sqrt(k**2 + 1.0) # sinegamma
#    cga = 1.0 / math.sqrt(1.0 + k**2)# cosgamma
#    sthe = math.sin(i)#sinetheta
#    cthe = math.cos(i)#costheta
#    cde = math.cos(l)#cosdelta
#    sde = math.sin(l)#sindelta
#    cbe = 1.0 / math.sqrt(1.0 + j**2)# cosbeta
#    sbe = j / math.sqrt(j**2 + 1.0) # sinebeta
#    ud1 = sga * cbe
#    ud2 = complex(- cthe * sbe * cde  - sthe * cga * cbe, cthe * sbe *sde)
#    ud3 = complex(sthe * sbe * cde - cthe * cga * cbe, sthe * sbe * sde)
#    uu1 = sga * sbe
#    uu2 = complex(cthe * cbe * cde - sthe * cga * sbe, cthe * cbe * sde)
#    uu3 = complex(- sthe * cbe * cde - cthe * cga * sbe, - sthe * cbe * sde)
#    ul1 = cga
#    ul2 = sthe * sga
#    ul3 = cthe * sga
#    return [[ud1,ud2,ud3],[uu1,uu2,uu3],[ul1,ul2,ul3]] 
#print(U(i,j,k,l))
#############################################################

##################################################################
#################################################################
#BOSONIC DECAY H+ > A*W
#KA=(ma/mhch)**2
#KW=(mw/mhch)**2
#Gam= 4.326 / mhch**2
#p2min=0.0
#p2max=1.0 - KA
#p2 = p2min+(p2max-p2min)* 
#def f(x):
#    dx2 = 0
##        dx2 += (x[d] - 0.5) ** 2
#    return math.exp(- dx2 * 100.) * 1013.2118364296088

#integ = vegas.Integrator([[-1, 1], [0, 1], [0, 1], [0, 1]])

#result = integ(f, nitn=10, neval=1000)
#print( result.summary())
#print( 'result = %s    Q = %.2f' % (result, result.Q))