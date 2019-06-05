#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 16:27:02 2018

@author: apple
"""
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
#import vegas
import eehh as fl
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from invariant import *

##########################################################################
#4*pi*SQRT(2); factor appears in partial width of 2 fermions
fac = 4.0 * np.sqrt(2.0) * PI 
#############################
#the = random.uniform(- math.pi/2, math.pi/2)#theta (i in loop)
#tbe = random.uniform(1.0,61.0) #tanbeta (j in loop)
#tga = random.uniform(1.0,61.0) #tangamma (k in loop)
#delta = random.uniform(0.0,2 * math.pi)#delta (l in loop)  #delta fixed#
the = np.arange(- PI / 2, PI/50 ,PI / 20)#theta (i in loop)
tbe = np.arange(1.0,61.0,1.0) #tanbeta (j in loop)
tga = np.arange(1.0,61.0,1.0) #tangamma (k in loop)
delta = np.arange(0.0,2.1 * PI , PI /6)#delta (l in loop)  #delta fixed#
A = []
B = []
read1 = str('')
read2 = str('')
i = - PI / 2.5 #theta
j = 40.0    #tangentbeta
k = 10.0 #tangamma
l = 0.0 # set delta (phase shift) to 0
x = np.arange(0.0,40.2,0.2) # x range
y = np.arange(0.0,0.62,0.02) # y range
z = np.arange(0.0,5.02,0.02) # z range
# set contour line value
linecs = np.arange(0.0,1.2,0.2) # percentage limit 
linecs1 = (0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01)   
linecs2 = (0.001,0.005,0.01,0.02,0.05,0.10) 
linecs3 = np.arange(-1.1,1.1,0.4) # REALXY Limit
linecs4 = np.arange(-0.1,0.11,0.01)#IMXY Limit
linecs5 = (0.0001,0.0005,0.001 ,0.002 ,0.005 ,0.01 )
print('--------------------------------------------------')   
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
####################################################
# BR : charged Higgs - tau , neutrino(tau)
def brtn (X,Y,Z):
       return lamtn(Z)/(lamtn(Z)+lamcs(X,Y)+lamcb(X,Y)+lamaw)
# BR : charged Higgs - charm , strange quark
def brcs (X,Y,Z):
       return lamcs(X,Y)/(lamtn(Z)+lamcs(X,Y)+lamcb(X,Y)+lamaw)
# BR : charged Higgs - charm, bottom quark
def brcb (X,Y,Z):
       return lamcb(X,Y)/(lamtn(Z)+lamcs(X,Y)+lamcb(X,Y)+lamaw)
# BR : charged Higgs - A ,W
#braw=LamAW/(LamTN+LamCS+LamCB+LamAW) 
############################################################################
#############################################################
###  Here evaluate the Partial  widths and Branching Ratio of t-> H+,b and t-> w,b
# Partial width for t-> H+,b
def lamtHb(X,Y):
       return (0.5/fac) * gf * mt * (mt**2 * Y**2 + mb**2 * X**2) * \
           (1.0 - mhch**2/mt**2)**2
# Partial width for t-> w,b
lamtWb = (0.5/fac) * gf * mt * (mt**2 + 2.0 * mw**2) * \
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
################################################
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
#unitary matrix (diagonalization matrix)
#### theta, tangentbeta, tangamma, delta for U^dagger
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
###################################################################
def start():#choose what 2 parameters in total 4 parameters (theta,tanbeta,tangamma,delta)
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

#############################################################################
def start1():# choose model
        global X2,Y2,Z2
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
#############################################################################
def complexyfunction(i,j,k,l):# (XY^*) function 
        return X2(i,j,k,l) * (np.conjugate(Y2(i,j,k,l)))
###############################################################################
def start3():
#        global X2,Y2,Z2 #
#        print('x2',X2)
        reference_array = [i,j,k,l]
        read1_int = int(read1)
        read2_int = int(read2)
        for var_b in B:
            for var_a in A:
                my_tuple = ()
                for counter in range(4):
                    if read1_int == counter:
                        my_tuple += (var_a,)
                    elif read2_int == counter:
                        my_tuple += (var_b,)
                    else:
                        my_tuple += (reference_array[counter],)
                absolutexlist.append(abs(X2(*my_tuple))) # unpacked bracket of my_tuple by * sign
                absoluteylist.append(abs(Y2(*my_tuple)))
                absolutezlist.append(abs(Z2(*my_tuple)))
                xyfun.append(complexyfunction(*my_tuple))
                yfun.append(Y2(*my_tuple))
                BRCBfinal.append(brcb(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple))))
#                BRCB2final.append(brcb(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple)))**2)
                BRCSfinal.append(brcs(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple))))
                BRTNfinal.append(brtn(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple))))
        print('lenbrcb',len(BRCBfinal))
        BRCBfinal_list.append(BRCBfinal)
        BRCSfinal_list.append(BRCSfinal)
        BRTNfinal_list.append(BRTNfinal)
        xyfun_list.append(xyfun)
        yfun_list.append(yfun)
        return
def max_valueposition(xxx):# x has to be np.array ; The max value postion of long array
    yyy = np.where(xxx == max(xxx))[0]
    return int(yyy)
def mhch_sig4jtag(x,y):
    y_value = []
    print('xlen',len(x))
    for n in np.arange(0,len(x)):
        fourjet = (fl.eeHHcbcb_2bsignal(x[n],x[n]) + fl.eeHHcbcb_1bsignal(x[n],x[n]) + \
                   fl.eeHHcbcs_2bsignal(x[n],y[n]) + fl.eeHHcbcs_1bsignal(x[n],y[n]) + \
                   fl.eeHHcbcb_0bsignal(x[n],x[n]) + fl.eeHHcbcs_0bsignal(x[n],y[n]) + \
                   fl.eeHHcscs_0bsignal(y[n],y[n])) * fl.epsilon / np.sqrt(fl.backgroundtagging()) 
#        print(len(fourjet))
        y_value.append(fl.eeHH_event() * fourjet)
    print('y_value',len(y_value))
    return y_value
#######################################################################
#######################################################################
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
#   print('BRCBfinal',BRCBfinal)
####################
BRCBfinal_list = []
BRCSfinal_list = []
BRTNfinal_list = []
xyfun_list =[]
yfun_list = []
start()
start1()
for n in np.arange(0,len(fl.mhch)):
    mhch = fl.mhch[n]
    strmhch = '%.2g'% mhch
#    mhch = np.array(np.arange(120,131.0,10.0)[n])
    alpmhch = alpmz /(1.0 + (alpmz/(12.0 * PI)) * (33.0 - 2 * nf) * \
                  math.log((mhch**2)/(mz**2))) 

    alpmb= alpmz /(1.0 + (alpmz/(12.0 * PI))*(33.0 - 2 * nf) * \
               math.log(4.89**2/(mz**2)))
    print('alpmb :',alpmb, ' alpmhch:',alpmhch)

# alpha_s at sale 2 GeV.Can read this off from Rev.Part.Properties
    alp2gev = 0.117 /(1.0 + (0.117 /(12.0 * PI)) * (33.0 - 2 * nf) *\
                  math.log(2.0**2 / mz**2))
#####################################################
    alp2gev = 0.310 # Take this value will have same figure as Andrew and Stefano' paper
#    [ PHYSICAL REVIEW D 85, 115002 (2012) ]
#####################################################
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
    print('ms_value',ms,mb)
#mb = 2.95# mass bottom quark
#ms =0.055
    print('MCCCCCCCCCCCCC',mc , mb, ms,fl.eeHH_event())
###################################################################################                                                                                                                                                                                                  
    if mhch > ma:
      lamaw=9.0 * gf**2 * mw**4 * mhch * gaw /(16.0 * PI**3)
    else:
      lamaw=0.00
###########################################################
### LOOP CALCULATION and  X, Y, Z dependence

    z_x = np.array([0.1] * len(x)) # create same length of z array for x array
    z_y = np.array([0.1] * len(y)) # create same length of z array for y array
    for Y in y:

        for X in x:
#for Z in np.arange(0.0,5.1,0.5):
                 Z = 1.0
#          for Z in z:          
#                 print(X,Y,Z,brtn(X,Y,Z),brcs(X,Y,Z),brcb(X,Y,Z),\
#                       brtn(X,Y,Z) + brcs(X,Y,Z) + brcb(X,Y,Z))

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
# Contour of Branching ratio for first H+ > c,s x-axis: |X|, y-axis = |Y|; |Z|= 0.1
    plt.figure()
    ContourBRCS = plt.contour(xarray,yarray,BRCS,levels = linecs,\
    colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
    plt.title('BR($H^{\pm} \\to $ cs)')#plot title
    plt.xlabel('|X|')
    plt.ylabel('|Y|')
    plt.text(12.0, 0.2, r'|Z|= '+ str(Z) +',$M_{H^{\pm}}$= '+ strmhch +' GeV')# text|Z|= fixed MH+= GeV 
    plt.axis([0.0, 30.0, 0.0, 0.5])# plot x,y-axis limit : ([xmin,xmax,ymin,ymax])
    plt.clabel(ContourBRCS, inline= 0.05, fontsize= 10)#plot
    plt.colorbar(ContourBRCS)
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('hcs.png')
    plt.close()
    print('--------------------------------------------------')      
# Contour of Branching ratio for first H+ > c,b x-axis: |X|, y-axis = |Y|; |Z|= 0.1
    plt.figure()
    ContourBRCB = plt.contour(xarray,yarray,BRCB,levels = linecs,\
     colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
    plt.title('BR($H^{\pm} \\to $ cb) ') #plot title
    plt.xlabel('|X|')
    plt.ylabel('|Y|')
    plt.text(12.0, 0.2, r'|Z|= '+ str(Z) +',$M_{H^{\pm}}$= '+ strmhch +' GeV')# text|Z|= fixed MH+= GeV 
    plt.axis([0.0, 30.0, 0.0, 0.5])# plot x,y-axis limit : ([xmin,xmax,ymin,ymax])
    plt.clabel(ContourBRCB, inline= 10, fontsize= 10)#plot
    plt.colorbar(ContourBRCB)
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('hcb.png')
#    plt.show()
    plt.close()
    print('--------------------------------------------------')
# Contour of Branching ratio for H+ > tau,tau-neutrino x-axis: |X|, y-axis = |Y|; |Z|= 0.1
    plt.figure()
    ContourBRTN = plt.contour(xarray,yarray,BRTN,levels = linecs,\
                colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
    plt.title('BR($H^{\pm} \\to $ $\\tau \\nu_\\tau $)')#plot title
    plt.xlabel('|X|')
    plt.ylabel('|Y|')
    plt.text(4.0, 0.2, r'|Z|= '+ str(Z) +',$M_{H^{\pm}}$= '+ strmhch +' GeV')# text|Z|= fixed MH+= GeV 
    plt.axis([0.0, 25.0, 0.0, 0.5])# plot x,y-axis limit : ([xmin,xmax,ymin,ymax])
    plt.clabel(ContourBRTN, inline=0.05, fontsize= 10)#plot 
    plt.colorbar(ContourBRTN)
    plt.grid(axis='y', linestyle='-', color='0.8') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.8') # show x-axis grid line
#   plt.savefig('htaunv.png')
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
    plt.title('BR($t \\to $ $H^{\pm}$b) X BR($H^{\pm} \\to$ cb)')#plot title
    plt.xlabel('|X|')
    plt.ylabel('|Y|')
    plt.text(25.0,0.4, r'|Z|='+ str(Z) +',\n $M_{H^+}$= '+ strmhch +' GeV')
    plt.colorbar(ContourBRTHBBRCB)
    plt.clabel(ContourBRTHBBRCB, inline=0.05, fontsize=10)#plot 
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('M{H^{\pm}}= '+ strmhch +' GeV,thbhcbXY.png')
#    plt.show()
    plt.close()
    print('--------------------------------------------------')
# Contour of Branching ratio for t >H+b - H+ > cs in X,Y plane
    BRTHBBRCS = csevents(xarray,yarray,0.1)
    plt.figure()
    plt.axis([0.0, 40.0, 0.0, 0.5]) 
    plt.title('BR($t \\to $ $H^{\pm}$b) X BR($H^{\pm} \\to $ cs)')#plot title
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
#   plt.savefig('M{H^{\pm}}= '+strmhch +' GeV,thbhcsXY.png')
#    plt.show()
    plt.close()
#    print('--------------------------------------------------')
# Contour of Branching ratio for t >H+b - H+ > tau,nv in Y,Z plane
#    for Y in y:
#          for Z in z: 
#              X = 5.0
#    ynum,znum = np.meshgrid(y,z)
#    BRTHBBRTN = tnevents(X,ynum,znum)
#    plt.figure()
#    plt.axis([0.0, 0.3, 0.0, 5.0]) 
#    ContourBRTHBBRTN = plt.contour(ynum,znum,BRTHBBRTN,\
#                    levels = linecs5,colors = ['black','royalblue','purple','darkgreen','brown','red','black']) 
#    plt.title('BR($t \\to $ $H^{\pm}$b) X BR($H^{\pm} \\to $ $\\tau \\nu_\\tau $)')#plot title
#    plt.plot(np.array([0.22] * len(z)),z,color='red',label = '|XY$^*$| < 1.1', linestyle='dashed',linewidth = 1.0)
#    plt.plot(np.array([0.14] * len(z)),z,color='blue',label = '|XY$^*$| < 0.7', linestyle='dotted',linewidth = 1.0)
#    plt.legend(loc= 'center right',frameon=False,fontsize = 'small')
#    plt.clabel(ContourBRTHBBRTN, fontsize= 8.0,fmt = '%.4f')#plot contour line
#    plt.xlabel('|Y|')
#    plt.ylabel('|Z|')
#    plt.text(0.22, 4.0, r'|X|='+ str(X))
#    plt.text(0.22, 3.0, str(mhch) +'GeV')
#    plt.colorbar(ContourBRTHBBRTN)
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,thbhtaunvYZ.png')
#   plt.show()
#    plt.close()
##############################################################################
# |X| value for theta, tangentbeta, tangamma, delta (only care about first H+ this time)
# Condition d = 1 u = 2 e =3

# create empty list to fill in all results of |X|,|Y|,|Z|
# from theta, tangentbeta, tangamma, delta
    absolutexlist =[]
    absoluteylist =[]
    absolutezlist =[]
    xyfun = []
    yfun = []
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
#def X2(i,j,k,l):
#    return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #X2
#def Y2(i,j,k,l):
#    return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
#def Z2(i,j,k,l) :
#    return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Z2
#############################################################################
    print('___________________________________________')
#print(xyfun)
# (4 parameters),for BRCS result
    BRCSfinal = []

# (4 parameters),for BRCB result
    BRCBfinal = []
# (4 parameters),for BRCB^2 result
    BRCB2final = []
# (4 parameters),for BRTN result
    BRTNfinal = []  
#############################################################################  
# paste here
    start3()
    BRCB_TN = np.array(BRCBfinal) * np.array(BRTNfinal)# product of BRCB and BRTN
    BRCB_CS = np.array(BRCBfinal) * np.array(BRCSfinal)# product of BRCB and BRCS
    BRCBPLUSCS = np.array(BRCBfinal) + np.array(BRCSfinal)# Sum of BRCB and BRCS
#    print('MAXCB + CS',max(BRCBPLUSCS),BRCBfinal[max_valueposition(BRCBPLUSCS)],\
#          BRCSfinal[max_valueposition(BRCBPLUSCS)])
#    print('MAXBRCBCS', len(BRCB_CS),max(BRCB_CS),max_valueposition(BRCB_CS),\
#         BRCBfinal[max_valueposition(BRCB_CS)],BRCSfinal[max_valueposition(BRCB_CS)])#
#    print('MAXBRCBTN',len(BRCB_TN),max(BRCB_TN), max_valueposition(BRCB_TN),\
#         BRCBfinal[max_valueposition(BRCB_TN)],BRTNfinal[max_valueposition(BRCB_TN)] )#
#    print('MAXBRTN',max(BRTNfinal), max_valueposition(BRTNfinal))
#    print('Y',len(absoluteylist),len(absolutexlist))
#4 parameters for t>H+b to H+ >cb result
    BRTHBBRCB1 = cbevents(np.array(absolutexlist),np.array(absoluteylist),np.array(absolutezlist))
#4 parameters for t>H+b to H+ >cs result
    BRTHBBRCS1 = csevents(np.array(absolutexlist),np.array(absoluteylist),np.array(absolutezlist))
#4 parameters for t>H+b to H+ >cs + t>H+b to H+ >cb result
    BRTHQUARK = totalcbcs(np.array(absolutexlist),np.array(absoluteylist),np.array(absolutezlist))
#################################
#PLOT OF CONTOUR WITH RANGE OS 4 PARAMETERS INTO |X|,|Y|,|Z = 0.1| WITH CS,CB,TAUNV
# (4 parameters):A,B, BRCS contour plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
#    plt.figure(1)
#    plt.subplot(221)
    Contourbrcs = plt.contourf(A,B, \
          np.resize(BRCSfinal,len(BRCSfinal)).reshape(len(B),len(A)),\
          colors = ['black','royalblue','purple','darkgreen','brown','red','black'],\
           levels = np.arange(0.15,0.40,0.05))
    plt.axis([min(A), max(A),min(B), max(B)])
#plt.clabel(Contourbrcs, inline= 0.01, fontsize=10)# contour level show
    plt.colorbar(Contourbrcs)
    plt.title('BR($H^{\pm} \\to $ cs), $M_{H^{\pm}}$: '+ strmhch +' GeV ')#,$M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,cs.png')
    plt.show()
    plt.close()
#    print('--------------------------------------------------')
# (4 parameters):A,B, BRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
#    plt.subplot(222)
    Contourbrcb = plt.contourf(A,B, \
           np.resize(BRCBfinal,len(BRCBfinal)).reshape(len(B),len(A)),\
           colors = ['black','royalblue','purple','darkgreen','brown','red','black'],\
           levels = np.arange(0.0,1.2,0.2))
#    plt.clabel(Contourbrcb, inline= 0.02, fontsize= 9)# contour level show
    plt.colorbar(Contourbrcb)
    plt.title('BR($H^{\pm} \\to $ cb), $M_{H^{\pm}}$: '+ strmhch +' GeV ')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,cb.png')
    plt.show()
    plt.close()
#    print('--------------------------------------------------')
#(4 parameters): A,B, BRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
#    plt.subplot(223)
    Contourbrtn = plt.contourf(A,B, \
           np.resize(BRTNfinal,len(BRTNfinal)).reshape(len(B),len(A)),\
           colors = ['black','royalblue','purple','darkgreen','brown','red','black'],\
             levels = np.arange(0.0,0.8,0.1))
#plt.clabel(Contourbrtn)# contour level show
    plt.colorbar(Contourbrtn)
    plt.title('BR($H^{\pm} \\to $ $\\tau \\nu_\\tau $), $M_{H^{\pm}}$: '+ strmhch +' GeV ')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,TN.png')
    plt.show()
    plt.close()
#    print('--------------------------------------------------')
#(4 parameters): A,B, BRCB + BRCS contour plot 
#    plt.subplot(224)
    ContourBRCBCS = plt.contourf(A,B, \
           np.resize(BRCBPLUSCS ,len(BRCBPLUSCS )).reshape(len(B),len(A)),\
           colors = ['black','royalblue','purple','yellow','brown','red','gray','green'])#, levels = np.arange(0.06, 0.16,0.02))
#    plt.clabel(ContourBRCBCS, inline= 0.02, fontsize= 9)# contour level show
    plt.colorbar(ContourBRCBCS)
    plt.title(' BR($H^{\pm} \\to $ hadrons), $M_{H^{\pm}}$: '+ strmhch +' GeV ')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.savefig('branchingratio'+ strmhch +'.png')
#    plt.subplots_adjust(top=0.99, bottom=0.06, left=0.1, right=0.99, hspace= 0.6,
#                    wspace=0.3)
#    plt.show()
#    plt.close()
#(4 parameters): A,B, BRTHBBRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
    ContourBRTHBBRCB1 = plt.contourf(A,B, \
           np.resize(BRTHBBRCB1,len(BRTHBBRCB1)).reshape(len(B),len(A)),\
           levels = linecs2,colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
    plt.title('BR($t \\to $ $H^{\pm}$b) X BR($H^{\pm} \\to $ cb),\n $M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
    plt.clabel(ContourBRTHBBRCB1, inline= 0.01, fontsize=10)# contour level show
    plt.colorbar(ContourBRTHBBRCB1)
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,thbhcb.png')
    plt.close()
    print('--------------------------------------------------')
#(4 parameters): A,B, BRTHBBRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
    plt.title('BR($t \\to $ $H^{\pm}$b) X BR($H^{\pm} \\to $ cs),\n $M_{H^{\pm}}$= '+ strmhch +' GeV')#plot title
    ContourBRTHBBRCS1 = plt.contourf(A,B, \
           np.resize(BRTHBBRCS1,len(BRTHBBRCS1)).reshape(len(B),len(A)),\
           levels = linecs2,colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
    plt.clabel(ContourBRTHBBRCS1, inline= 0.01, fontsize=10)# contour level show
    plt.colorbar(ContourBRTHBBRCS1)
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('M{H^{\pm}}= '+ str(strmhch) +' GeV,thbhcs.png')
    plt.close()
    print('--------------------------------------------------')
#(4 parameters): A,B, BRTHBBRCB + BRTHBBRCS contour plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
    plt.title('BR($t \\to $ $H^{\pm}$b) X [BR($H^{\pm} \\to $ cb)\
+ BR($H^{\pm} \\to $ cs)],\n $M_{H^{\pm}}}$= '+strmhch +' GeV')#plot title
    ContourTOTAL = plt.contourf(A,B, \
           np.resize(BRTHQUARK,len(BRTHQUARK)).reshape(len(B),len(A)),\
           levels = linecs2,colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.clabel(ContourTOTAL, inline= 0.01, fontsize=10)# contour level show
    plt.colorbar(ContourTOTAL)
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#   plt.savefig('M{H^{\pm}}= '+ strmhch +' GeV,thbtotal.png')
    plt.close()
    print('--------------------------------------------------')
#(4 parameters): A,B, |X| plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
    plt.title('|X| in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane')#plot title
    plt.axis([0.0, max(A),min(B), max(B)])
    contourx = plt.contour(A,B,\
            np.resize(np.array(absolutexlist),len(np.array(absolutexlist))).reshape(len(B),len(A)))
    plt.clabel(contourx, inline= 0.01, fontsize=10)# contour level show
    plt.colorbar(contourx)
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('x.png')
#   plt.show()
    plt.close()
    print('--------------------------------------------------')
#(4 parameters): A,B, REAL(XY^{*}) plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
    plt.title('REAL(XY^*) in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane,$M_{H^+}$= '+ str(strmhch) +' GeV')#plot title
    contourxy = plt.contour(A,B,\
            np.resize(np.array(np.real(xyfun)),len(np.array(np.real(xyfun)))).reshape(len(B),len(A)),\
            levels = np.arange(-1.1,0.9,0.2),cmap = 'brg')#,levels = linecs3)#
    plt.clabel(contourxy, inline= 0.01, fontsize=10)# contour level show
#   plt.colorbar(contourxy)
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.show()
#   plt.savefig('M{H^{\pm}}= '+ strmhch +' GeV,realxy.png')
    plt.close()
    print('--------------------------------------------------')
#(4 parameters): A,B, IM(XY^{*}) plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
    plt.title('IM(XY^*) in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane,$M_{H^+}$= '+ strmhch +' GeV')#plot title
    contourxy2 = plt.contour(A,B,np.resize(np.array(np.imag(xyfun)),len(np.array(np.imag(xyfun)))).reshape(len(B),len(A)),\
                            levels = linecs4,cmap = 'brg')
#   plt.clabel(contourxy2, inline= 0.01, fontsize=10)# contour level show
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.show()
#   plt.savefig('M{H^{\pm}}= '+ strmhch +' GeV,imxy.png')
    plt.close()
    print('--------------------------------------------------')
#(4 parameters): A,B, |Y| plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
    plt.title('|Y| in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane')#plot title
    contoury = plt.contour(A,B,\
            np.resize(np.array(absoluteylist),len(np.array(absoluteylist))).reshape(len(B),len(A)))
    plt.clabel(contoury, inline= 0.01, fontsize=10)# contour level show
    plt.colorbar(contoury)
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('y.png')
#   plt.show()
    plt.close()
    print('--------------------------------------------------')
#(4 parameters): A,B, |Z| plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
    plt.title('|Z| in '+ readlist[int(read1)] +','+ readlist[int(read2)] +' plane')#plot title
    contourz = plt.contour(A,B,\
            np.resize(np.array(absolutezlist),len(np.array(absolutezlist))).reshape(len(B),len(A)))
    plt.clabel(contourz, inline= 0.01, fontsize=10)# contour level show
    plt.colorbar(contourz)
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('z.png')
#   plt.show()
    plt.close()
    print('--------------------------------------------------')
#    print('CCCc',fl.eeHH_event()[n] * 0.8**2 * 0.7**2 /np.sqrt(fl.backgroundtagging()* 0.1**2))
#    print(len(fl.sig(fl.eeHH_event() * np.array(BRCBfinal)**2 * 0.7**2,fl.backgroundtagging()* 0.1**2)),fl.sig(fl.eeHH_event() * np.array(BRCBfinal)**2 * 0.7**2,fl.backgroundtagging()* 0.1**2))
    print('___________________________________________________')
#    print((fl.eeHHcbcb_2bsignal(BRCBfinal,BRCBfinal) + fl.eeHHcbcb_1bsignal(BRCBfinal,BRCBfinal) + \
#                   fl.eeHHcbcs_2bsignal(BRCBfinal,BRCSfinal) + fl.eeHHcbcs_1bsignal(BRCBfinal,BRCSfinal) + \
#                   fl.eeHHcbcb_0bsignal(BRCBfinal,BRCBfinal) + fl.eeHHcbcs_0bsignal(BRCBfinal,BRCSfinal) + fl.eeHHcscs_0bsignal(BRCSfinal,BRCSfinal)) * fl.epsilon)
    print('---------------------------------------------------')
########################################################################
#### 2 jet tagged plot array
    twotagsig = fl.eeHH_event()[n] * (fl.eeHHcbtn_1bsignal(BRCBfinal,BRTNfinal) +\
                fl.eeHHcbtn_0bsignal(BRCBfinal,BRTNfinal) + \
                fl.eeHHcstn_1bsignal(BRCSfinal,BRTNfinal) + fl.eeHHcstn_0bsignal(BRCSfinal,BRTNfinal) \
                ) * fl.selection_2j / np.sqrt(fl.backgroundtagging2())
#### 4 jet 2b-tagged plot array
    fourtag2bsig = fl.eeHH_event()[n] * (fl.eeHHcbcb_2bsignal(BRCBfinal,BRCBfinal) +\
                   fl.eeHHcbcb_1bsignal(BRCBfinal,BRCBfinal) + \
                   fl.eeHHcbcs_2bsignal(BRCBfinal,BRCSfinal) + fl.eeHHcbcs_1bsignal(BRCBfinal,BRCSfinal) + \
                   fl.eeHHcbcb_0bsignal(BRCBfinal,BRCBfinal) + fl.eeHHcbcs_0bsignal(BRCBfinal,BRCSfinal) + \
                   fl.eeHHcscs_0bsignal(BRCSfinal,BRCSfinal)) \
                              * fl.epsilon / np.sqrt(fl.backgroundtagging())
#### 4 jet 1b-tagged plot array
    fourtag1bsig = fl.eeHH_event()[n] * (fl.real_b_cbcb(BRCBfinal,BRCBfinal) +\
                   fl.fake_b_cbcb(BRCBfinal,BRCBfinal) + fl.real_b_cbcs(BRCBfinal,BRCSfinal) + \
                   fl.fake_b_cbcs(BRCBfinal,BRCSfinal) + fl.real_b_cscs(BRCSfinal,BRCSfinal) + \
                   fl.fake_b_cscs(BRCSfinal,BRCSfinal)) \
                              * fl.epsilon / np.sqrt(fl.backgroundtagging())
#########################################################################
#    plt.figure()
#    Contoursignal0 = plt.contourf(A,B, \
#         np.resize(fl.sig(fl.eeHH_event()[n] * np.array(BRCBfinal)**2 * 0.7**2 * fl.epsilon , 1.0/np.sqrt(fl.backgroundtagging()) ),\
#               len(fl.sig(fl.eeHH_event()[n] * np.array(BRCBfinal)**2 * 0.7**2 * fl.epsilon , 1.0/np.sqrt(fl.backgroundtagging()) ))).\
#        reshape(len(B),len(A)),colors = ['black','royalblue','purple','orange','brown','red','gray'])# ,levels = np.arange(1.0,6.0,1.0))
#    plt.title('Significance of H_char>cbcb in '+\
#             readlist[int(read1)] +','+ readlist[int(read2)])
#    plt.xlabel(readlist[int(read1)])# x-axis label
#    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.colorbar(Contoursignal0)
#    plt.savefig('significance_cbcb'+ str(mhch) +'.png')
#    plt.show()
#    plt.close()
#############################################
#############################################
# (4 parameters):A,B, BRCB^2 contour plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure()
    Contourbrcb_2 = plt.contourf(A,B, \
           np.resize(np.array(BRCBfinal)**2,\
                 len(np.array(BRCBfinal)**2) ).reshape(len(B),len(A)),\
           colors = ['black','royalblue','purple','yellow','brown','red','gray','green'],\
           levels = np.arange(0.0, 0.8,0.1))
#    plt.clabel(Contourbrcb_2)# contour level show
    plt.colorbar(Contourbrcb_2)
    plt.title('BR($H^{\pm} \\to $ cb)^2,$M_{H^{\pm}}$= '+ strmhch +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('CharHcbsquared'+ strmhch +'.png')
    plt.show()
    plt.close()
#######################################################
# (4 parameters):A,B,plot (BR(H+> cb + cs)  * BR(H+> tn) )
    plt.figure()
    Contourbrcbtn = plt.contourf(A,B, \
           np.resize(BRCBPLUSCS * BRTNfinal,len(BRCBPLUSCS * BRTNfinal)).reshape(len(B),len(A)),\
           colors = ['black','royalblue','purple','yellow','brown','red','gray','green'])#, levels = np.arange(0.06, 0.16,0.02))
#    plt.clabel(Contourbrcbtn, inline= 0.02, fontsize= 9)# contour level show
    plt.colorbar(Contourbrcbtn)
    plt.title(' BR($H^{\pm} \\to $ hadrons) X BR($H^{\pm} \\to $ $\\tau \\nu_\\tau $) , $M_{H^{\pm}}$: '+ strmhch +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('CharHcbpcsmtn'+ strmhch +'.png')
    plt.show()
    plt.close()
#######################################################
# (4 parameters):A,B,plot (BR(H+> cb)  * BR(H+> cs) )
#    plt.figure()
#    Contourbrcbcs = plt.contourf(A,B, \
#           np.resize(BRCBPLUSCS,len(BRCBPLUSCS)).reshape(len(B),len(A)),\
#           colors = ['black','royalblue','purple','yellow','brown','red','gray','green'])#, levels = np.arange(0.06, 0.16,0.02))
#    plt.clabel(Contourbrcbcs, inline= 0.02, fontsize= 9)# contour level show
#    plt.colorbar(Contourbrcbcs)
#    plt.title(' BR($H^{\pm} \\to $ cb) * BR($H^{\pm} \\to $ cs) , $M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
#    plt.xlabel(readlist[int(read1)])# x-axis label
#    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.savefig('CharHcbcs'+ str(mhch) +'.png')
#    plt.show()
#    plt.close()
###############################################
# (4 parameters):A,B, 4jet tagged plots
    plt.figure()
    signal4jet = plt.contourf(A,B, \
        np.resize( fourtag2bsig ,\
              len( fourtag2bsig)).\
        reshape(len(B),len(A)),\
        colors = ['black','royalblue','purple','yellow','brown','red','gray','green'] ,levels = np.arange(0.0,4.5,0.5))
    plt.title('S/$\sqrt{B}$ 4jet2b '+\
             ', $M_{H^{\pm}}$= '+ strmhch +' GeV')
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.colorbar(signal4jet)
    plt.savefig('sig_4jet2b'+ strmhch +'.png')
    plt.show()
    plt.close()
###############################################
# (4 parameters):A,B, 4jet tagged plots 1-b
    plt.figure()
    signaloneb4jet = plt.contourf(A,B, \
        np.resize( fourtag1bsig ,\
              len( fourtag1bsig)).\
        reshape(len(B),len(A)),\
        colors = ['black','royalblue','purple','yellow','brown','red','gray','green'] ,levels = np.arange(0.0,3.0,0.5))
    plt.title('S/$\sqrt{B}$ 4jet1b '+\
             ', $M_{H^{\pm}}$= '+ strmhch +' GeV')
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.colorbar(signaloneb4jet)
    plt.savefig('sig_4jet1b'+ strmhch +'.png')
    plt.show()
    plt.close()
################################
#(4 parameters):A,B,4jet notagging plot (BR($H^{\pm} \\to $ cb + cs))
    plt.figure()
    signalnotag4jet = plt.contourf(A,B, \
        np.resize(fl.eeHH_event()[n] * (BRCBPLUSCS**2)\
                  * fl.epsilon / np.sqrt(fl.backgroundnotagging()) ,\
              len(fl.eeHH_event()[n] * (BRCBPLUSCS**2)\
                  * fl.epsilon / np.sqrt(fl.backgroundnotagging()))).\
        reshape(len(B),len(A)),colors = ['black','royalblue','purple','yellow','brown','red','gray','green'],\
        levels = np.arange(0.0,2.5,0.5))
    plt.colorbar(signalnotag4jet)
    plt.title('S/$\sqrt{B}$ 4jet0b,$M_{H^{\pm}}$= '+ strmhch +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('sig_4jet0b'+ strmhch +'.png')
    plt.show()
    plt.close()
####################################
#    plt.figure()
#    Contourbrcstn = plt.contour(A,B, \
#        np.resize(np.array(BRCSfinal) * np.array(BRTNfinal),\
#              len(np.array(BRCSfinal) * np.array(BRTNfinal))).reshape(len(B),len(A)),\
#        cmap = 'brg')#,levels = np.arange(0.0,5.0,1.0))
#    plt.colorbar(Contourbrcstn)
#    plt.title(' BR($H^{\pm} \\to $ cs) * BR($H^{\pm} \\to $ tn),$M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
#    plt.xlabel(readlist[int(read1)])# x-axis label
#    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.savefig('CharHcstn'+ str(mhch) +'.png')
#    plt.close()
###################################
#    BR($H^{\pm} \\to $ cb * $\\tau \\nu_\\tau$) with tagging efficiencies
    plt.figure()
    Contoursignal7 = plt.contourf(A,B, \
        np.resize(twotagsig ,\
              len(twotagsig )).reshape(len(B),len(A)),\
        colors = ['black','royalblue','purple','yellow','brown','red','gray','green'],levels = np.arange(0.0,3.0,0.5))
    plt.colorbar(Contoursignal7)
    plt.title('S/$\sqrt{B}$ 2jet1b,$M_{H^{\pm}}$= '+ strmhch +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('sig_2jet1b'+ strmhch +'.png')
    plt.show()
    plt.close()
#    BR($H^{\pm} \\to  (cb+cs) * \\tau \\nu_\\tau $) with no-tagging efficiencies
    plt.figure()
    Contoursignal8 = plt.contourf(A,B, \
        np.resize(fl.eeHH_event()[n] * 2.0 * BRCBPLUSCS * np.array(BRTNfinal) / np.sqrt(fl.backgroundnotagging2()) ,\
              len(fl.eeHH_event()[n] * 2.0 * BRCBPLUSCS * np.array(BRTNfinal) / np.sqrt(fl.backgroundnotagging2()))).\
                  reshape(len(B),len(A)),colors = ['black','royalblue','purple','orange','brown','red','gray','green'],\
                  levels = np.arange(0.0,2.5,0.5))
    plt.colorbar(Contoursignal8)
    plt.title('S/$\sqrt{B}$ 2jet0b,$M_{H^{\pm}}$= '+ strmhch +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('sig_2jet0b'+ strmhch +'.png')
    plt.show()
    plt.close()
####################################
#    print(len(twotagsig),A)
#    np.savetxt('fourtag2bsig.txt',(fourtag2bsig),header="fourtag2bsig",fmt = '%10.5f')
#    np.savetxt('twotagsig.txt',(twotagsig),header="twotagsig",fmt = '%10.5f')
####################################
def start4():#Scan-plot
        start1()
        for n in np.arange(0,len(fl.mhch)):
            mhch = fl.mhch[n]
            strmhch = '%.2g'% mhch
            alpmhch = alpmz /(1.0 + (alpmz/(12.0 * PI)) * (33.0 - 2 * nf) * \
                  math.log((mhch**2)/(mz**2))) 

            alpmb= alpmz /(1.0 + (alpmz/(12.0 * PI))*(33.0 - 2 * nf) * \
               math.log(4.89**2/(mz**2)))
#            print('alpmb :',alpmb, ' alpmhch:',alpmhch)

# alpha_s at sale 2 GeV.Can read this off from Rev.Part.Properties
            alp2gev = 0.117 /(1.0 + (0.117 /(12.0 * PI)) * (33.0 - 2 * nf) *\
                  math.log(2.0**2 / mz**2))
#####################################################
            alp2gev = 0.310 # Take this value will have same figure as Andrew and Stefano' paper
#    [ PHYSICAL REVIEW D 85, 115002 (2012) ]
#####################################################
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
#            print('ms_value',ms,mb)
#mb = 2.95# mass bottom quark
#ms =0.055
            print('MC MB MS',mc , mb, ms,fl.eeHH_event())
###################################################################################                                                                                                                                                                                                  
            if mhch > ma:
                lamaw=9.0 * gf**2 * mw**4 * mhch * gaw /(16.0 * PI**3)
            else:
                lamaw=0.00
            BRCBfinal = []
            BRCSfinal = []
            BRTNfinal = []
            list_a =[]
            list_b =[]
            list_c =[]
            for var_i in the:
                for var_j in tbe:
                    for var_k in tga:
#            for var_j,var_k in (tbe,tga):
                        my_tuple = (var_i,var_j, var_k,0.0) 
                        list_a.append(var_i)
                        list_b.append(var_j)
                        list_c.append(var_k)
                        
#                    for counter in range(4):
#                            my_tuple += (reference_array[counter],)
#                    absolutexlist.append(abs(X2(*my_tuple))) # unpacked bracket of my_tuple by * sign
#                    absoluteylist.append(abs(Y2(*my_tuple)))
#                    absolutezlist.append(abs(Z2(*my_tuple)))
#                    xyfun.append(complexyfunction(*my_tuple))
                        BRCBfinal.append(brcb(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple))))
                        BRCSfinal.append(brcs(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple))))
                        BRTNfinal.append(brtn(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple)))) 
                        arrayhadron = np.array(BRCBfinal) + np.array(BRCSfinal)
                        arrayhadrontn =  arrayhadron * np.array(BRTNfinal)
                        twotagsig = fl.eeHH_event()[n] * (fl.eeHHcbtn_1bsignal(BRCBfinal,BRTNfinal) +\
                fl.eeHHcbtn_0bsignal(BRCBfinal,BRTNfinal) + \
                fl.eeHHcstn_1bsignal(BRCBfinal,BRTNfinal) + fl.eeHHcstn_0bsignal(BRCBfinal,BRTNfinal) \
                ) * fl.selection_2j  / np.sqrt(fl.backgroundtagging2())
                        fourtag2bsig = fl.eeHH_event()[n] * (fl.eeHHcbcb_2bsignal(BRCBfinal,BRCBfinal) +\
                   fl.eeHHcbcb_1bsignal(BRCBfinal,BRCBfinal) + \
                   fl.eeHHcbcs_2bsignal(BRCBfinal,BRCSfinal) + fl.eeHHcbcs_1bsignal(BRCBfinal,BRCSfinal) + \
                   fl.eeHHcbcb_0bsignal(BRCBfinal,BRCBfinal) + fl.eeHHcbcs_0bsignal(BRCBfinal,BRCSfinal) + \
                   fl.eeHHcscs_0bsignal(BRCSfinal,BRCSfinal)) \
                   * fl.epsilon / np.sqrt(fl.backgroundtagging())
                        fourtag1bsig = fl.eeHH_event()[n] * (fl.real_b_cbcb(BRCBfinal,BRCBfinal) +\
                   fl.fake_b_cbcb(BRCBfinal,BRCBfinal) + fl.real_b_cbcs(BRCBfinal,BRCSfinal) + \
                   fl.fake_b_cbcs(BRCBfinal,BRCSfinal) + fl.real_b_cscs(BRCSfinal,BRCSfinal) + \
                   fl.fake_b_cscs(BRCSfinal,BRCSfinal))  * \
                   fl.epsilon / np.sqrt(fl.backgroundtagging())\
                   * fl.epsilon / np.sqrt(fl.backgroundtagging())
                        fourjnontag = fl.eeHH_event()[n] *  (arrayhadron**2)\
                  * fl.epsilon / np.sqrt(fl.backgroundnotagging() )
                        twojnontag = fl.eeHH_event()[n] * 2.0 * arrayhadrontn  / np.sqrt(fl.backgroundnotagging2())
            d1 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           'BRCB': [*BRCBfinal]})
            d2 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           'BRCS': [*BRCSfinal]})
            d3 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           'BRTN': [*BRTNfinal]})
            d4 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           'BRCB+CS': [*fourjnontag]})    
            d5 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           'BRCBTN': [*twojnontag]})      
            d6 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           '2tagsig': [*twotagsig]}) 
            d7 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           '4tag2bsig': [*fourtag2bsig]}) 
            d8 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           '4tag1bsig': [*fourtag1bsig]}) 
            d9 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           '4tag1bsig': [*arrayhadrontn]})  
            d10 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'tangamma': [*list_c],
                           '4tag1bsig': [*arrayhadron]})  

            d1.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/dcb'+ strmhch +'.csv',\
                   index = None, header=True) 
            d2.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/dcs'+ strmhch +'.csv',\
                   index = None, header=True)
            d3.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/dtn'+ strmhch+'.csv',\
                   index = None, header=True)
            d4.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/dcbpcs'+ strmhch +'.csv',\
                   index = None, header=True)
            d5.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/dcbmtn'+ strmhch +'.csv',\
                   index = None, header=True)
            d6.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/d2jtag'+ strmhch +'.csv',\
                   index = None, header=True)
            d7.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/d4j2btag'+ strmhch +'.csv',\
                   index = None, header=True)
            d8.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/d4j1btag'+ strmhch +'.csv',\
                   index = None, header=True)
            d9.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/dhadrontn'+ strmhch +'.csv',\
                   index = None, header=True)
            d10.to_csv (r'/Volumes/Backup/Allen/PHDmeeting/3HDM codes/significanceplot/dhadron'+ strmhch +'.csv',\
                   index = None, header=True)
            print('alpmb :',alpmb, ' alpmhch:',alpmhch)
        print('mhchfinal', mhch)
#        xtbe,ytga = np.meshgrid(tbe,tga)
        #zzz = np.reshape(np.array(BRCBfinal),( len(the) , len(tbe) * len(tga) ))
        #break whole list to [ len(tga) * len(tbe) ], * len(the) 
#        list_bx = np.split(np.array(list_b),len(the)) 
#        list_cx = np.split(np.array(list_c),len(the))
#        cbcb = np.split(np.array(BRCBfinal), len(the))
#        cscs = np.split(np.array(BRCSfinal), len(the))
#        tntn = np.split(np.array(BRTNfinal), len(the))
#        print(cbcb,type(cbcb),len(cbcb),np.array(cbcb).ndim)
#        sum_cbcb = []
#        sum_bx = []
#        for i in np.arange(0,len(the)):
            # break single [ len(tga) * len(tbe) ] to [len(tga), len(tbe)]
#            list_bx[i] = np.split(np.array(list_bx[i]),len(tbe) )
#            list_cx[i] = np.split(np.array(list_cx[i]),len(tbe) )
#            cbcb[i] = np.split(np.array(cbcb[i]),len(tbe))
#            cscs[i] = np.split(np.array(cscs[i]),len(tbe))
#            tntn[i] = np.split(np.array(tntn[i]),len(tbe))
#            print(len(cbcb[i]),len(np.array(tntn[i])))
#            sum_cbcb = sum_cbcb + cbcb[i]
#            print('-',len(np.array(sum_cbcb).flatten()))
#            fig = plt.figure()
#            ax = Axes3D(fig)
#            surf = ax.plot_surface(xtbe,ytga,np.array(tntn[i]) ,\
#                rstride=1, cstride=1, cmap=cm.hsv,linewidth=0, antialiased=True)
#            surf1 = ax.plot_surface(xtbe,ytga,np.array(cbcb[i]) ,\
#                 rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=True)
#            surf2 = ax.plot_surface(xtbe,ytga,np.array(cscs[i]) ,\
#                 rstride=1, cstride=1, cmap=cm.cool,linewidth=0, antialiased=True)
#            fig.colorbar(surf, shrink= 0.5, aspect=5)
#            fig.colorbar(surf1, shrink=0.5, aspect=5)
#            fig.colorbar(surf2, shrink=0.5, aspect=5)
#        ax.plot_trisurf(d1['tangamma'],d1['tanbeta'],d1['BRCB'],cmap=cm.jet, linewidth=0.2)
#        ax.plot_trisurf(d1['tangamma'],d1['tanbeta'],d1['BRCS'],cmap=cm.jet, linewidth=0.2)
#        ax.plot_trisurf(d1['tangamma'],d1['tanbeta'],d1['BRTN'],cmap=cm.jet, linewidth=0.2)
        
#            plt.show()
#        cbcb_3600 = np.split(np.array(BRCBfinal), len(tga) * len(tbe))
#        print('--',len(cbcb_3600))
#        print(len(np.array(sum_cbcb).reshape(len(tbe)* len(tga),len(the) )) )
#        sum_cbcb1 = np.array(sum_cbcb).reshape(len(tbe)* len(tga),len(the) )
#        print(len(sum_cbcb1),len(xxx),np.array(xxx).ndim, sum_cbcb1.ndim )
#        fig = plt.figure()
#        ax = Axes3D(fig)
#        surfsum = ax.plot_surface(xxx, yyy,np.array(cbcb_3600) ,\
#                rstride=1, cstride=1, cmap=cm.hsv,linewidth=0, antialiased=True)
#        fig.colorbar(surfsum, shrink= 0.5, aspect=5)
#        plt.show()

#start4()
#####################################################
#print('len',len(np.array(BRCBfinal_list).flatten()) ,np.array(BRCBfinal_list).flatten().ndim)
#plt.figure
#plt.plot(fl.mhch,np.array(mhch_sig4jtag(list(np.array(BRCBfinal_list).flatten()),\
#    list(np.array(BRCSfinal_list).flatten()) )).reshape((len(fl.mhch),\
#         len(np.array(mhch_sig4jtag(list(np.array(BRCBfinal_list).flatten()),\
#    list(np.array(BRCSfinal_list).flatten()) )) )) )  )
#plt.show()
#plt.close