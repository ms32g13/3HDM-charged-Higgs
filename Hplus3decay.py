#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 17:44:28 2018

@author: apple
"""
###     Program which calculates the branching ratios of the charged              
###     Higgs to the channels cs, tn, cb, AW*,hW*. Valid for the 2HDM and              
###     MHDM. It also calculate BR(t -> H+b) x BR(H+ -> cb.. 
###     This mainly focus on AW* and hW*  decay
###

from invariant import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from invariant import *
import vegas
import gvar as gv
import random
from scipy.integrate import dblquad
from scipy import integrate
#####################################
#some constants
ma = 51
mhch = 120
mw = 80.33
mh = 125
PI = np.pi
gf = 1.1663787e-05#0.000011663787 # fermi constant   
def kx(x,y):#a = mA, c = mM+
    return (x/y)**2
ka = kx(ma,mhch)   #Ma^2 / Mhch^2 = k_i or k_j
kw = kx(mw,mhch)    #Mw^2 / Mhch^2 = k_i or k_j
kh = kx(mh,mhch)   #Mh^2 / Mhch^2 = k_i or k_j
print('ka,kw,kh',ka,kw,kh)
pwwsquared = 4.326 #  partial width of w boson ^2 
def gam_2(b):
    return pwwsquared / b**2 # (small gamma)_j 
print('gam_2',gam_2(mhch))
############################################################
##################################### Analytical formula
#   Function of G_ij
def g(i,j): #G_ij function (ex. i for ka, j for kw)
    lam_ij = - 1.0 + 2 * i + 2 * j - (i - j)**2   # lamda_ij
    numerator = j * (1 - j + i) - lam_ij
    denomator = (1 - i) * np.sqrt(lam_ij)
    chunk_1 = PI / 2 + np.arctan(numerator/denomator)    
    chunk_2 = 1 / 3 * (1.0 - i) * (5 * (1.0 + i) - 4 * j - 2 * lam_ij / j)
#    print('ka,kw,lamij',i,j,lam_ij)
#    print('numerator',numerator)
#    print('denomator',denomator)
#    print('chunk_1', chunk_1, 'chunk2',chunk_2) 
    return 1 / 4 * (2 * (j - i - 1.0) * np.sqrt(lam_ij) * chunk_1 + \
                    (lam_ij - 2 * i) * np.log(i) + chunk_2) 

#######################################
def lam3decay(a,b):# Partial width with g(i,j) a= ma, b=mhch
    return 9.0 * gf**2 * mw**4 * b * g(kx(a,b),kx(mw,b)) /(16.0 * PI**3)
#print('analytical',ma,mhch,g(kx(ma,mhch),kx(mw,mhch)),lam3decay(ma,mhch))
############################################################
############################################################
#Boundary Interval
# For H+>AW approach
#x2_min = 0.0
#x2_max = 1 - kx(ma,mhch)
#x2 = random.uniform(x2_min, x2_max)
#print('x2',x2)
#x1_min = 1 - x2 - kx(ma,mhch)
#x1_max = 1 - kx(ma,mhch) / (1 - x2)
############################
# For H+>hW approach
#x1_min = 1 - x2 - kh
#x1_max = 1 - kh / (1 - x2)
#x2_min = 0.0
#x2_max = 1 - kh
#print(x1_min,x1_max,x2_min,x2_max)
#############################
###################################################
def lamaw(a,b):# Partial width with integration method a = ma , b = mhch 

#Intergrand
    def integrand(x,y):#Faw for the formula 59 in https://arxiv.org/pdf/hep-ph/9511342.pdf
# For H+>AW approach   
#        x2_min = 0.0
#        x2_max = 1 - kx(ma,mhch)
#        x2 = (x2_max-x2_min) * x[1]
#        x1_min = 1 - x2 - kx(ma,mhch)
#        x1_max = 1 - kx(ma,mhch) / (1 - x2)
#        x1 = x1_min +(x1_max-x1_min) * x[0]
#        jac = (x1_max-x1_min)*(x2_max-x2_min)
#        print(x[0],x[1])
        numer = (1.0 - x ) * (1.0 - y ) - kx(a,b)
        denom = (1.0 - x - y - kx(a,b) + kx(mw,b))**2 \
                + kx(mw,b) * gam_2(b)  
        return  numer/denom #* jac
    def bounds_y():
        return [0, 1 - kx(a,b) ]
    def bounds_x(y):
        return [1 - y - kx(a,b), 1 - kx(a,b) / (1 - y)]
    FHpaw = integrate.nquad(integrand, [bounds_x, bounds_y])[0]
#            vectorize_Faw = np.vectorize(FHpaw)# Vectorize integrand
#FHpaw = integ( Faw,nitn=100, neval=1000) 
#integral = 0.0
#for x, wgt in integ.random():
#    integral += wgt * Faw(x)
#print(integral)
    
    return 9.0 * gf**2 * mw**4 * b * FHpaw / (16 * PI**3)

print('lamaw()', lamaw(25,120))
for i in np.arange(45 ,91 ):
    for j in np.arange(75,125):
        print('123',i,j,lamaw(i,j))
############################################################
########################################PLOT
           
def ma_mhch_lam():
#    ma_list = []
#    mhch_list = []
    ax = np.arange(10 ,117 )
    hchy = np.arange(120,121)
    xarray,yarray = np.meshgrid(ax,hchy) 
    plt.figure()
    result = plt.contourf(xarray,yarray,lam3decay(xarray,yarray),levels = np.array([0,0.001,0.002,0.008]))
    plt.colorbar(result)
    plt.xlabel('$M_A$')
    plt.ylabel('$M_{H^{\pm}}$')
    plt.show()
    plt.close()
    return
#ma_mhch_lam() # MA_MHCH_LamdaH+>AW*
##########################################################
def lam_mhch():
    lam_aw_mhch = []
    lam_aw_mhch50 = []
    lam_aw_ma = []
    mamass = 45
    chmass = 90
    A = np.arange(25 ,91 ) # Pseudoscalar A mass range
    Hplus = np.arange(75,200
                      )#HPlus H+ mass range
    
    for m in Hplus:
            lam_aw_mhch.append(lamaw( mamass,m))
    for i in A:
            lam_aw_ma.append(lamaw(i,chmass))
            
    plt.figure()
#    plt.plot(Hplus,lam_aw_mhch)
    plt.plot(Hplus,lam3decay( mamass,Hplus))
    plt.plot(Hplus,lam_aw_mhch)
    plt.title('$\Gamma \\to AW^*$ against $M_{\pm}$ with $M_A$ = ' + str(mamass) +' GeV ')
    plt.legend(('Analytical', 'Integration'),\
               loc='upper left', shadow=True,prop={'size': 7.5} )  
    plt.xlabel('$M_{\pm}$')
    plt.ylabel('$\Gamma \\to AW^*$')
    plt.show()
    plt.close()
    
    plt.figure()
    plt.plot(A,lam3decay( A,chmass))
    plt.plot(A,lam_aw_ma)
    plt.title('$\Gamma \\to AW^*$ against $M_{A}$ with $M_{\pm}$ = ' + str(chmass) +' GeV ')
    plt.legend(('Analytical', 'Integration'),\
               loc='upper left', shadow=True,prop={'size': 7.5} )  
    plt.xlabel('$M_{A}$')
    plt.ylabel('$\Gamma \\to AW^*$')
    plt.show()
    plt.close()
    return 

#lam_mhch()# Plot lamaw against mhch        
            