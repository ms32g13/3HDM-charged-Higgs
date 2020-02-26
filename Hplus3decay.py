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
from matplotlib.ticker import MultipleLocator
from invariant import *
import vegas
import gvar as gv
import random
#####################################
#some constants
ma = 85
mhch = 125
print('mhch,mh,ma,mw',mhch,mh,ma,mw)
def kx(x,y):#a = mA, c = mM+
    return (x/y)**2
ka = kx(ma,mhch)   #Ma^2 / Mhch^2 = k_i or k_j
kw = kx(mw,mhch)    #Mw^2 / Mhch^2 = k_i or k_j
kh = kx(mh,mhch)   #Mh^2 / Mhch^2 = k_i or k_j
print('ka,kw,kh',ka,kw,kh)

pww = 2.085  #  partial width of w boson 
gam_2 = (pww / mhch)**2 # (small gamma)_j 
#####################################
#   Function of GAW

# naive thought
def g(i,j): #G_ij function (ex. i for ka, j for kw)
    
    lam_ij = - 1.0 + 2 * i + 2 * j - (i - j)**2   # lamda_ij
    chunk_1 = PI / 2 + np.arctan((j * (1.0 - j + i) - lam_ij) \
                                 / ((1.0 - i) * np.sqrt(lam_ij)))
    chunk_2 = 1 / 3 * (1.0 - i) * (5 * (1.0 + i) - 4 * j - 2 * lam_ij / j)
    return 1 / 4 * (2 * (j - i - 1.0)* np.sqrt(lam_ij) * chunk_1 + \
                    (lam_ij - 2 * i) * np.log(i) + chunk_2) 


print('g()',g(ka,kw),g(kw,kh)) 
#####################################
# F_ij(x1,x2) Density function 
#def F_ij(i,j,x1,x2):
#    numer = (1.0 - x1 ) * (1.0 - x2 ) - i
#    denom = (1.0 - x1 - x2 - i + j)**2 + j * gam_2    
#    return numer / denom

#####################################
# parital width of H+>AW* and H+>hW*
#def lamaw(i,j):
#    return 9.0 * gf**2 * mw**4 * mhch * g(i,j) /(16.0 * PI**3)
#print(lamaw(kw,ka))
#def lamaw(x1,x2):
    
#    return 9.0 * gf**2 * mw**4 * mhch * 


def lamhw(i,j):
    return 9.0 * gf**2 * mw**4 * mhch * g(i,j) /(16.0 * PI**3)
print('||',lamhw(kw,kh))
######################################
###########################
#x2 = 2.0 * 80/mhch
#x2 = np.round(random.uniform(0.0, 1- ka), 100)
# For H+>AW approach
x2_min = 0.0
x2_max = 1 - ka
x1_min = 1 - x2 - ka
x1_max = 1 - ka / (1 - x2)
############################
# For H+>hW approach
#x1_min = 1 - x2 - kh
#x1_max = 1 - kh / (1 - x2)
#x2_min = 0.0
#x2_max = 1 - kh
#print(x1_min,x1_max,x2_min,x2_max)
#############################

        

#print( result.summary())
#print( 'result = %s    Q = %.2f' % (result, result.Q))

#########################
#Intergrand
def Faw(x):
        # For H+>AW approach       
#        x2_min = 0.0
#        x2_max = 1 - ka
        x2_min = 0.0
        x2_max = 1 - ka
        x2 = (x2_max-x2_min) * x[1]
        x1_min = 1 - x2 - ka
        x1_max = 1 - ka / (1 - x2)
        x1 = x1_min +(x1_max-x1_min) * x[0]
#        jac = (x1_max-x1_min)*(x2_max-x2_min)
        numer = (1.0 - x1 ) * (1.0 - x2 ) - ka
        denom = (1.0 - x1 - x2 - ka + kw)**2 + kw * gam_2  
        return  numer/denom #* jac
integ = vegas.Integrator([[x1_min, x1_max], [x2_min,x2_max]])
#FHpaw = integ( Faw,nitn=100, neval=5000) 
integral = 0.0
for x, wgt in integ.random():
    integral += wgt * Faw(x)

print(9.0 * gf**2 * mw**4 * mhch * integral / (16 * PI**3) )
#print(9.0 * gf**2 * mw**4 * mhch * FHpaw / (16 * PI**3) )
#print(9.0 * gf**2 * mw**4 * mhch * g(ka,kw) / (16 * PI**3) )
