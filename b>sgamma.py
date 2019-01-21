#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 14:58:53 2019
@author: muyuansong0419
"""

import math
import random
import numpy as np
import matplotlib.pyplot as plt
import gammaeffective as gaeff
from scipy import special as sp
#import vegas
### INVARIANT VALUE (PHYSICS VALUE OF SM PARAMATERS)
mt = 171.2    # mass of top quark
mb = 4.89     # mass of bottom quark
mc = 1.64     # mass of charm quark
mz = 91.18    # mass of z boson
nf = 5.0      # origin 5
PI = math.pi  # pi number 
mw = 80.385    # mass of w boson 80.33
mtau = 1.7771 # mass of tau neutrino
mh = 125 #  mass of higgs
gf = 0.0000116639 # fermi constant
#mhch = 130.00  # mass of charged higgs
ma = 200.00 # mass of Pseudoscalar 
# CKM elements
vcs = 0.97
vcb = 0.04
x = mt**2 / mw**2 # mt^2 /mw^2
y = mt**2 / mh**2 # mt^2 / mh^2
a_i = np.array([14 /23, 16 /23, 6 /23, - 12/23,0.4086,-0.4230,-0.8994,0.1456])#{a_i}
h_i = np.array([626126/272277 , - 56281/51730, - 3/7, - 1/14, - 0.6494, - 0.0380,\
                - 0.0186, - 0.0057])#{h_i}
a2_i = np.array([14/23, 0.4086, - 0.4230, - 0.8994, 0.1456])#{a'_i}
h2_i = np.array([313063/363036, - 0.9135, 0.0873, - 0.0571, 0.0209])#{h'_i}
e_i = np.array([4661194/816831, - 8516/2217, 0.0,0.0,- 1.9043, - 0.1008,\
                0.1216,0.0183])#{e_i}
f_i = np.array([- 17.3023, 8.5027, 4.5508, 0.7519, 2.0040, 0.7476, - 0.5385, 0.0914])#{f_i}
k_i = np.array([9.9372, - 7.4878, 1.2688, - 0.2925, -2.2923, - 0.1461, 0.1239, 0.0812])#{k_i}
l_i = np.array([0.5784, - 0.3921, -0.1429, 0.0476, - 0.1275, 0.0317, 0.0078, - 0.0031])#{l_i}

print(gaeff.gamma0eff())#gamma_0_eff_ji matrix values
print(gaeff.gamma1eff())#gamma_1_eff_ji matrix values
############
#QCD running coupling constant (alp) at energy scale MH. Relevant                                                                                                                                                                        
#for QCD corrections to Higgs decay widths.
############
coeffmc = 12.0/25.0
coeffmb = 12.0/23.0
alpmz = 0.1185 #alpha at z 
mu_w = 0.119 # scale at mu_w
mu_b = 5.0 # scale of mu_b 
def run_quark_bar(q):#running quark mass at scale mu_w under minimal subtract scheme
    c1 = np.log(q**2 / mw**2)
    c2 = mu_w / PI
    return q * (1 + c2 * c1 - 4 / 3 * c2 )
print(run_quark_bar(mb))

##########################################################################
#4*pi*SQRT(2); factor appears in partial width of 2 fermions
fac = 4.0 * np.sqrt(2.0) * PI 
############################
#############################
### Wilson coefficient at matching scale(mu_w)
####################### LO
def c0_7sm(mu_w):
    print('x',x)
    chunk = - 8 * x**3 + 3 * x**2 + 12 * x - 7
    chunk1 = (18 * x**2 - 12 * x) * np.log(x)
    return x / 24 * ( (chunk + chunk1) / (x - 1)**4 )
def c0_8sm(mu_w):
    chunk = - x**3 + 6 * x**2 - 3 * x - 2 - 6 * x * np.log(x)
    return x / 8 * ( chunk / (x - 1)**4 )
def c0_7yy(mu_w):
    chunk = - 8 * y**3 + 3 * y**2 + 12 * y - 7
    chunk1 = (18 * y**2 - 12 * y) * np.log(y)
    return y / 24 * ( (chunk + chunk1) / (y - 1)**4 )
def c0_8yy(mu_w):
    chunk = - y**3 + 6 * y**2 - 3 * y - 2 - 6 * y * np.log(y)
    return y / 8 * ( chunk / (y - 1)**4 )
def c0_7xy(mu_w):
    chunk = - 5 * y**2 + 8 * y  - 3
    chunk1 = (6 * y - 4) * np.log(y)
    return y / 12 * ( (chunk + chunk1) / (y - 1)**3)
def c0_8xy(mu_w):
    chunk = - y**2 + 4 * y - 3 - 2 * np.log(y)
    return y / 4 * ( chunk / (y - 1)**3 )
###
def c0_7eff(mu_w,i,j):# i = Y^2 j = (XY^*) 
    return c0_7sm(mu_w) + np.array(abs(i)**2)**2 * c0_7yy(mu_w) + np.array(j) * \
c0_7xy(mu_w)
def c0_8eff(mu_w,i,j):# i = Y^2 j = (XY^*) 
    return c0_8sm(mu_w) + np.array(abs(i)**2)**2 * c0_8yy(mu_w) + np.array(j) * \
c0_8xy(mu_w)
# LO Effective Wilson coefficient  # i = Y^2 j = (XY^*) 
def c0_eff(mu_w,i,j): # C0_2 effective = 1.0 C0_(1,3,4,5,6) = 0.0 
    c0_eff = []
    for n in np.arange(1.0,9.0,1.0):
        c0_eff.append(0.0)
    c0_eff[1] = 1.0 # c0_2eff
    c0_eff[6] = c0_7eff(mu_w,i,j) # c0_7eff
    c0_eff[7] = c0_8eff(mu_w,i,j) # c0_8eff
    return np.array(c0_eff)
print('c0_eff(mu_w,i,j)',c0_eff(mu_w,1.0,20))
####################### NLO
def E_0(mu_w):#NLO
    chunk_1 = x * (x**2 + 11 * x - 18) / (12 * (x - 1)**3)
    chunk_2 = x**2 * (4 * x**2 - 16 * x + 15) * np.log(x) / (6 * (x - 1)**4)
    chunk_3 = - 2 / 3 * np.log(x) - 2 / 3
    return chunk_1 + chunk_2 + chunk_3
def w7_sm(mu_w):#NLO
    chunk_1 = (- 16 * x**4 - 122 * x**3 + 80 * x**2 - 8 * x) / (9 * (x - 1)**4)
    chunk_2 = sp.spence(1.0 - 1.0/x)
    chunk_3 = (6 * x**4 + 46 * x**3 - 28 * x**2) * (np.log(x))**2 / (3 * (x - 1)**5) 
    chunk_4 = (- 102 * x**5 - 588 *x**4 - 2262 * x**3 + 3244 * x**2 - 1364 * x + 208 ) * \
    np.log(x) / (81 * (x - 1)**5)
    chunk_5 = (1646 * x**4 + 12205 * x**3 - 10740 * x**2 + 2509 * x - 436) / (486 * (x - 1)**4)
    return chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5
def w8_sm(mu_w):#NLO
    chunk_1 = (- 4 * x**4 + 40 * x**3 + 41 * x**2 + x)/(6 * (x - 1)**4)
    chunk_2 = sp.spence(1.0 - 1.0/x)
    chunk_3 = (- 17 * x**3 - 31 * x**2) * (np.log(x))**2 / (2 * (x - 1)**5)
    chunk_4 = (- 210 * x**5 + 1086 * x**4 + 4893 * x**3 + 2857 * x**2 - 1994 * x + 280) * \
    np.log(x) / (216 * (x - 1)**5)
    chunk_5 = (737 * x**4 - 14102 * x**3 - 28209 * x**2 + 610 * x - 508) \
    / (1296 * (x - 1)**4)
    return chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5
def m7_sm(mu_w):#NLO
    chunk_1 = 82 * x**5 + 301 * x**4 + 703 * x**3 - 2197 * x**2 + 1319 * x - 208
    chunk_2 = - (162 * x**4 + 1242 * x**3 + 756 * x**2) * np.log(x)
    return (chunk_1 + chunk_2) / (81 * (x - 1)**5)
def m8_sm(mu_w):#NLO
    chunk_1 = 77 * x**5 - 475 * x**4 - 1111 * x**3 + 607 * x**2 + 1042 * x - 140
    chunk_2 = ( 918 * x**3 + 1674 * x**2) * np.log(x)
    return (chunk_1 + chunk_2) / (108 * (x - 1)**5)
def t7_sm(mu_w):#NLO
    chunk_1 = 47 * x**3 - 63 * x**2 + 9 * x + 7 
    chunk_2 = - (18 * x**3 + 30 * x**2 - 24 * x) * np.log(x)
    return x / 3 * ( (chunk_1 + chunk_2) / (x - 1)**5 )
def t8_sm(mu_w):#NLO
    chunk_1 = - x**3 - 9 * x**2 + 9 * x + 1 + (6 * x**2 + 6 * x) * np.log(x)
    chunk_2 = (x - 1)**5
    return 2 * x * (chunk_1 / chunk_2)
#coupling Y^2
def E_H(mu_w):#NLO
    chunk_1 = 7 * y**3 - 36 * y**2 + 45 * y - 16
    chunk_2 = (18 * y - 12) * np.log(y)
    return 1 * y / 36 * ( (chunk_1 + chunk_2) / ((y - 1)**4) )
def w7_yy(mu_w):#NLO
    chunk_1 = (8 * y**3 - 37 * y**2 + 18 * y) / ((y - 1)**4)
    chunk_2 = sp.spence(1.0 - 1.0 / y)
    chunk_3 = (3 * y**3 + 23 * y**2 - 14 * y) * (np.log(y))**2 / ((y - 1)**5)
    chunk_4 = (21 * y**4 - 192 * y**3 - 174 * y**2 + 251 * y - 50) \
    * np.log(y) / (9 * (y - 1)**5)
    chunk_5 = (- 1202 * y**3 + 7569 * y**2 - 5436 * y + 797) / (108 * (y - 1)**4)
    return 2 * y / 9 * (chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5) - \
4 / 9 * E_H(mu_w)
def w8_yy(mu_w):#NLO
    chunk_1 = (13 * y**3 - 17 * y**2 + 30 * y) / ((y - 1)**4)
    chunk_2 = sp.spence(1.0 - 1.0 / y)
    chunk_3 = - (17 * y**2 + 31 * y) / ((y - 1)**5) * (np.log(y))**2
    chunk_4 = (42 * y**4 + 318 * y**3 + 1353 * y**2 + 817 * y - 226) * \
np.log(y) / (36 * (y - 1)**5)
    chunk_5 = (- 4451 * y**3 + 7650 * y**2 - 18153 * y + 1130) / (216 * (y - 1)**4)
    return 1 / 6 * y * (chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5) - \
1 / 6 * E_H(mu_w)
def m7_yy(mu_w):#NLO
    chunk_1 = - 14 * y**4 + 149 * y**3 - 153 * y**2 - 13 * y + 31
    chunk_2 = - (18 * y**3 + 138 * y**2 - 84 * y) * np.log(y)
    chunk_3 = (y - 1)**5
    return 1 / 27 * y * ( (chunk_1 + chunk_2) / chunk_3 )
def m8_yy(mu_w):#NLO
    chunk_1 = - 7 * y**4 + 25 * y**3 - 279 * y**2 + 223 * y + 38 
    chunk_2 = (102 * y**2 + 186 * y) * np.log(y)
    chunk_3 = (y - 1)**5
    return 1 / 36 * y * ( (chunk_1 + chunk_2) / chunk_3)
def t7_yy(mu_w):#NLO
    chunk_1 = 47 * y**3 - 63 * y**2 + 9 * y + 7 
    chunk_2 = - (18 * y**3 + 30 * y**2 - 24 * y) * np.log(y)
    return 1 / 3 * (y / 3 * ( (chunk_1 + chunk_2) / (y - 1)**5 ) )
def t8_yy(mu_w):#NLO
    chunk_1 = - y**3 - 9 * y**2 + 9 * y + 1 + (6 * y**2 + 6 * y) * np.log(y)
    chunk_2 = (y - 1)**5
    return 1 / 3 * ( 2 * y * (chunk_1 / chunk_2) )

#Coupling (XY*)
def w7_xy(mu_w):#NLO
    chunk_1 = (8 * y**2 - 28 * y + 12) / (3 * (y - 1)**3)
    chunk_2 = sp.spence(1.0 - 1.0 / y)
    chunk_3 = (3 * y**2 + 14 * y - 8) * (np.log(y))**2 / (3 * (y - 1)**4)
    chunk_4 = (4 * y**3 - 24 * y**2 + 2 * y + 6 ) * np.log(y) / (3 * (y - 1)**4)
    chunk_5 = (- 2 * y**2 + 13 * y  - 7) / ((y - 1)**3)
    return 4 / 3 * y * (chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5)
def w8_xy (mu_w):#NLO
    c1 = (17 * y**2 - 25 * y + 36) / (2 * (y - 1)**3)
    c2 = sp.spence(1.0 - 1.0 / y)
    c3 = - (17 * y + 19) * (np.log(y))**2 / ((y - 1)**4)
    c4 = (14 * y**3 - 12 * y**2 + 187 * y + 3) * np.log(y) / (4 * (y - 1)**4)
    c5 = - (3 * (29 * y**2 - 44 * y + 143)) / (8 * (y - 1)**3)
    return 1 / 3 * y * (c1 * c2 + c3 + c4 + c5)
def m7_xy(mu_w):#NLO
    c1 = - 8 * y**3 + 55 * y**2 - 68 * y + 21 
    c2 = - (6 * y**2 + 28 * y - 16) * np.log(y)
    c3 = (y - 1)**4
    return  2 / 9 * y * ( (c1 + c2) / c3 )
def m8_xy(mu_w):#NLO
    c1 = - 7 * y**3 + 23 * y**2 - 97 * y  + 81
    c2 = (34 * y + 38) * np.log(y)
    c3 = (y - 1)**4
    return 1 / 6 * y * ( (c1 + c2) / c3 )
def t7_xy(mu_w):#NLO
    c1 = 13 * y**2 - 20 * y + 7 - (6 * y**2 + 4 * y - 4) * np.log(y)
    c2 = (y - 1)**4
    return 2 / 3 * y * (c1 / c2)
def t8_xy(mu_w):#NLO
    c1 = - y**2 - 4 * y + 5 + (4 * y + 2) * np.log(y)
    c2 = (y - 1)**4
    return 2 * y * (c1 / c2)

####
def Wi_sm(mu_w):#w7,8_sm values as array at NLO
    return np.array([w7_sm(mu_w),w8_sm(mu_w)])
def Mi_sm(mu_w):#m7,8_sm values as array at NLO
    return np.array([m7_sm(mu_w),m8_sm(mu_w)])
def Ti_sm(mu_w):#t7,8_sm values as array at NLO
    return np.array([t7_sm(mu_w),t8_sm(mu_w)])
def Wi_yy(mu_w):##w7,8_yy values as array at NLO
    return np.array([w7_yy(mu_w),w8_yy(mu_w)])
def Mi_yy(mu_w):##m7,8_yy values as array at NLO
    return np.array([m7_yy(mu_w),m8_yy(mu_w)])
def Ti_yy(mu_w):##t7,8_yy values as array at NLO
    return np.array([t7_yy(mu_w),t8_yy(mu_w)])
def Wi_xy(mu_w):##w7,8_xy values as array at NLO
    return np.array([w7_xy(mu_w),w8_xy(mu_w)])
def Mi_xy(mu_w):##m7,8_xy values as array at NLO
    return np.array([m7_xy(mu_w),m8_xy(mu_w)])
def Ti_xy(mu_w):##t7,8_xy values as array at NLO
    return np.array([t7_xy(mu_w),t8_xy(mu_w)])
############
# NLO Effective Wilson coefficient  # i = Y^2 j = (XY^*) 
def c1_eff(mu_w,i,j):#c1,eff,i,sm with Y^2 and XY*
    ratio1 = np.log(mt**2 / mu_w**2)
    ratio = np.log(mu_w**2 / mw**2)
    list1 = []
    for n in np.arange(0.0,6.0,1.0):
        list1.append(0.0)
    list1[0] = Wi_sm(mu_w)[0] + Mi_sm(mu_w)[0] * ratio + \
        Ti_sm(mu_w)[0] * (ratio1  - 4 / 3)#7
    list1[1] = Wi_sm(mu_w)[1] + Mi_sm(mu_w)[1] * ratio + \
        Ti_sm(mu_w)[1] * (ratio1  - 4 / 3)#8
    list1[2] = Wi_yy(mu_w)[0] + Mi_yy(mu_w)[0] * ratio + \
        Ti_yy(mu_w)[0] * (ratio1  - 4 / 3)#7
    list1[3] = Wi_yy(mu_w)[1] + Mi_yy(mu_w)[1] * ratio + \
        Ti_yy(mu_w)[1] * (ratio1  - 4 / 3)#8
    list1[4] = Wi_xy(mu_w)[0] + Mi_xy(mu_w)[0] * ratio + \
        Ti_xy(mu_w)[0] * (ratio1  - 4 / 3)#7
    list1[5] = Wi_xy(mu_w)[1] + Mi_xy(mu_w)[1] * ratio + \
        Ti_xy(mu_w)[1] * (ratio1  - 4 / 3)#8
    print('list1', list1)
    c1_eff = []
    for n in np.arange(1.0,9.0,1.0):
        c1_eff.append(0.0)
    c1_eff[0] = 15 + 6 * ratio #c1,eff,1,sm
    c1_eff[3] = E_0(mu_w) + 2 / 3 * ratio + np.array(abs(i)**2) * E_H(mu_w)
    c1_eff[6] = list1[0] + np.array(abs(i)**2) * list1[2] + np.array(j) * list1[4]
    c1_eff[7] = list1[1] + np.array(abs(i)**2) * list1[3] + np.array(j) * list1[5]
    return np.array(c1_eff)
print('c1_eff(mu_w,i,j)',c1_eff(mu_w,1.0,20))
######Total Ceffective_i (mu_w) at matching scale 
def c_i_eff_muw(mu_w,i,j):
    return c0_eff(mu_w,i,j) + alpmz / (4.0 * PI) * c1_eff(mu_w,i,j)
print('c_i_eff_muw(mu_w,i,j)',c_i_eff_muw(mu_w,1.0,20))
#################################################################
##################################################################
#### Wilson coefficient at low scale(mu_b)

