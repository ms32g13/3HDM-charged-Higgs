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
from scipy import special as sp
#import vegas
### INVARIANT VALUE (PHYSICS VALUE OF SM PARAMATERS)
mt = 171.2    # mass of top quark
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
############
#QCD running coupling constant (alp) at energy scale MH. Relevant                                                                                                                                                                        
#for QCD corrections to Higgs decay widths.
############
coeffmc = 12.0/25.0
coeffmb = 12.0/23.0
alpmz = 0.1185 #alpha at z 
mu_w = 246.0 # scale at mu_w
mu_b = 5.0 # scale of mu_b 
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
    return c0_7sm(mu_w) + np.array(i)**2 * c0_7yy(mu_w) + np.array(j) * \
c0_7xy(mu_w)
def c0_8eff(mu_w,i,j):# i = Y^2 j = (XY^*) 
    return c0_8sm(mu_w) + np.array(i)**2 * c0_8yy(mu_w) + np.array(j) * \
c0_8xy(mu_w)
# LO Effective Wilson coefficient  # i = Y^2 j = (XY^*) 
def c0_eff(mu_w,i,j): # C0_2 effective = 1.0 C0_(1,3,4,5,6) = 0.0 
    c0_eff = []
    for i in np.arange(1.0,9.0,1.0):
        c0_eff.append(0.0)
    c0_eff[1] = 1.0 # c0_2eff
    c0_eff[6] = c0_7eff(mu_w,i,j) # c0_7eff
    c0_eff[7] = c0_8eff(mu_w,i,j) # c0_8eff
    return c0_eff
print(c0_eff(mu_w,20,1))


# NLO Effective Wilson coefficient  # i = Y^2 j = (XY^*) 
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
####
def Wi_sm(mu_w):#w7,8_sm values as array at NLO
    wi_sm = []
    for i in np.arange(1.0,9.0,1.0):
        wi_sm.append(0.0)
    wi_sm[6] = w7_sm(mu_w)
    wi_sm[7] = w8_sm(mu_w)
    return wi_sm
print(Wi_sm(mu_w))
def Mi_sm(mu_w):#m7,8_sm values as array at NLO
    mi_sm = []
    for i in np.arange(1.0,9.0,1.0):
        mi_sm.append(0.0)
    mi_sm[6] = m7_sm(mu_w)
    mi_sm[7] = m8_sm(mu_w)
    return mi_sm
def Ti_sm(mu_w):#t7,8_sm values as array at NLO
    ti_sm = []
    for i in np.arange(1.0,9.0,1.0):
        ti_sm.append(0.0)
    ti_sm[6] = t7_sm(mu_w)
    ti_sm[7] = t8_sm(mu_w)
    return ti_sm
def Wi_yy(mu_w):##w7,8_yy values as array at NLO
    wi_yy = []
    for i in np.arange(1.0,9.0,1.0):
        wi_yy.append(0.0)
    wi_yy[6] = w7_yy(mu_w)
    wi_yy[7] = w8_yy(mu_w)
    return wi_yy
def Mi_yy(mu_w):##m7,8_yy values as array at NLO
    mi_yy = []
    for i in np.arange(1.0,9.0,1.0):
        mi_yy.append(0.0)
    mi_yy[6] = m7_yy(mu_w)
    mi_yy[7] = m8_yy(mu_w)
    return mi_yy
def Ti_yy(mu_w):##t7,8_yy values as array at NLO
    ti_yy = []
    for i in np.arange(1.0,9.0,1.0):
        ti_yy.append(0.0)
    ti_yy[6] = t7_yy(mu_w)
    ti_yy[7] = t8_yy(mu_w)
    return ti_yy
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
def Wi_xy(mu_w):##w7,8_xy values as array at NLO
    wi_xy = []
    for i in np.arange(1.0,9.0,1.0):
        wi_xy.append(0.0)
    wi_xy[6] = w7_xy(mu_w)
    wi_xy[7] = w8_xy(mu_w)
    return wi_xy
def Mi_xy(mu_w):##m7,8_xy values as array at NLO
    mi_xy = []
    for i in np.arange(1.0,9.0,1.0):
        mi_xy.append(0.0)
    mi_xy[6] = m7_xy(mu_w)
    mi_xy[7] = m8_xy(mu_w)
    return mi_xy
def Ti_xy(mu_w):##t7,8_xy values as array at NLO
    ti_xy = []
    for i in np.arange(1.0,9.0,1.0):
        ti_xy.append(0.0)
    ti_xy[6] = t7_xy(mu_w)
    ti_xy[7] = t8_xy(mu_w)
    return ti_xy
#############################
#### Wilson coefficient at low scale(mu_b)