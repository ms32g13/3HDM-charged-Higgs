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
#######################
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
#######################
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
def Wi_sm(mu_w):#w7,8_sm function at NLO
    wi_sm = []
    for i in np.arange(1.0,9.0,1.0):
        wi_sm.append(0.0)
    wi_sm[6] = w7_sm(mu_w)
    wi_sm[7] = w8_sm(mu_w)
    return wi_sm
print(Wi_sm(mu_w))
# NLO Effective Wilson coefficient  # i = Y^2 j = (XY^*) 

#############################
#### Wilson coefficient at low scale(mu_b)