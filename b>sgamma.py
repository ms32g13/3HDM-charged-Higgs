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
import matplotlib.ticker as ticker
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
mu_w = 0.119 # alpha at w
##########################################################################
#4*pi*SQRT(2); factor appears in partial width of 2 fermions
fac = 4.0 * np.sqrt(2.0) * PI 
############################
#############################
#LO wilson coefficient at matching scale(mu_w)
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
# LO Effective Wilson coefficient 
c0_2eff  = 1.0 # C_2 effective = 1.0 C_(1,3,4,5,6) = 0.0
def c0_7eff(mu_w,i,j):
    return c0_7sm(mu_w) + np.array(i)**2 * c0_7yy(mu_w) + np.array(j) * \
c0_7xy(mu_w)

def c0_8eff(mu_w,i,j):
    return c0_8sm(mu_w) + np.array(i)**2 * c0_8yy(mu_w) + np.array(j) * \
c0_8xy(mu_w)
#############################
#NLO wilson coefficient at low scale(mu_b)