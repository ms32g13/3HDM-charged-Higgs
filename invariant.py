#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:41:45 2019

@author: muyuansong0419
"""
import math
###Invariant parameters 
### INVARIANT VALUE (PHYSICS VALUE OF SM PARAMATERS)
mt = 171.2    # mass of top quark
mc = 1.29     # mass of charm quark
mb = 4.18     # mass of bottom quark
nf = 5.0      # origin 5
PI = math.pi  # pi number 
mw = 80.385   # mass of w boson 80.33
mz = 91.18 # 95.0 #98.14954576223639# 101.59443179342345  #91.18    # mass of z boson
pwz = 2.4952  # partial width of z boson
pww = 2.085  #  partial width of w boson 
mtau = 1.7771 # mass of tau neutrino
gf = 0.0000116639 # fermi constant
#mhch = 130.00  # mass of charged higgs
mh = 125.0 #  mass of higgs
ma = 130.00 # mass of Pseudoscalar 
# CKM elements
vcs = 0.97
vcb = 0.04
vts = 0.0404
vtb = 0.999146
############
#QCD running coupling constant (alp) at energy scale MH. Relevant                                                                                                                                                                        
#for QCD corrections to Higgs decay widths.
############
######################################################
coeffmc = 12.0/25.0
coeffmb = 12.0/23.0
alpmz = 0.1185 #coupling of z 
alpha_electroweak = 1.0 / 128.0 # alpha for EM coupling
lambda_t = vtb * vts # lambda_t = V_tb * (V_ts).conjugate  CKM products
