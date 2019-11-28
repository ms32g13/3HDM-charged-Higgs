#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:41:45 2019

@author: muyuansong0419
"""
import math
###Invariant parameters 
### INVARIANT VALUE (PHYSICS VALUE OF SM PARAMATERS)
mt = 175 # mass of top quark  174.6
mc = 1.39  # mass of charm quark 1.29 / 1.279/1.275 # 1.3846 2hdm paper
mb = 4.78  # mass of bottom quark  #4.78 2hdm paper
md = 0.0048 # mass of down quark
mu = 0.0023 # mass of up quark
nf = 5.0      # num of active flavours
nc = 3        # repesentation of SU(3)
PI = math.pi  # pi number 
mw = 80.33 # mass of w boson #80.33 2hdm paper
mz = 91.18 # 95.0 #98.14954576223639 #91.18    # mass of z boson
pwz = 2.4952  # partial width of z boson
pww = 2.085  #  partial width of w boson 
mtau = 1.7771 # mass of tau lepton
gf = 1.1663787e-05#0.000011663787 # fermi constant
mhch = 100.00  # mass of charged higgs
mh = 125.0 #  mass of higgs
ma = 300.00 # mass of Pseudoscalar 
# CKM elements
s12_ckm = 0.2229  
s23_ckm = 0.0412 
s13_ckm = 0.0036 
c12_ckm = math.sqrt(1.0 - s12_ckm**2) 
c23_ckm = math.sqrt(1.0 - s23_ckm**2) 
c13_ckm = math.sqrt(1.0 - s13_ckm**2) 
Vud = c12_ckm * c13_ckm 	     
Vus = s12_ckm * c13_ckm
vub = s13_ckm
vud = -s12_ckm * c23_ckm - c12_ckm * s23_ckm * s13_ckm
vcs = c12_ckm * c23_ckm - s12_ckm * s23_ckm * s13_ckm
vcb = s23_ckm * c13_ckm
vtd = s12_ckm * s23_ckm - c12_ckm * c23_ckm * s13_ckm
vts = -c12_ckm * s23_ckm - s12_ckm * c23_ckm * s13_ckm
vtb = c23_ckm * c13_ckm
############
#QCD running coupling constant (alp) at energy scale MH. Relevant                                                                                                                                                                        
#for QCD corrections to Higgs decay widths.
############
######################################################
coeffmc = 12.0/25.0
coeffmb = 12.0/23.0
alpmz = 0.119 #coupling of z #0.119 2hdm paper
alpha_electroweak = 1.0 / 130.3 # alpha for EM coupling
# CKM matrix elements observables from CKMfitter
#https://indico.cern.ch/event/684284/contributions/2952455/attachments/1719296/2774804/Vale_Silva_3.pdf
rho_bar_ckm = 0.1577
eta_bar_ckm = 0.3493
rho_ckm = 0.135
eta_ckm = 0.349
lanmda_ckm = 0.224747
A_ckm = 0.8403
lambda_t = vtb * vts # lambda_t = V_tb * (V_ts).conjugate  CKM products
squared_vtsvtbovervcb = 0.95 # 0.95 in Borzumati's paper / 0.9626 Andrew,Stefano
#################################################################
#Measured Semi- leptonic branching ratio B_SL for B>Xs + gamma result
B_SL = 0.1049 # Phys. Rev. Lett. 76, 1570 – Published 4 March 1996 =  0.1049