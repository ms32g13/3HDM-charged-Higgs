#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:05:15 2018

@author: apple
"""
import math
import random
from math import *
import numpy as np
import matplotlib.pyplot as plt
import vegas
import gvar as gv
########################################################
### INVARIANT VALUE (PHYSICS VALUE OF SM PARAMATERS)
mt = 171.2    # mass of top quark
mz = 91.18    # mass of z boson
nf = 5.0      # origin 5
PI = math.pi  # pi number 
mw = 80.33    # mass of w boson
mtau = 1.7771 # mass of tau neutrino
gf = 0.0000116639 # fermi constant
pwz = 2.4952  # partial width of z boson
pww = 2.085  #  partial width of w boson 
alpha_weak = 1.0 / 137 # alpha for sigma_0 
#mhch = 130.00  # mass of charged higgs
ma = 200.00 # mass of Pseudoscalar 
e_b = 0.55 # b-tag efficiency ratio
e_HHchar = 0.82# charged_higgs selection efficiency
############
# CKM elements
vcs = 0.97
vcb = 0.04
###########
massrange = np.arange(10.00,11.00,1.00)# charged Higgs ranged values
print('charH_mass:',massrange)
costhetaw = mw / mz                     # cos (weinberg angle : thetaw)
sinthetaw = math.sqrt(1 - (costhetaw) ** 2)      # sin(weinberg angle: thetaw)
#S = np.arange(180.00 , 219.00 , 10.00) ** 2  # Centre of Mass Energy ^ 2
########################################################
# DELPHI: (CENTRE OF MASS ENERGY= sqrt(S))^2
#S = np.array(180.0,192.0,196.0,200.0,202.0,205.0,206.3,206.6)**2
squared_s = np.arange(180.0,206.0,1.0)
S = (np.array(squared_s)**2)
print(squared_s)

#######################################################################
Lumin_h = (158.0,25.9,76.9,84.3,41.1,75.6,87.8,60.8) #DELPHI Hadronic Luminosicity (pb^(-1))
print(sum(Lumin_h))
Lumin_l = (153.8,24.5,72.4,81.8,39.4,69.1,79.8,50.0) #DELPHI Leptonic Luminosicity (pb^(-1))
########################################################
#e+e- > H+H- Cross-section     
#A.G.Akeroyd Nuclear Physics B447(1995)3-17  EQUATIONS(15 - 17)
C_A = - 1.0 / (4 * sinthetaw * costhetaw)  # C_A function
C_V = (1.0 - 4 * (sinthetaw ** 2)) / (4 * sinthetaw * costhetaw) #C_V function
C1_V = (- 1.0 + 2 * (sinthetaw ** 2 )) / (2 * sinthetaw * costhetaw) #C'_V function
print('CACVC"V',C_A,C_V,C1_V)
#F(s, mz,partialwidth_z,weinberg angle) function
def F_HIGGS(i):
    # denominator of F(s, mz,partialwidth_z,weinberg angle)
    denom_SZ = (i - mz ** 2) ** 2 + mz ** 2 * pwz ** 2
    return 1 - 2 * C_V * C1_V * (i * (i - mz ** 2)) / denom_SZ + \
           C1_V ** 2 * (C_V ** 2 + C_A ** 2) * (i ** 2)/ denom_SZ

#######################################################################
##################Cross-section for e+e- > w+w-
# Nuclear Physics B124 (1977) 511- 520 (Production of ww from e+e-)
# Cross-section of e+e- > w+w- 
def cross_section_eeww():
    E_squared = S /4 #E^2 
    print('E_squared',E_squared)
    beta = np.sqrt(1.0 - mw**2 / E_squared)#beta
    print('beta',beta)
    x = sinthetaw**2 # x = sin^2(weinberg angle)
#    x = 0.30
    Y =  np.array(S) / mw**2 #Y
    print('Y',Y)
#    Y = [ i / mw**2  for i in S ] #Y
    for i in np.arange(0,len(S)):
         L = np.log((1.0 + beta[i]) / (1.0 - beta[i])) #L
         
         L_over_beta = L / beta[i]
         #print('L',L,L_over_beta)
# np.asarray(): this can let list of numbers multiply float 
         chunk_1 =  PI * (alpha_weak)**2 * beta / (2.0 * x**2 * np.array(S))
#         print('chunk_1',chunk_1)
         chunk_2 = (1.0 + 2.0/Y + 2.0/(Y**2)) * L_over_beta
        # print('chunk_2',chunk_2)
         chunk_3 = (mz**2 * (1.0 - 2.0 * x)/(np.array(S) - mz**2)) * (2.0/(Y**2) * \
    (1.0 + 2 * Y) * L_over_beta - Y/12 - 5.0/3.0 - 1.0/Y)
        # print('chunk_3',chunk_3)
         chunk_4 = (mz**4 * (8.0 * x**2 - 4.0 * x + 1) * beta**2 / \
    (48.0 * (np.array(S) - mz**2 )**2)) * (Y**2 + 20.0 * Y + 12.0)
         #print('chunk_4',chunk_4)
    return chunk_1 * ( chunk_2 - 5.0/4.0 + chunk_3 + chunk_4)           
#################################
#cross-section for e+e->w+w- 
print('Single cross-section of ww:',cross_section_eeww() / ((2.56819 * 10 **(-9))))
print('----------------')
print('Total cross-section of e+e->w+w-',sum(cross_section_eeww())) 
############################################################
##################    Partial width  and Branching ratio for charged Higgs 
##################    based on different Models
#Higgs-fermion-fermion couplings in MHDM:
#couptn = Z**2
#coupc = Y**2
#coups = X**2
#coupb = X**2
#KKKKK = PI * alpha_weak**2  * np.sqrt(1.0 - mw**2  / np.array(S)) / (np.array(S) * sinthetaw**4) *\
#       (2 * np.log(np.array(S) / mw**2) - 5 / 2 - 1 /(3 * costhetaw**2) + 5/(24 *costhetaw**4))
#print('KKKKK',KKKKK)
plt.figure()
cross_CMS = plt.plot(np.sqrt(S), cross_section_eeww() / ((2.56819 * 10 **(-9))))
plt.show()
# Partial Width charged H - tau, neutrino(tau-type)
def lamtn (j,Z): # j for mhch
        couptn = Z**2
        return gf * j * (mtau**2) * couptn / fac
# Charged H - charm, strange quark
def lamcs (j,X,Y):   # j for mhch
        coupc = Y**2
        coups = X**2 
        return (3 * gf * j * (vcs**2)) * \
     (mc**2 * coupc + ms**2 * coups) *(1.0 + 17.0 * alpmz/(3.0 * PI)) / fac

# Charged H - charm, bottom quark
def lamcb (j,X,Y):
        coupc = Y**2
        coupb = X**2
        return (3 * gf * j * (vcb**2)) * \
     (mc**2 * coupc + mb**2 * coupb) * (1.0 + 17.0 * alpmz/(3.0 * PI)) / fac
# BR : charged Higgs - A ,W
#braw=LamAW/(LamTN+LamCS+LamCB+LamAW) 
#############################################################
# BR : charged Higgs - tau , neutrino(tau)
def brtn (j,X,Y,Z):
        return lamtn(j,Z)/(lamtn(j,Z)+lamcs(j,X,Y)+lamcb(j,X,Y)+lamaw)
# BR : charged Higgs - charm , strange quark
def brcs (j,X,Y,Z):
        return lamcs(j,X,Y)/(lamtn(j,Z)+lamcs(j,X,Y)+lamcb(j,X,Y)+lamaw)
# BR : charged Higgs - charm, bottom quark
def brcb (j,X,Y,Z):
        return lamcb(j,X,Y)/(lamtn(j,Z)+lamcs(j,X,Y)+lamcb(j,X,Y)+lamaw)
########################################################
#Charged Higgs cross-section produced from e+e- >H+H- based on CoM Energy^2(S).
for j in massrange:#char_H mass
    ListeechaHH = [] # a list for same H+ mass, differ COM crosssection H+H-
    print('_______________',j)
#    Summ_HH = 0.0 #Total(cross-section * Luminosity)
    for i in S:# COM Energy^2
#        print('i',i)
        BETA = math.sqrt(1.0 - 4 * (j)**2 / i)# betafunction with mhch,S
#        print('BETA:',j,BETA)
        # eechaHH , Production cross section of e+e- > H+H- 
        eechaHH = (1.0 / 4) * (4.0 * PI * alpha_weak **2 ) / (3.0 * i) * BETA**3 * F_HIGGS(i)
        # eechaHH / (2.56819 * 10 **(-9)) convert pb^(-1) to GeV
        sigma_eeHH = eechaHH / (2.56819 * 10 **(-9)) #crosssection of e+e->H+H-,convert to GeV
#        print('Production cross section of e+e- > H+H- :',j,i, sigma_eeHH,'GeV**2')
        ListeechaHH.append(eechaHH/ (2.56819 * 10 **(-9))) #put same H+ mass, differ COM crosssection H+H- together
#    print(len(ListeechaHH),ListeechaHH)#  
    Summ_HH = []#Total(cross-section * Luminosity)
    print('ListeechaHH',ListeechaHH)
#    for k in np.arange(0,len(S)):
#            Product_HH = ListeechaHH[k] * Lumin_h[k]  #Single(cross-section * Luminosity)
#            print('||||','m_charH',k,j,'S',S[k],'Luminosity',Lumin_h[k],Product_HH)
#            Summ_HH.append(Product_HH) # Sum_(cross-section * Luminosity)
#    print('Total(cross-section * Luminosity_h)_charH:', sum(Summ_HH),'events')
###########################################################################
    
###########################################################################
############
#QCD running coupling constant (alp) at energy scale MH. Relevant                                                                                                                                                                        
#for QCD corrections to Higgs decay widths.
############
    coeffmc = 12.0/25.0
    coeffmb = 12.0/23.0
    alpmz = 0.1185 #coupling of z 
    alpmhch = alpmz /(1.0 + (alpmz/(12.0 * PI)) * (33.0 - 2 * nf) * \
              math.log((j**2)/(mz**2))) 

    alpmb = alpmz /(1.0 + (alpmz/(12.0 * PI))*(33.0 - 2 * nf) * \
           math.log(4.89**2/(mz**2)))
#        print('alpmb :',alpmb, ' alpmhch:',alpmhch)

# alpha_s at sale 2 GeV.Can read this off from Rev.Part.Properties
    alp2gev = 0.117 /(1.0 + (0.117 /(12.0 * PI)) * (33.0 - 2 * nf) *\
              math.log(2.0**2 / mz**2))
# alp2gev = 0.310
# for running of ms at 2 GeV to ms at 4.88 GeV (mb pole mass)                                                                                                                                                                                                                          
    fc2gev = 1.0 + 1.101 * alp2gev/PI + 1.39 * (alp2gev/PI)**2 + 1.09 * \
            (alp2gev/PI)**3
    fcmb = 1.0 + 1.101 * alpmb/PI + 1.39 * (alpmb/PI)**2 + 1.09 * \
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
#        print(mc , mb, ms)
##########################################################################
#4*pi*SQRT(2); factor appears in partial width of 2 fermions
    fac = 4.0 * math.sqrt(2.0) * PI 

###################################################################################                                                                                                                                                                                                  
    if j > ma:
       lamaw=9.0 * gf**2 * mw**4 * j * gaw /(16.0 * PI**3)
    else:
       lamaw=0.00


####################################################################################
#        eventHH = eechaHH * 500 * 10**9 / 2.6 
#                 print()#BETA**3,F_HIGGS(i))
#        print('centre of Mass',i**(1 / 2),'chargehiggs mass',j,\
#              'Production crosssection ',eechaHH,\
#              'Event rate for cbcb',eventHH )
####################################################################################
#        ### LOOP CALCULATION and  X, Y, Z dependence
#    x = np.arange(0.0,40.2,0.2)
#    y = np.arange(0.0,0.62,0.02)
#    z = np.arange(0.0,5.02,0.02)
#    for Y in y:
#       for X in x:
#for Z in np.arange(0.0,5.1,0.5):

    
#            Z = 1.0
#          for Z in z:          
#                 print('|',brcb(j,X,Y,0.1))
                 # expected number event for cbcb from H+H-
#            print('brcb',j,brcb(j,Y,X,Z))     
#############################################################
