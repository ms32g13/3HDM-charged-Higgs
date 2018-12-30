#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:05:15 2018

@author: apple
"""
import math
import random
import numpy as np
import matplotlib.pyplot as plt
#import vegas
#import gvar as gv
import itertools
########################################################
### INVARIANT VALUE (PHYSICS VALUE OF SM PARAMATERS)
mt = 171.2    # mass of top quark
mz = 91.18 # 95.0 #98.14954576223639# 101.59443179342345  #91.18    # mass of z boson
nf = 5.0      # origin 5
PI = math.pi  # pi number 
mw = 80.385   # 80.385 # mass of w boson
mtau = 1.7771 # mass of tau neutrino
gf = 0.0000116639 # fermi constant
pwz = 2.4952  # partial width of z boson
pww = 2.085  #  partial width of w boson 
alpha_electroweak = 1.0 / 128.0 # alpha for sigma_0 
#mhch = 130.00  # mass of charged higgs
ma = 130.00 # mass of Pseudoscalar 
e_blist = np.arange(0.5,0.75,0.1) # b-tag efficiency 0.7
e_b = float(input('e_b value:'))
e_clist = np.arange(0.01,0.11,0.01) # c-tag efficiency 0.1
e_c = float(input('e_c value:'))
e_4jet = 0.841 # charged Higgs after 4-jet selection
e_mass = 0.718 / 0.841 #charged Higgs after mass-selection
e_antiww = 0.536 / 0.578 # charged Higgs after anti-ww background
e_l = 0.01 #light jet tag efficiency
e_isignal = 0.9  # invariant mass cut on signal
e_iback = 0.1 # invariant mass cut on background
epsilon = e_mass * e_4jet * e_isignal # total epsilon
############
# CKM elements
vcs = 0.97
vcb = 0.04
###########
massrange = np.arange(90.00,91.00,1.00)# charged Higgs ranged values
print('charH_mass:',massrange)
costhetaw = mw / mz                     # cos (weinberg angle : thetaw)
sinthetaw = math.sqrt(1 - (costhetaw) ** 2)      # sin(weinberg angle: thetaw)
#S = np.arange(180.00 , 219.00 , 10.00) ** 2  # Centre of Mass Energy ^ 2
########################################################
# DELPHI: (CENTRE OF MASS ENERGY= sqrt(S))^2
S_tuple = (189.0**2,192.0**2,196.0**2,200.0**2,202.0**2,205.0**2,206.3**2,206.6**2)
S = np.array(S_tuple)
#squared_s = np.arange(180.0,301.0,1.0)
#S = np.array(squared_s)**2
#print(squared_s)
#######################################################################
Lumin_h = (158.0,25.9,76.9,84.3,41.1,75.6,87.8,60.8) #DELPHI Hadronic Luminosicity (pb^(-1))
print(sum(Lumin_h))
Lumin_l = (153.8,24.5,72.4,81.8,39.4,69.1,79.8,50.0) #DELPHI Leptonic Luminosicity (pb^(-1))
print(sum(Lumin_l))
#######################################################################
def matrixmul(x,y):#matrix multiplication
    return np.array(list(itertools.chain(*np.outer(x,y))))
##############################################################
##################Cross-section for e+e- > w+w-
# Nuclear Physics B124 (1977) 511- 520 (Production of ww from e+e-)
# Cross-section of e+e- > w+w- 
def cross_section_eeww():
    x = 0.20
    E = np.sqrt(np.array(S) / 4.0)
    beta = np.sqrt(1.0 - mw**2 / E**2)
    y = np.array(S) / mw**2
    L = np.log((1.0 + beta) / (1.0 - beta))
    chunk_1 = PI * alpha_electroweak**2 * beta / (2.0 * x**2 * np.array(S))
    chunk_2 = (1.0 + 2.0/y + 2.0 /y**2) * L / beta
    chunk_3 = mz**2 * (1 - 2.0 * x) / (np.array(S) - mz**2)
    chunk_4 = 2.0 / y**2 * (1.0 + 2.0 * y) * L / beta - y/12.0 - 5/3.0 - 1.0/y
    chunk_5 = mz**4 * (8.0*x**2 - 4.0*x + 1.0) * beta**2 / (48.0 * (np.array(S) - mz**2)**2)
    chunk_6 = y**2 + 20.0 * y + 12.0
    return chunk_1 *(chunk_2 - 5.0/4.0 + chunk_3 * chunk_4 + chunk_5 * chunk_6)
#################################
#cross-section for e+e->w+w- , /(2.56819 * 10 **(-9)) : convert to pb^(-1)
#print('Single cross-section of ww:',cross_section_eeww() / (2.56819 * 10 **(-9)))
#print('Total cross-section of e+e->w+w-',sum(cross_section_eeww())) 
def eeww_event():# e+e->ww total bakcground events
    summ_ww =[]
    for k in np.arange(0,len(S)):
        # Single(cross-section * luminosicity)
        product_ww = cross_section_eeww()[k] / (2.56819 * 10 **(-9)) * Lumin_h[k]
        #Put all Single(cross-section * luminosicity) into list
        summ_ww.append(product_ww)
        print('|| e+e- > w w background:',sum(summ_ww),'events')
    return sum(summ_ww)# SIGMA(cross-section * luminosity)
#print('total eeww event',eeww_event())# total eeww background events
#print('||||After multiply 0.7^2:',eeww_event() * 0.7**2,'events') #5666.2838797003815
#################################################################
def wwcscs_background():# ww>cscs background
    chunk_1 = 11779 * 0.67**2 #  11779 total ww with hadronic-decay background
    chunk_2 = chunk_1 * 0.77 # after 4-jet preselection
    chunk_3 = chunk_2 * 2506.2 / 4076.9 # after anti-qqbar
    chunk_4 = chunk_3 / 4 # after ww>cscs
    return chunk_4 
print('|ww background:',wwcscs_background())
def backgroundtagging():# ww>cscs
    #wwcscs with e_c , or 3389.6
    result = 3389.6 * 1855.5 / 2506.2 / 4  * e_c**2
    # Z>b b-bar after b-tag .0.15 from z>bbbar oit of Z-decay ; e_iback = 0.1 invariant mass cut
    qqbar = 964.6 * 0.15 * e_b**2 * e_iback
    return result + qqbar
def backgroundnotagging():
    #wwcscs with e_c , or 3389.6
    result_1 = 3389.6 * 1855.5 / 2506.2 
    # Z>b b-bar after b-tag .0.15 from z>bbbar oit of Z-decay ; e_iback : 0.1 invariant mass cut
    qqbar_1 = 964.6 * e_iback
    return result_1 + qqbar_1
print( backgroundnotagging())
print('WW>cscs and qq_bar background:',backgroundtagging())
def eewwcscs_2bsignal():
    brcs = 0.8
    return wwcscs_background() * brcs**2 * 0.0
#print('||||| Singal of W_cscs after 2b-tagging:',eewwcscs_2bsignal(),'events')
def eewwcscs_1bsignal():
    brcs = 0.2
    return wwcscs_background() * brcs**2 * 0.0
#print('||||| Singal of W_cscs after 1b-tagging:',eewwcscs_1bsignal(),'events')
def eewwcscs_0bsignal():
    brcs = 0.8
    e_l = 0.01
    chunk2 = 2 * e_c * e_l * (1 - e_c) * (1 - e_l) + e_c**2 * (1 - e_l)**2 + \
    e_l**2 * (1 - e_c)**2
    return wwcscs_background() * brcs**2 * chunk2
#print('||||| Singal of W_cscs after 0b-tagging:',eewwcscs_0bsignal(),'events')
def eewwcstn_1bsignal():
    brtn = 0.2
    brcs = 0.8
    return wwcscs_background() * brcs * brtn * 0.0
#print('||||| Singal of W_cstn after 1b-tagging:',eewwcstn_1bsignal(),'events')
def eewwcstn_0bsignal():
    brtn = 0.3
    brcs = 0.6
    e_l = 0.01
    chunk1 = e_c * (1 - e_l) + e_l * (1 - e_c)
    return wwcscs_background() * brcs * brtn * chunk1
#print('||||| Singal of W_cstn after 0b-tagging:',eewwcstn_0bsignal(),'events')
#################################################################
#################################################################
#e+e- > H+H- Cross-section     
#A.G.Akeroyd Nuclear Physics B447(1995)3-17  EQUATIONS(15 - 17)
#F_HIGGS(s, mz,partialwidth_z,weinberg angle) function
def F_HIGGS(i):
    C_A = - 1.0 / (4 * sinthetaw * costhetaw)  # C_A function
    C_V = (1.0 - 4 * (sinthetaw ** 2)) / (4 * sinthetaw * costhetaw) #C_V function
    C1_V = (- 1.0 + 2 * (sinthetaw ** 2 )) / (2 * sinthetaw * costhetaw) #C'_V function
#    print('CACVC"V',C_A,C_V,C1_V)
    # denominator of F(s, mz,partialwidth_z,weinberg angle)
    denom_SZ = (i - mz ** 2) ** 2 + mz ** 2 * pwz ** 2
    return 1 - 2 * C_V * C1_V * (i * (i - mz ** 2)) / denom_SZ + \
           C1_V ** 2 * (C_V ** 2 + C_A ** 2) * (i ** 2)/ denom_SZ
########################################################
#Charged Higgs cross-section produced from e+e- >H+H- based on CoM Energy^2(S).
def cross_section_eeHH():
    ListeechaHH = [] # a list for same H+ mass, differ COM crosssection H+H-
    for j in massrange:#char_H mass
#    Summ_HH = 0.0 #Total(cross-section * Luminosity)
        for i in np.array(S):# COM Energy^2
#            print('i',i)
            BETA = math.sqrt(1.0 - 4 * (j)**2 / i)# betafunction with mhch,S
#        print('BETA:',j,BETA)
        # eechaHH , Production cross section of e+e- > H+H- 
            eechaHH = (1.0 / 4) * (4.0 * PI * alpha_electroweak **2 ) / (3.0 * i) * BETA**3 * F_HIGGS(i)
        # eechaHH / (2.56819 * 10 **(-9)) convert pb^(-1) to GeV
            sigma_eeHH = eechaHH / (2.56819 * 10 **(-9)) #crosssection of e+e->H+H-,convert to GeV
#            print('Production cross section of e+e- > H+H- :',j,i, sigma_eeHH,'GeV**2')
            ListeechaHH.append(sigma_eeHH) #put same H+ mass, differ COM crosssection H+H- together
    return ListeechaHH
########################################################
#calcualte e+e- > H+H- signal events from cross-section
def eeHH_event():
    Summ_HH = []
    eventHH = []
    for k in np.arange(0,len(np.array(S))*len(massrange)):
        #expand luminosity to fit length of chargedHiggs mass length
        Lumin_h_array = np.array(Lumin_h *len(massrange))
        #Single(cross-section * Luminosity)
        Product_HH = np.array(cross_section_eeHH())[k] * Lumin_h_array [k]  
        Summ_HH.append(Product_HH) # List_(cross-section * Luminosity)
    #split whole List_(cross-section * Luminosity) to different charged Higgs mass 
    #region in order to calculate signal events for single chargedHiggs mass
    eecharH_seperate = np.split(np.array(Summ_HH),len(massrange))
    #calculate single charged_higgs mass total signals from SIGMA(cross-section * luminosity)
#    print('SIGMA(cross-section * luminosity) events for each charged Higgs mass:')
    for i in np.arange(0,len(massrange)):
    #SIGMA(cross-section * luminosity) of each charged Higgs mass produced signals
        eecharH = sum(eecharH_seperate[i]) 
        eventHH.append(eecharH)
        print(i,'charH_mass:',massrange[i],'GeV',eecharH,'events')
    return np.array(eventHH)
print('SIGNAL after H+H- > cbcb',eeHH_event() * 0.8**2 ,'events')
print('signal/sqrt(background) FOR cbcb decay',eeHH_event() * 0.8**2 * 0.7**2 / \
      np.sqrt(backgroundtagging() * 0.1**2))
print('signal/sqrt(background) FOR no-tagging ',eeHH_event()* 1.0**2* 0.536 / \
      np.sqrt(1855.5))
# epsilon_signal = 0.536 * 71.8 / 57.8 = 0.6658269896193771
######################################################################
#with b signal of charged Higgs event:
def eeHHcbcb_2bsignal(x,y): #2real b 0 fake b
    chunk1 = eeHH_event() * x * y * \
    e_b**2 * (1 - e_c)**2# before compare with table 
    return chunk1 # * epsilon  # after 3 selections chosen
#print('||||| Singal of char_H_cbcb after 2b-tagging:',eeHHcbcb_2bsignal(),'events')
######################################################################
def eeHHcbcb_1bsignal(x,y):#1real b 1 fake b
    chunk1 = eeHH_event() * x * y * 4.0 * \
    e_b * e_c * (1 - e_b) * (1 - e_c)
      # before compare with table 
    return chunk1 #* epsilon # * e_antiww # after 3 selections chosen
#print('||||| Singal of char_H_cbcb after 1b-tagging:',eeHHcbcb_1bsignal())
######################################################################
def eeHHcbcb_0bsignal(x,y):#0real b 2 fake b
    chunk1 = eeHH_event() * x * y * e_c**2 * (1 - e_b)**2# before compare with table 
    return chunk1 # * epsilon # * e_antiww # after 3 selections chosen
#print('||||| Singal of char_H_cbcb after 0b-tagging:',eeHHcbcb_0bsignal())
######################################################################
def eeHHcbcs_2bsignal(x,y):#2real b 0 fake b
    return eeHH_event() * x * y * 0.0 #* epsilon #* e_antiww # after 3 selections chosen
#print('||||| Singal of char_H after 2b-tagging:',eeHHcbcs_2bsignal())
######################################################################
def eeHHcbcs_1bsignal(x,y):#1real b 1 fake b
    e_l = 0.01
    chunk2 = e_b * e_c * (1 - e_c) * (1 - e_l) + e_b * e_l * (1 - e_c)**2
    chunk1 = eeHH_event() * x * y * chunk2  # before compare with table 
    return chunk1 #* epsilon #* e_antiww # after 3 selections chosen
#print('||||| Singal of char_H_cbcs after 1b-tagging:',eeHHcbcs_1bsignal())
######################################################################
def eeHHcbcs_0bsignal(x,y):#0real b 2 fake b
    e_l = 0.01
    chunk2 = 2 * e_c * e_l * (1 - e_b) * (1 - e_c) + e_c**2 * (1 - e_b) * (1 - e_l)
    chunk1 = eeHH_event() * x * y * chunk2  # before compare with table 
    return chunk1# * epsilon #* e_antiww # after 3 selections chosen
#print('||||| Singal of char_H_cbcs after 0b-tagging:',eeHHcbcs_0bsignal())
    
######################################################################
def eeHHcscs_0bsignal(x,y): #0real b 2 fake b
    e_l = 0.01
    chunk2 = 2 * e_c * e_l * (1 - e_c) * (1 - e_l) + e_c**2 * (1 - e_l)**2 + \
    e_l**2 * (1 - e_c)**2
    chunk1 = eeHH_event() * x * y * chunk2  # before compare with table 
    return chunk1 #* epsilon #* e_antiww  after 3 selections chosen
#print('||||| Singal of char_H_cscs after 0b-tagging:',eeHHcscs_0bsignal())
######################################################################
def eeHHcbtn_1bsignal(x,y):
#    brcb = 0.8
#    brtn = 0.2
    chunk1 = eeHH_event() * x * y * \
    e_b * (1 - e_c)  # before compare with table 
    return chunk1 #* epsilon #* e_antiww  after 3 selections chosen
#print('||||| Singal of char_H_cbtn after 1b-tagging:',eeHHcbtn_1bsignal())
######################################################################
def eeHHcbtn_0bsignal(x,y):
#    brcb = 0.8
#    brtn = 0.2
    chunk1 = eeHH_event() * x * y * \
    e_b * (1 - e_c) # before compare with table 
    return chunk1 #* epsilon #* e_antiww  after 3 selections chosen
#print('||||| Singal of char_H_cbtn after 0b-tagging:',eeHHcbtn_0bsignal())
######################################################################
def eeHHcstn_1bsignal(x,y):
#    brcs = 0.8
#    brtn = 0.2
    chunk1 = eeHH_event() * x * y * 0.0  # before compare with table 
    return chunk1 #* epsilon# * e_antiww  after 3 selections chosen 
#print('||||| Singal of char_H_cstn after 1b-tagging:',eeHHcstn_1bsignal())
######################################################################
def eeHHcstn_0bsignal(x,y):
#    brcs = 0.3
#    brtn = 0.4
    e_l = 0.01 # e_s
    chunk2 = e_c * (1 - e_l) + e_l * (1 - e_c)
    chunk1 = eeHH_event() * x * y * chunk2  # before compare with table 
    return chunk1 #* epsilon# * e_antiww  after 3 selections chosen 
#print('||||| Singal of char_H_cstn after 0b-tagging:',eeHHcstn_0bsignal())
################################
###########################################################################
#KKKKK = PI * alpha_weak**2 / (np.array(S) * sinthetaw**4) *\
#       (2 * np.log(np.array(S) / mw**2) - 5 / 2 - 1 /(3 * costhetaw**2) + 5.0/(24.0 *costhetaw**4))
#print('KKKKK',KKKKK)
plt.figure()
plt.plot(np.sqrt(np.array(S)), cross_section_eeww()/ (2.56819 * 10 **(-9)))
plt.ylabel('Cross-section in pb^-1')
plt.xlabel('CoM energy Squared root of S')
plt.savefig('picture.png')
plt.show()
####################################################################################
################################################################################### 
def boson_decay_exist():  #check bosonic decay of charged Higgs exist or not and get lamaw value
    for i in massrange:                                                                                                                                                                                              
        if i > ma:
           gaw = 0.1
           lamaw=9.0 * gf**2 * mw**4 * i * gaw /(16.0 * PI**3)
        else:
           lamaw=0.00
    return lamaw
# BR : charged Higgs - A ,W
#braw=LamAW/(LamTN+LamCS+LamCB+LamAW) 
############################################################################
####################################################################################
#############################################################
def sig(x,y):# Signal/sqrt(Background)
    return matrixmul(x,1/np.sqrt(y))
print(sig(eeHH_event() * 0.8**2,backgroundtagging()))
print(matrixmul(e_blist**2,(1 - e_clist)**2),len(e_blist),len(e_clist))
print(sig(eeHH_event() * 0.4*0.6 * matrixmul(e_blist - e_blist**2,e_clist - e_clist**2),backgroundtagging()* e_clist**2),\
      len(sig(eeHH_event() * 0.4*0.6 * matrixmul(e_blist**2,(1 - e_clist)**2),backgroundtagging()* e_clist**2)))
plt.figure()
Contoursigma1 = plt.contour(e_clist,e_blist, \
           np.resize(sig(eeHH_event() * 0.8**2 * e_blist**2,backgroundtagging()* e_clist**2),\
            len(sig(eeHH_event() * 0.8**2 * e_blist**2,backgroundtagging() * e_clist**2))).reshape(len(e_blist),len(e_clist))\
                    )# ,levels = np.arange(1.0,6.0,1.0))
plt.title('Significance of H+H->cbcb / sqrt(background) with e_b and e_c')
plt.xlabel('e_c')# x-axis label
plt.ylabel('e_b')# y-axis label
plt.colorbar(Contoursigma1)

plt.figure()
Contoursigma2 = plt.contour(e_clist,e_blist, \
           np.resize(sig(eeHH_event() * 0.8**2 * matrixmul(e_clist**2 , (1 - e_blist)**2),backgroundtagging()* e_clist**2),\
            int(len(sig(eeHH_event() * 0.8**2  * matrixmul(e_clist**2 , (1 - e_blist)**2),backgroundtagging()* e_clist**2))/len(e_clist))).\
                     reshape(len(e_blist),len(e_clist)))# ,levels = np.arange(1.0,6.0,1.0))
plt.title('Significance of H+H->cbcb / sqrt(background) with e_b and e_c')
plt.xlabel('e_c')# x-axis label
plt.ylabel('e_b')# y-axis label
plt.colorbar(Contoursigma2)
#plt.figure()
#Contoursignal2 = plt.contour(e_c,e_b, \
#           np.resize(eeHHcbcb_1bsignal(),len(eeHHcbcb_1bsignal())).reshape(len(e_b),len(e_c)))
#plt.colorbar(Contoursignal2)