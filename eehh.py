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
import vegas
import gvar as gv
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
#charH_mass = float(input('chaH mass value:'))
ma = 130.00 # mass of Pseudoscalar 
e_blist = np.arange(0.5,0.75,0.023) # b-tag efficiency 0.7
#e_b = float(input('e_b value:')) #
e_b = 0.7 # set b_tagg = 0.7
e_clist = np.arange(0.005,0.05,0.0045) # c-tag efficiency

while True: 
    inputec = input('e_c single or range (s or r)?:') 
    if inputec == 'r' or inputec == 'R':
        e_c = e_clist
        break
    elif inputec =='s' or inputec == 'S':
        e_c = 0.01
#        e_c = float(input('e_c value prefered:'))
        break
    else:
        inputec = inputec
e_4jet = 0.841 # charged Higgs after 4-jet selection
e_mass = 0.718 / 0.841 #charged Higgs after mass-selection
e_antiww = 0.536 / 0.578 # charged Higgs after anti-ww background
e_l = 0.01 #light jet tag efficiency
e_isignal = 0.9  # invariant mass cut on signal
e_iback = 0.1 # invariant mass cut on background
e_ibacklist = np.arange(0.1,1.08,0.09)
print(len(e_ibacklist),e_ibacklist)
epsilon = e_mass * e_4jet * e_isignal # total epsilon
############
# CKM elements
vcs = 0.97
vcb = 0.04
###########
mhch = np.arange(80.0,91.0 ,1.0)# charged Higgs ranged values
print('charH_mass:',mhch)
costhetaw = mw / mz                     # cos (weinberg angle : thetaw)
sinthetaw = math.sqrt(1 - (costhetaw)**2)      # sin(weinberg angle: thetaw)
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
#def backgroundtagging():# ww>cscs 4jet DELPHI
    #wwcscs with e_c , or 3389.6
#    result = 3389.6 * 1855.5 / 2506.2 / 4  * e_c**2
    # Z>b b-bar after b-tag .0.15 from z>bbbar oit of Z-decay ; e_iback = 0.1 invariant mass cut
#    qqbar = 964.6 * 0.15 * e_b**2 * e_iback
#    return result + qqbar
def backgroundtagging():# OPAL 4jets tagged background
    # total 1117.8 events. 90% is ww.
    ww = 1117.8 * 0.9
    ww_cscs = ww / 4 * e_c**2
    qq_bar = 1117.8 * 0.1 * 0.15 * e_b**2  #z>bb_bar
    return  ww_cscs + qq_bar
#print('OPAL',backgroundtagging())
def backtag_invarmasscut():# OPAL 4jets tagged background with invariant mass cut
    # total 1117.8 events. 90% is ww.
    ww = 1117.8 * 0.9
    ww_cscs = ww / 4 * e_c**2
    qq_bar = 1117.8 * 0.1 * 0.15 * e_b**2 * e_ibacklist  #z>bb_bar
    return  ww_cscs + qq_bar
#print('invariantmass cut 4jet tag',backtag_invarmasscut(),len(backtag_invarmasscut()))
def backgroundtagging2():# 2jet tagged
    e_l = 0.01
    background = 300 / 2.0  * (e_c * (1 - e_l) + e_l * (1 - e_c))
    return background
def backgroundnotagging2():# 2jet untagged
    background = 300 / 2.0   
    return background
#def backgroundnotagging():# ww>4jet no specific
    #wwcscs with e_c , or 3389.6
#    result_1 = 3389.6 * 1855.5 / 2506.2 
    # Z>b b-bar after b-tag .0.15 from z>bbbar oit of Z-decay ; e_iback : 0.1 invariant mass cut
#    qqbar_1 = 964.6 * e_iback
#    return result_1 + qqbar_1
def backgroundnotagging():# OPAL 4jets untagged background
    # total 1117.8 events. 90% is ww.
    ww = 1117.8 * 0.9
    ww_cscs = ww / 4 
    qq_bar = 1117.8 * 0.1 * 0.15   #z>bb_bar
    return  ww_cscs + qq_bar
#print('WW>cscs and qq_bar background:',backgroundtagging())
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
    for j in mhch:#char_H mass
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
    for k in np.arange(0,len(np.array(S))*len(mhch)):
        #expand luminosity to fit length of chargedHiggs mass length
        Lumin_h_array = np.array(Lumin_h *len(mhch))
        #Single(cross-section * Luminosity)
        Product_HH = np.array(cross_section_eeHH())[k] * Lumin_h_array [k]  
        Summ_HH.append(Product_HH) # List_(cross-section * Luminosity)
    #split whole List_(cross-section * Luminosity) to different charged Higgs mass 
    #region in order to calculate signal events for single chargedHiggs mass
    eecharH_seperate = np.split(np.array(Summ_HH),len(mhch))
    #calculate single charged_higgs mass total signals from SIGMA(cross-section * luminosity)
#    print('SIGMA(cross-section * luminosity) events for each charged Higgs mass:')
    for i in np.arange(0,len(mhch)):
    #SIGMA(cross-section * luminosity) of each charged Higgs mass produced signals
        eecharH = sum(eecharH_seperate[i]) 
        eventHH.append(eecharH)
#        print(i,'charH_mass:',mhch[i],'GeV',eecharH,'events',eventHH)
    return np.array(eventHH)
print('SIGNAL after H+H- > cbcb',eeHH_event() * 0.8**2 ,'events')
print('signal/sqrt(background) FOR cbcb decay',eeHH_event()[0] * 0.8**2 * 0.7**2 / \
      np.sqrt(backgroundtagging() * 0.1**2))
print('signal/sqrt(background) FOR no-tagging ',eeHH_event()[0] * 1.0**2* 0.536 / \
      np.sqrt(1855.5))
print('|Background for 2jetTN',backgroundtagging2())
# epsilon_signal = 0.536 * 71.8 / 57.8 = 0.6658269896193771
######################################################################
#with b signal of charged Higgs event:
def eeHHcbcb_2bsignal(x,y): #2real b 0 fake b
    chunk1 = np.array(x) * np.array(y) * \
    e_b**2 * (1 - e_c)**2# before compare with table 
    return np.array(chunk1) # * epsilon  # after 3 selections chosen
#print('||||| Singal of char_H_cbcb after 2b-tagging:',eeHHcbcb_2bsignal(),'events')
######################################################################
def eeHHcbcb_1bsignal(x,y):#1real b 1 fake b
    chunk1 =  np.array(x) * np.array(y) * 4.0 * \
    e_b * e_c * (1 - e_b) * (1 - e_c)
      # before compare with table 
    return np.array(chunk1) #* epsilon # * e_antiww # after 3 selections chosen
#print('||||| Singal of char_H_cbcb after 1b-tagging:',eeHHcbcb_1bsignal())
######################################################################
def eeHHcbcb_0bsignal(x,y):#0real b 2 fake b
    chunk1 =  np.array(x) * np.array(y) * e_c**2 * (1 - e_b)**2# before compare with table 
    return chunk1 # * epsilon # * e_antiww # after 3 selections chosen
#print('||||| Singal of char_H_cbcb after 0b-tagging:',eeHHcbcb_0bsignal())
######################################################################
def eeHHcbcs_2bsignal(x,y):#2real b 0 fake b
    return  2.0 * np.array(x) * np.array(y) * 0.0 #* epsilon #* e_antiww # after 3 selections chosen # 2.0 for cbcs and cscb 
#print('||||| Singal of char_H after 2b-tagging:',eeHHcbcs_2bsignal())
######################################################################
def eeHHcbcs_1bsignal(x,y):#1real b 1 fake b
    e_l = 0.01
    chunk2 = e_b * e_c * (1 - e_c) * (1 - e_l) + e_b * e_l * (1 - e_c)**2
    chunk1 =  2.0 * np.array(x) * np.array(y) * chunk2  # before compare with table  # 2.0 for cbcs and cscb 
    return np.array(chunk1) #* epsilon #* e_antiww # after 3 selections chosen
#print('||||| Singal of char_H_cbcs after 1b-tagging:',eeHHcbcs_1bsignal())
######################################################################
def eeHHcbcs_0bsignal(x,y):#0real b 2 fake b
    e_l = 0.01
    chunk2 =  2 * e_c * e_l * (1 - e_b) * (1 - e_c) + e_c**2 * (1 - e_b) * (1 - e_l)
    chunk1 =  2.0 * np.array(x) * np.array(y) * chunk2  # before compare with table  # 2.0 for cbcs and cscb 
    return np.array(chunk1) # * epsilon #* e_antiww # after 3 selections chosen
#print('||||| Singal of char_H_cbcs after 0b-tagging:',eeHHcbcs_0bsignal())
    
######################################################################
def eeHHcscs_0bsignal(x,y): #0real b 2 fake b
    e_l = 0.01
    chunk2 = 2 * e_c * e_l * (1 - e_c) * (1 - e_l) + e_c**2 * (1 - e_l)**2 + \
    e_l**2 * (1 - e_c)**2
    chunk1 =  np.array(x) * np.array(y) * chunk2  # before compare with table 
    return np.array(chunk1) #* epsilon #* e_antiww  after 3 selections chosen
######################################################################
######################################################################
def eeHHcbtn_1bsignal(x,y):
#    brcb = 0.8
#    brtn = 0.2
    chunk1 = 2.0 * np.array(x) * np.array(y) *  e_b * (1 - e_c) # before compare with table  # permutations  
    return np.array(chunk1) #* epsilon #* e_antiww  after 3 selections chosen
#print('||||| Singal of char_H_cbtn after 1b-tagging:',eeHHcbtn_1bsignal())
######################################################################
def eeHHcbtn_0bsignal(x,y):
    chunk1 =  2.0 * np.array(x) * np.array(y) * e_b * (1 - e_c) # 2.0 for cbtn and tncb permutations
     # before compare with table 
    return np.array(chunk1) #* epsilon #* e_antiww  after 3 selections chosen
#print('||||| Singal of char_H_cbtn after 0b-tagging:',eeHHcbtn_0bsignal())
######################################################################
def eeHHcstn_1bsignal(x,y):
#    brcs = 0.8
#    brtn = 0.2
    chunk1 = 2.0 * np.array(x) * np.array(y) * 0.0  # before compare with table # 2.0 for cstn and tncs permutations
    return np.array(chunk1)#* epsilon# * e_antiww  after 3 selections chosen 
#print('||||| Singal of char_H_cstn after 1b-tagging:',eeHHcstn_1bsignal())
######################################################################
def eeHHcstn_0bsignal(x,y):
#    brcs = 0.3
#    brtn = 0.4
    e_l = 0.01 # e_s
    chunk2 = e_c * (1 - e_l) + e_l * (1 - e_c)
    chunk1 =  2.0 * np.array(x) * np.array(y) * chunk2  # before compare with table # 2.0 for cstn and tncs 
    return np.array(chunk1) #* epsilon# * e_antiww  after 3 selections chosen 
################################
###########################################################################
plt.figure()
plt.plot(np.sqrt(np.array(S)), cross_section_eeww()/ (2.56819 * 10 **(-9)))
plt.ylabel('Cross-section in pb^-1')
plt.xlabel('CoM energy Squared root of S')
plt.savefig('picture.png')
plt.show()
####################################################################################
################################################################################### 
def boson_decay_exist():  #check bosonic decay of charged Higgs exist or not and get lamaw value
    for i in mhch:                                                                                                                                                                                              
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
    return matrixmul(x,y)
#plt.figure()
#Contoursigma1 = plt.contour(e_clist,e_blist, \
#    np.resize(sig(eeHH_event() * 0.8**2 * e_blist**2,backgroundtagging()* e_clist**2),\
#    len(sig(eeHH_event() * 0.8**2 * e_blist**2,backgroundtagging() * e_clist**2))).reshape(len(e_blist),len(e_clist))\
#)# ,levels = np.arange(1.0,6.0,1.0))
#plt.title('Significance of H+H->cbcb / sqrt(background) with e_b and e_c')
#plt.xlabel('e_c')# x-axis label
#plt.ylabel('e_b')# y-axis label
#plt.colorbar(Contoursigma1)

#plt.figure()
#Contoursigma2 = plt.contour(e_clist,e_blist, \
    #np.resize(sig(eeHH_event() * 0.8**2 * matrixmul(e_clist**2 , (1 - e_blist)**2),backgroundtagging()* e_clist**2),\
    #int(len(sig(eeHH_event() * 0.8**2  * matrixmul(e_clist**2 , (1 - e_blist)**2),backgroundtagging()* e_clist**2))/len(e_clist))).\
    #reshape(len(e_blist),len(e_clist)))# ,levels = np.arange(1.0,6.0,1.0))
#plt.title('Significance of H+H->cbcb / sqrt(background) with e_b and e_c')
#plt.xlabel('e_c')# x-axis label
#plt.ylabel('e_b')# y-axis label
#plt.colorbar(Contoursigma2)
#plt.figure()
#Contoursignal2 = plt.contour(e_c,e_b, \
#           np.resize(eeHHcbcb_1bsignal(),len(eeHHcbcb_1bsignal())).reshape(len(e_b),len(e_c)))
#plt.colorbar(Contoursignal2)

##############################################################
def massH_ec_plane4jet(x,y):
    print('wlw',mhch,eeHH_event())
#     for  c in e_c:
#        print(eeHHcbcb_2bsignal(0.8,0.8))
#    for j in mhch:#char_H mass
    tagging_4jet = (eeHHcbcb_2bsignal(x,x) + eeHHcbcb_1bsignal(x,x) + \
                   eeHHcbcs_2bsignal(x,y) + eeHHcbcs_1bsignal(x,y) + \
                   eeHHcbcb_0bsignal(x,x) + eeHHcbcs_0bsignal(x,y) + \
                   eeHHcscs_0bsignal(y,y)) * epsilon
#    eeHH_eventarry,tagging_4jetarry = np.meshgrid(eeHH_event(),tagging_4jet)
#    print(len(eeHH_event()), len(tagging_4jet),len(eeHH_eventarry), len(tagging_4jetarry))
#    for j in eeHH_event():
#        print (j)
    print('tagging_4jet',tagging_4jet,len(tagging_4jet))
    print('backgroundtagging()',backgroundtagging(),len(backgroundtagging()))
    print('eeHH_event()',eeHH_event(),len(eeHH_event()))
    print('sig',sig(eeHH_event() , tagging_4jet / np.sqrt(backgroundtagging() )),\
          len(sig(eeHH_event() , tagging_4jet / np.sqrt(backgroundtagging() ))) )
    plt.figure()
    signal4jet_tag = plt.contourf(mhch,e_c,\
                   np.resize(sig(eeHH_event() , tagging_4jet / np.sqrt(backgroundtagging() )),\
                          len(sig(eeHH_event() , tagging_4jet / np.sqrt(backgroundtagging() )))).\
               reshape(len(e_c),len(mhch)),\
#                levels = np.arange(0.0, 8.0,1.0), \
               colors = ['black','royalblue','purple','darkgreen','brown','red','gray'])
    plt.title('S/$\sqrt{B}$ of $H^{\pm}$ 4jet_tag with max BR ')
    plt.xlabel('mhch')# x-axis label
    plt.ylabel('e_c')# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.colorbar(signal4jet_tag)
    plt.savefig('sig_4jetecmhch.png')
    plt.show()
    plt.close()
##########################################################################    
##########################################################################    
def massH_soverb4jetag(x,y):#background not tagged with e_b and e_c
    fourjet = (eeHHcbcb_2bsignal(x,x) + eeHHcbcb_1bsignal(x,x) + \
                   eeHHcbcs_2bsignal(x,y) + eeHHcbcs_1bsignal(x,y) + \
                   eeHHcbcb_0bsignal(x,x) + eeHHcbcs_0bsignal(x,y) + \
                   eeHHcscs_0bsignal(y,y)) * epsilon 
    
    plt.plot(mhch,eeHH_event() * fourjet / np.sqrt(backgroundtagging()) )
    plt.title('Relation between mhch and S/$\sqrt{B}$ in 4jettag')
    plt.xlabel('mhch')# x-axis label
    plt.ylabel('S/$\sqrt{B}$')# y-axis label
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('sig4jettag_mhch.png')
    plt.show()
    plt.close()
##########################################################################
def massH_soverb4jetnotag(x,y):#background not tagged with e_b and e_c
    fourjet = (np.array(x) * np.array(y)) * epsilon
    
    plt.plot(mhch,eeHH_event() * fourjet / np.sqrt(backgroundnotagging()) )
    plt.title('Relation between mhch and S/$\sqrt{B}$ in 4jetnotag')
    plt.xlabel('mhch')# x-axis label
    plt.ylabel('S/$\sqrt{B}$')# y-axis label
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('sig4jetnotag_mhch.png')
    plt.show()
    plt.close()
##################################################################
def massH_ec_plane2jet(x,y):
#    BR($H^{\pm} \longrightarrow $ cb * tn) with tagging efficiencies
    tagging_2jet = eeHHcbtn_1bsignal(x,y) + eeHHcbtn_0bsignal(x,y) + \
                                eeHHcstn_1bsignal(x,y) + eeHHcstn_0bsignal(x,y)
                                 
    plt.figure()
    signal2jet_tag = plt.contourf(mhch,e_c, \
        np.resize(sig(eeHH_event() , tagging_2jet * 0.3 / np.sqrt(backgroundtagging2())) ,\
              len(sig(eeHH_event() , tagging_2jet * 0.3 / np.sqrt(backgroundtagging2())) )).\
        reshape(len(e_c),len(mhch)),cmap = 'brg')# 0.3 signal selection efficiency
    plt.colorbar(signal2jet_tag)
    plt.title('S/$\sqrt{B}$ 2jet_tag with max BR')#plot title
    plt.xlabel('mhch')# x-axis label
    plt.ylabel('e_c')# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('sig_2jetecmhch.png')
    plt.show()
    plt.close()
##################################################################
##################################################################
def massH_soverb2jetag(x,y):#background not tagged with e_b and e_c
    twojet = (eeHHcbtn_1bsignal(x,y) + eeHHcbtn_0bsignal(x,y) + \
             eeHHcstn_1bsignal(x,y) + eeHHcstn_0bsignal(x,y) ) * 0.3
    
    plt.plot(mhch,eeHH_event() * twojet / np.sqrt(backgroundtagging2()) )
    plt.title('Relation between mhch and S/$\sqrt{B}$ in 2jettag')
    plt.xlabel('mhch')# x-axis label
    plt.ylabel('S/$\sqrt{B}$')# y-axis label
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('sig2jettag_mhch.png')
    plt.show()
    plt.close()
##################################################################
def massH_soverb2jetnotag(x,y): #background not tagged with e_b and e_c
#    BR($H^{\pm} \longrightarrow $ cb * tn) with tagging efficiencies
    twojet = np.array(x) * np.array(y)
                                 
    plt.figure()
    plt.plot(mhch,eeHH_event() * \
                               twojet / np.sqrt(backgroundnotagging2()))# 0.3 signal selection efficiency
#    plt.colorbar(signal2jetnotag)
    plt.title('Relation between mhch and S/$\sqrt{B}$ in 2jetnotag')#plot title
    plt.xlabel('mhch')# x-axis label
    plt.ylabel('S/$\sqrt{B}$')# y-axis label
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('sig2jetnotag_mhch.png')
    plt.show()
    plt.close()
###################################################################
# 4jet tagging with plane in [e_c ,e_b] for mhch from 80 to 90 GeV
def ec_eb_plane4jet(x,y):
    global e_b,e_c
    for n in np.arange(0,len(mhch)): # for each time, eeHH_event()[n] gives specific charH event
        totalcbcb_2bsig =[]
        totalcbcs_2sig =[]
        totalcbcb_1sig =[]
        totalcbcs_1sig =[]
        totalcbcb_0sig =[]
        totalcbcs_0sig =[]
        totalcscs_0sig =[]
        background4tag = []
        for b in e_blist:
            e_b = b
            for c in e_clist:
                e_c = c
                totalcbcb_2bsig.append(float(eeHHcbcb_2bsignal(x,x)))
                totalcbcs_2sig.append(float(eeHHcbcs_2bsignal(x,y)))
                totalcbcb_1sig.append(float(eeHHcbcb_1bsignal(x,x)))
                totalcbcs_1sig.append(float(eeHHcbcs_1bsignal(x,y)))
                totalcbcb_0sig.append(float(eeHHcbcb_0bsignal(x,x)))
                totalcbcs_0sig.append(float(eeHHcbcs_0bsignal(x,y)))
                totalcscs_0sig.append(float(eeHHcscs_0bsignal(y,y)))
                background4tag.append(backgroundtagging())
        tagging_4jet = np.array(totalcbcb_2bsig) + np.array(totalcbcs_2sig) + np.array(totalcbcb_1sig) +\
            np.array(totalcbcs_1sig) + np.array(totalcbcb_0sig) + np.array(totalcbcs_0sig) + np.array(totalcscs_0sig)
#        print(tagging_4jet,type(tagging_4jet),len(background4tag),len(tagging_4jet))
        one = plt.contourf(e_clist,e_blist,\
            np.resize(eeHH_event()[n] * tagging_4jet * epsilon / np.sqrt(np.array(background4tag)) ,\
                  len(eeHH_event()[n] * tagging_4jet * epsilon / np.sqrt(np.array(background4tag)) )).\
                reshape(len(e_blist),len(e_clist)),\
#                 levels = np.arange(0.0, 8.0,1.0), \
                colors = ['black','royalblue','purple','darkgreen','brown','red','gray'])
        plt.colorbar(one)
        plt.title('S/$\sqrt{B}$ 4jet_tag Max BR with charH mass: '+ str(mhch[n]))#plot title
        plt.xlabel('e_c')# x-axis label
        plt.ylabel('e_b')# y-axis label
        plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
        plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
        plt.savefig('sig_4jetebec'+ str(mhch[n]) + 'GeV.png')
        plt.show()
        plt.close()
########################################
# 2jet tagging with plane in [e_c ,e_b] for mhch from 80 to 90 GeV
def ec_eb_plane2jet(x,y):
    global e_b,e_c
    for n in np.arange(0,len(mhch)):# for each time, eeHH_event()[n] gives specific charH event
        totalcbtn_1bsig =[]
        totalcbtn_0bsig =[]
        totalcstn_1bsig =[]
        totalcstn_0bsig =[]
        background2tag = []
        for b in e_blist:
            e_b = b
            for c in e_clist:
                e_c = c
                totalcbtn_1bsig.append(float(eeHHcbtn_1bsignal(x,y)))
                totalcbtn_0bsig.append(float(eeHHcbtn_0bsignal(x,y)))
                totalcstn_1bsig.append(float(eeHHcstn_1bsignal(x,y)))
                totalcstn_0bsig.append(float(eeHHcstn_0bsignal(x,y)))
                background2tag.append(backgroundtagging2())
#                print('ebeccbtn1b',e_b,e_c,eeHHcbtn_1bsignal(0.65,0.25))
        tagging_2jet = np.array(totalcbtn_1bsig) + np.array(totalcbtn_0bsig) + np.array(totalcstn_1bsig) +\
            np.array(totalcstn_0bsig)
#        print(tagging_2jet,type(tagging_2jet),len(background2tag),len(tagging_2jet))
        one1 = plt.contourf(e_clist,e_blist,\
            np.resize(eeHH_event()[n] * tagging_2jet * 0.3 / np.sqrt(np.array(background2tag)) ,\
                  len(eeHH_event()[n] * tagging_2jet * 0.3 / np.sqrt(np.array(background2tag)) )).\
             reshape(len(e_blist),len(e_clist)),\
#             levels = np.arange(0.0, 8.0,1.0), \
             colors = ['black','royalblue','purple','darkgreen','brown','red','gray'])
        plt.colorbar(one1)
        plt.title('S/$\sqrt{B}$ 2jet_tag Max BR with charH mass: '+ str(mhch[n]))#plot title
        plt.xlabel('e_c')# x-axis label
        plt.ylabel('e_b')# y-axis label
        plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
        plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
        plt.savefig('sig_2jetebec'+ str(mhch[n]) + 'GeV.png')
        plt.show()
        plt.close()
###################################################################
def invariantmsscut_ec(x,y):#4jet case specific
    global e_c
    for n in np.arange(0,len(mhch)):# for each time, eeHH_event()[n] gives specific charH event
         totalcbcb_2bsig =[]
         totalcbcs_2sig =[]
         totalcbcb_1sig =[]
         totalcbcs_1sig =[]
         totalcbcb_0sig =[]
         totalcbcs_0sig =[]
         totalcscs_0sig =[]
         background4tag = []
         for c in e_clist:
                e_c = c
                totalcbcb_2bsig.append(float(eeHHcbcb_2bsignal(x,x)))
                totalcbcs_2sig.append(float(eeHHcbcs_2bsignal(x,y)))
                totalcbcb_1sig.append(float(eeHHcbcb_1bsignal(x,x)))
                totalcbcs_1sig.append(float(eeHHcbcs_1bsignal(x,y)))
                totalcbcb_0sig.append(float(eeHHcbcb_0bsignal(x,x)))
                totalcbcs_0sig.append(float(eeHHcbcs_0bsignal(x,y)))
                totalcscs_0sig.append(float(eeHHcscs_0bsignal(y,y)))
                background4tag.append(backgroundtagging())
         e_ix,e_cy = np.meshgrid(e_ibacklist,e_clist)
         tagging_4jet = np.array(totalcbcb_2bsig) + np.array(totalcbcs_2sig) + np.array(totalcbcb_1sig) +\
np.array(totalcbcs_1sig) + np.array(totalcbcb_0sig) + np.array(totalcbcs_0sig) + np.array(totalcscs_0sig)
         s4jettagmasscut = plt.contourf(e_ibacklist,e_clist,\
np.resize(sig(eeHH_event()[n] * tagging_4jet, 1.0 / np.sqrt(backtag_invarmasscut() ) ),\
len(sig(eeHH_event()[n] * tagging_4jet, 1.0 / np.sqrt(backtag_invarmasscut() ) ))).\
reshape(len(e_clist),len(e_ibacklist)),\
     # levels = np.arange(0.0, 8.0,1.0), 
colors = ['black','royalblue','purple','darkgreen','brown','red','gray'])
         plt.title('S/$\sqrt{B}$ of $H^{\pm}$ 4jet_tag under mass cut '+ str(mhch[n]))
         plt.xlabel('invarmasscut')# x-axis label
         plt.ylabel('e_c')# y-axis label
         plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
         plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
         plt.colorbar(s4jettagmasscut)
         plt.savefig('sig4jetecmasscut'+ str(mhch[n]) + 'GeV.png')
         plt.show()
         plt.close()
###################################################################
def mhch_invariantmsscut(x,y):#4jet case specific  
     tagging_4jet = (eeHHcbcb_2bsignal(x,x) + eeHHcbcb_1bsignal(x,x) + \
                   eeHHcbcs_2bsignal(x,y) + eeHHcbcs_1bsignal(x,y) + \
                   eeHHcbcb_0bsignal(x,x) + eeHHcbcs_0bsignal(x,y) + \
                   eeHHcscs_0bsignal(y,y)) * epsilon    
     print('wwwwwwwwwwww',sig(eeHH_event() ,tagging_4jet / np.sqrt(backtag_invarmasscut() ) ) )
     s4jettagmasscut = plt.contourf(mhch,e_ibacklist,\
np.resize(sig(eeHH_event() , tagging_4jet / np.sqrt(backtag_invarmasscut() )),\
len(sig(eeHH_event() , tagging_4jet / np.sqrt(backtag_invarmasscut() )))).\
reshape(len(e_ibacklist),len(mhch)),\
                levels = np.arange(0.0,31.0,5.0), \
colors = ['black','royalblue','purple','darkgreen','brown','red','gray'])
     plt.title('S/$\sqrt{B}$ of $H^{\pm}$ 4jet_tag with max BR ')
     plt.xlabel('mhch')# x-axis label
     plt.ylabel('invarmasscut')# y-axis label
     plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
     plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
     plt.colorbar(s4jettagmasscut)
     plt.savefig('sig4jetmhch_masscut.png')
     plt.show()
     plt.close()                 
###################################################################
print('e-b,e-c',e_b,e_c)
#print(type(sig(eeHH_event(), tagging_4jet / np.sqrt(backtag_invarmasscut() ))),\
#      sig(eeHH_event(), tagging_4jet / np.sqrt(backtag_invarmasscut() )))
if type(e_c) == type(e_clist):
    massH_ec_plane4jet(0.8,0.19)
    massH_ec_plane2jet(0.65,0.35)
    ec_eb_plane4jet(0.8,0.19)
    ec_eb_plane2jet(0.65,0.35)
#    ec_eb_plane2jetnotag(0.40,0.33)
    invariantmsscut_ec(0.8,0.19)
else:
    mhch_invariantmsscut(0.8,0.19)
    massH_soverb2jetnotag(0.65,0.35)#0.40,0.33
    massH_soverb2jetag(0.65,0.35)
    massH_soverb4jetnotag(0.8,0.19)
    massH_soverb4jetag(0.8,0.19)
# mhch[n] : the n th charH mass, eeHH_events()[n] : the n th charH events produced
for n in np.arange(0,len(mhch)):
    print('#',n,mhch[n],eeHH_event()[n])
    print('eb,ec',e_b,e_c)
########################
#ms32g13@soton.ac.uk
