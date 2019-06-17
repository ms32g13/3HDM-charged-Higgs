#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 12:43:28 2019

@author: muyuansong0419
"""

###     Program which calculates the Production of the charged              
###     Higgs from e+e- collider CEPC. Followed by H+ decay to to cb,cs,tn plots
###

import math
import random
import numpy as np
import matplotlib.pyplot as plt
from invariant import *
import eehh as fl

mhch_list_CEPC = np.arange(80.0,122.0,4.0) # charged Higgs for CEPC search
#########################################################
# CEPC with integrated Luminosity 5000 fb^-1 to 5000 fb^-1 with COM : sqrt(s) = 240 GeV
def cross_section_eeHH_CEPC():# Total crosssection of e+e- > H+H- in COM = 240 GeV
    ListeechaHH = [] 
    S = 240**2 #COM^2
    for j in mhch_list_CEPC:
        BETA = math.sqrt(1.0 - 4 * (j)**2 / S)
        eechaHH = (1.0 / 4) * (4.0 * PI * alpha_electroweak **2 ) / (3.0 * S) * BETA**3 * fl.F_HIGGS(S)
        sigma_eeHH = eechaHH / (2.56819 * 10 **(-9)) #convert pb^(-1) to GeV
        ListeechaHH.append(sigma_eeHH)
    return ListeechaHH
#0.001ab^-1 = 1fb^-1 = 1000 pb^-1
def eeHH_event_CEPC(i):# i will = pb^-1 luminoscity
    event = np.array(cross_section_eeHH_CEPC()) * i
    for n in np.arange(0,len(mhch_list_CEPC)):
        print('event',event[n],mhch_list_CEPC[n])
    return event
#################################################
def massH_ec_plane4jet(x,y):
    tagging_4jet = (fl.eeHHcbcb_2bsignal(x,x) + fl.eeHHcbcb_1bsignal(x,x) + \
                   fl.eeHHcbcs_2bsignal(x,y) + fl.eeHHcbcs_1bsignal(x,y) + \
                   fl.eeHHcbcb_0bsignal(x,x) + fl.eeHHcbcs_0bsignal(x,y) + \
                   fl.eeHHcscs_0bsignal(y,y)) * fl.epsilon
    sig42b = fl.sig(CEPCevent, tagging_4jet/ np.sqrt(fl.backgroundtagging() * 1600 ) ) 
    print('tagging_4jet',tagging_4jet,len(tagging_4jet))
    print('CEPCevent',CEPCevent,len(mhch_list_CEPC), )
    print('backgroundtagging()',fl.backgroundtagging() * 1600,len(fl.backgroundtagging()))
    print('sig',sig42b,len(sig42b),len(fl.e_c),len(mhch_list_CEPC),len(fl.e_c) )
    plt.figure()
    signal4jet_tag = plt.contourf(mhch_list_CEPC,fl.e_c,\
                   np.resize(sig42b ,\
                          len(sig42b )).\
               reshape(len(fl.e_c),len(mhch_list_CEPC)),\
#                levels = np.arange(0.0, 8.0,1.0), \
               colors = ['black','royalblue','purple','darkgreen','brown','red','gray'])
    plt.title('S/$\sqrt{B}$ of $H^{\pm}$ 4jet2b'+ ':cb = '+ str(x))
    plt.xlabel('$M_{H^{\pm}}$')# x-axis label
    plt.ylabel('$e_c$')# y-axis label
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.colorbar(signal4jet_tag)
    plt.savefig('cepc4j2becmhch'+':cb:'+ str(x) +'.png')
    plt.show()
    plt.close()
#############################################
def massH_ec_planeone4jet(x,y):
    tagging_one4jet = ( fl.real_b_cbcb(x,x) + fl.fake_b_cbcb(x,x) + \
                        fl.real_b_cbcs(x,y) + fl.fake_b_cbcs(x,y) + \
                        fl.real_b_cscs(y,y) + fl.fake_b_cscs(y,y) ) * fl.epsilon
    sig41b = fl.sig(CEPCevent,tagging_one4jet/ np.sqrt(fl.backgroundtag4j1b() * 1600 ) )        
    plt.figure()
    sig4jet1tag = plt.contourf(mhch_list_CEPC,fl.e_c,\
                   np.resize(sig41b ,\
                          len(sig41b )).\
               reshape(len(fl.e_c),len(mhch_list_CEPC)),\
#                levels = np.arange(0.0, 8.0,1.0), \
               colors = ['black','royalblue','purple','darkgreen','brown','red','gray'])
    plt.title('S/$\sqrt{B}$ of $H^{\pm}$ 4jet1b'+ ':cb = '+ str(x))
    plt.xlabel('$M_{H^{\pm}}$')# x-axis label
    plt.ylabel('$e_c$')# y-axis label
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.colorbar(sig4jet1tag)
    plt.savefig('cepc4j1becmhch'+':cb:'+ str(x) +'.png')
    plt.show()
    plt.close()
####################################################
def massH_ec_plane2jetag(x,y,z): #2jtagged
    twojet = (fl.eeHHcbtn_1bsignal(x,y) + fl.eeHHcbtn_0bsignal(x,y) + \
             fl.eeHHcstn_1bsignal(y,z) + fl.eeHHcstn_0bsignal(y,z) ) * fl.selection_2j
    sig21b = fl.sig(CEPCevent,twojet/ np.sqrt(fl.backgroundtagging2() * 1600 ) )   
    sig2jet1tag = plt.contourf(mhch_list_CEPC,fl.e_c,\
                   np.resize(sig21b ,\
                          len(sig21b )).\
               reshape(len(fl.e_c),len(mhch_list_CEPC)),\
#                levels = np.arange(0.0, 8.0,1.0), \
               colors = ['black','royalblue','purple','darkgreen','brown','red','gray'])
    plt.title('S/$\sqrt{B}$ of $H^{\pm}$2jet1b'+ ':cb = '+ str(x) )
    plt.xlabel('$M_{H^{\pm}}$')# x-axis label
    plt.ylabel('$e_c$')# y-axis label
    plt.colorbar(sig2jet1tag)
    plt.savefig('cepc2j1becmhch'+':cb:'+ str(x) +'.png')
    plt.show()
    plt.close()
#######################################
CEPCevent = eeHH_event_CEPC(float(input('Luminoscity required (1000 * 10^3 pb^-1):')))#  pb^-1 values
if type(fl.e_c) == type(fl.e_clist):
    massH_ec_plane4jet(0.02,0.5)#
    massH_ec_planeone4jet(0.02,0.5)
    massH_ec_plane2jetag(0.02,0.5,1 - 0.02 - 0.5)
else:
    print('CEPCevent',CEPCevent, mhch_list_CEPC)