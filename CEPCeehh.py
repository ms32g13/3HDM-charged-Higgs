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
import vegas
import gvar as gv
import itertools
from invariant import *
import eehh as fl

mhch_list_CEPC = np.arange(80.0,130.0,10.0) # charged Higgs for CEPC search
#########################################################
# CEPC with integrated Luminosity 5000 fb^-1 to 5000 fb^-1 with COM : sqrt(s) = 240 GeV
def cross_section_eeHH_CEPC():# Total crosssection of e+e- > H+H- in COM = 240 GeV
    ListeechaHH = [] 
    S = 241**2 #COM^2
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
CEPCevent = eeHH_event_CEPC(float(input('Luminoscity required (1000 * 10^3 pb^-1):')))#  pb^-1 values
print(CEPCevent)


    