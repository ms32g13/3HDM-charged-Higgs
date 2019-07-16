#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 14:58:53 2019
@author: muyuansong0419
"""

import numpy as np
import matplotlib.pyplot as plt
import gammaeffective as gaeff
from scipy.integrate import quad
from scipy import special as sp
from invariant import *
from matplotlib.ticker import MultipleLocator
#from eehh import mhch
#from exercise import xyfun_list,yfun_list,xyfun,yfun
#import exercise as ex
# I make two arrays to pretend solution of XY* and Y
#ll1 = np.array(xyfun) #j = XY^*
#ll2 = np.array(yfun)#i = Y
#############
############
#QCD running coupling constant (alp) at energy scale MH. Relevant                                                                                                                                                                        
#for QCD corrections to Higgs decay widths.
################ LO for strong coupling constant
def LOalpha_s(i): #alpha_s(mu) at LO
#    print(i,zz,alpmz,mz)
    beta_0 = 23/3
    beta_1 = 0
    v = 1 - beta_0 * (alpmz / (2.0 * PI)) * np.log(mz / i)
    ratio1 = np.log(v) / v
    return alpmz / v * (1 - (beta_1 / beta_0) * (alpmz / (4 * PI)) * ratio1)

################ NLO for strong coupling constant
def NLOalpha_s(i): #alpha_s(mu) at NLO
    beta_0 = 23 / 3
    beta_1 = 116 / 3
    v = 1 - beta_0 * (alpmz / (2.0 * PI)) * np.log(mz / i)
    ratio1 = np.log(v) / v
    return alpmz / v * (1 - (beta_1 / beta_0) * (alpmz / (4 * PI)) * ratio1)
print('mu_w scale at LO and NLO:',mw,LOalpha_s(mw),NLOalpha_s(mw))
print('mu_z scale at LO and NLO:',mz,LOalpha_s(mz),NLOalpha_s(mz))
print('mu_b scale at LO and NLO:',mb,LOalpha_s(mb),NLOalpha_s(mb))
print('mu_t scale at LO and NLO:',mt,LOalpha_s(mt),NLOalpha_s(mt))
#################running quark mass at scale mu_w under minimal subtract scheme
def run_quark_bar(q):# 
    c1 = np.log(q**2 / mw**2)
    c2 = LOalpha_s(mw) / PI
    return q * (1 + c2 * c1 - 4 / 3 * c2 )
#############
charHm_100 = np.array([100,300,500,1000]) #for MH+2 = 100GeV
charHm_large = np.array([300,500,1000,1200]) #for MH+3 > 100GeV
delta_cp = 0.3 #delta_cp fraction of Energy cut in CP-asymmetry 
zz = (run_quark_bar(mc) / run_quark_bar(mb) )**2  # mc^2 / mb^2
LOzz = (mc / mb)**2
xx = mt**2 / mw**2 # mt^2 /mw^2
NLOxx = run_quark_bar(mt)**2 / mw**2
PI = np.pi
# mt^2 / mh^2
def yy(mass):
    return mt**2 / mass**2
def NLOyy(mass):
    return run_quark_bar(mt)**2 / mass**2
print('zz',zz,1/zz)
a_i = np.array([14 /23, 16 /23, 6 /23, - 12/23,0.4086,-0.4230,-0.8994,0.1456])#{a_i}
h_i = np.array([626126/272277 , - 56281/51730, - 3/7, - 1/14, - 0.6494, - 0.0380,\
                - 0.0186, - 0.0057])#{h_i}
a2_i = np.array([14/23, 0.4086, - 0.4230, - 0.8994, 0.1456])#{a'_i}
h2_i = np.array([313063/363036, - 0.9135, 0.0873, - 0.0571, 0.0209])#{h'_i}
e_i = np.array([4661194/816831, - 8516/2217, 0.0,0.0,- 1.9043, - 0.1008,\
                0.1216,0.0183])#{e_i}
f_i = np.array([- 17.3023, 8.5027, 4.5508, 0.7519, 2.0040, 0.7476, - 0.5385, 0.0914])#{f_i}
k_i = np.array([9.9372, - 7.4878, 1.2688, - 0.2925, -2.2923, - 0.1461, 0.1239, 0.0812])#{k_i}
l_i = np.array([0.5784, - 0.3921, -0.1429, 0.0476, - 0.1275, 0.0317, 0.0078, - 0.0031])#{l_i}
y_i = np.array([0.0,0.0,- 1/3.0,- 4/9.0,- 20/3.0, - 80/9.0])# {y_i} 
z_i = np.array([0,0,1.0,- 1/6.0, 20, - 10/3.0])#{z_i}
##########################################################
# PHYSICAL REVIEW D, VOLUME 58,074004-7 EQUATION (36,39)
# PHASE- SPACE FUNCTION g(z)
g_z = 1 - 8 * zz + 8 * zz**3 - zz**4 - 12 * zz**2 * np.log(zz)
# QCD-RADIATION FUNCTION f(z)
f_z = (PI**2 - 31 / 4) * (1 - np.sqrt(zz))**2 + 3 / 2
#NON-Perturbative kronc-delta (GeV^2)
delta_NP_ga = - 0.5 / 2 - 9 * (- 0.12) / 2
delta_NP_c = - (- 0.12) / 2
delta_NP_SL = - 0.5 /2 + 3 * (- 0.12) / 2 * (1 - 4 * (1 - zz)**4 / g_z) 
#########################################################
print(gaeff.gamma0eff())#gamma_0_eff_ji matrix values
print(gaeff.gamma1eff())#gamma_1_eff_ji matrix values
print('run-bquark-at-mw-scale',run_quark_bar(4.18))
print('LOalpha_s(run_m_b)',LOalpha_s(run_quark_bar(mb)),\
      'LO a_s(mw) / LO a_s(run_mb)', LOalpha_s(mw) / LOalpha_s(run_quark_bar(mb)))
#########################################################################
# Function G(t) in E2 of PHYSICAL REVIEW D, Vol58, 074004-19
def G_(t):# t < 4 and t > 4
    if t < 4:
       return - 2 * np.arctan(np.sqrt(t /(4 - t) ) )**2
    else:
       log_value = np.log((np.sqrt(t) + np.sqrt(t - 4))/2)
       repart = - PI**2 / 2 + 2 * log_value**2
       impart =  - 2 * PI * log_value
       return complex(repart , impart)
def grand1(t): # f_22 integrand
    return (1 - zz * t)**2 * np.abs(G_(t) / t + 1 / 2)**2
def grand2(t): # f_27 integrand
    return (1 - zz * t) * (G_(t).real  + t / 2)
def grand3(t): # f_27 integrand imaginary
    return (1 - zz * t) * (G_(t).imag)
print('quad3',G_(1.0 / zz), zz,quad(grand3, 0, 1.0 / zz ), type(quad(grand3, 0, 1.0/ zz )))        
##########################################################################
#4*pi*SQRT(2); factor appears in partial width of 2 fermions
fac = 4.0 * np.sqrt(2.0) * PI 
############################
#############################
####Virtual Correction Functions r_i

##########################################################################
### Wilson coefficient at matching scale(mu_w)
####################### LO
def c0_7sm():
    x = xx
    chunk = - 8 * x**3 + 3 * x**2 + 12 * x - 7
    chunk1 = (18 * x**2 - 12 * x) * np.log(x)
    return x / 24 * ( (chunk + chunk1) / (x - 1)**4 )
def c0_8sm():
    x = xx
    chunk = - x**3 + 6 * x**2 - 3 * x - 2 - 6 * x * np.log(x)
    return x / 8 * ( chunk / (x - 1)**4 )
def c0_7yy(mass):
    y = yy(mass)
    chunk = - 8 * y**3 + 3 * y**2 + 12 * y - 7
    chunk1 = (18 * y**2 - 12 * y) * np.log(y)
    return y / 24 * ( (chunk + chunk1) / (y - 1)**4 ) * 1 / 3
def c0_8yy(mass):
    y = yy(mass)
    chunk = - y**3 + 6 * y**2 - 3 * y - 2 - 6 * y * np.log(y)
    return y / 8 * ( chunk / (y - 1)**4 ) * 1 / 3
def c0_7xy(mass):
    y = yy(mass)
    chunk = - 5 * y**2 + 8 * y  - 3
    chunk1 = (6 * y - 4) * np.log(y)
    return y / 12 * ( (chunk + chunk1) / (y - 1)**3)
def c0_8xy(mass):
    y = yy(mass)
    chunk = - y**2 + 4 * y - 3 - 2 * np.log(y)
    return y / 4 * ( chunk / (y - 1)**3 )
##########################
########################## NLO
def E_0():#NLO
    x = NLOxx
    chunk_1 = x * (x**2 + 11 * x - 18) / (12 * (x - 1)**3)
    chunk_2 = x**2 * (4 * x**2 - 16 * x + 15) * np.log(x) / (6 * (x - 1)**4)
    chunk_3 = - 2 / 3 * np.log(x) - 2 / 3
    return chunk_1 + chunk_2 + chunk_3
def w7_sm():#NLO
    x = NLOxx
    chunk_1 = (- 16 * x**4 - 122 * x**3 + 80 * x**2 - 8 * x) / (9 * (x - 1)**4)
    chunk_2 = sp.spence(1.0/x)
    chunk_3 = (6 * x**4 + 46 * x**3 - 28 * x**2) * (np.log(x))**2 / (3 * (x - 1)**5) 
    chunk_4 = (- 102 * x**5 - 588 * x**4 - 2262 * x**3 + 3244 * x**2 - 1364 * x + 208 ) * \
    np.log(x) / (81 * (x - 1)**5)
    chunk_5 = (1646 * x**4 + 12205 * x**3 - 10740 * x**2 + 2509 * x - 436) / (486 * (x - 1)**4)
    return chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5
def w8_sm():#NLO
    x = NLOxx
    chunk_1 = (- 4 * x**4 + 40 * x**3 + 41 * x**2 + x)/(6 * (x - 1)**4)
    chunk_2 = sp.spence(1.0 - 1.0/x)
    chunk_3 = (- 17 * x**3 - 31 * x**2) * (np.log(x))**2 / (2 * (x - 1)**5)
    chunk_4 = (- 210 * x**5 + 1086 * x**4 + 4893 * x**3 + 2857 * x**2 - 1994 * x + 280) * \
    np.log(x) / (216 * (x - 1)**5)
    chunk_5 = (737 * x**4 - 14102 * x**3 - 28209 * x**2 + 610 * x - 508) \
    / (1296 * (x - 1)**4)
    return chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5
def m7_sm():#NLO
    x = NLOxx
    chunk_1 = 82 * x**5 + 301 * x**4 + 703 * x**3 - 2197 * x**2 + 1319 * x - 208
    chunk_2 = - (162 * x**4 + 1242 * x**3 - 756 * x**2) * np.log(x)
    return (chunk_1 + chunk_2) / (81 * (x - 1)**5)
def m8_sm():#NLO
    x = NLOxx
    chunk_1 = 77 * x**5 - 475 * x**4 - 1111 * x**3 + 607 * x**2 + 1042 * x - 140
    chunk_2 = ( 918 * x**3 + 1674 * x**2) * np.log(x)
    return (chunk_1 + chunk_2) / (108 * (x - 1)**5)
def t7_sm():#NLO
    x = NLOxx
    chunk_1 = 47 * x**3 - 63 * x**2 + 9 * x + 7 
    chunk_2 = - (18 * x**3 + 30 * x**2 - 24 * x) * np.log(x)
    return x / 3 * ( (chunk_1 + chunk_2) / (x - 1)**5 )
def t8_sm():#NLO
    x = NLOxx
    chunk_1 = - x**3 - 9 * x**2 + 9 * x + 1 + (6 * x**2 + 6 * x) * np.log(x)
    chunk_2 = (x - 1)**5
    return 2.0 * x * (chunk_1 / chunk_2)
#coupling Y^2
def E_H(mass):#NLO
    y = NLOyy(mass)
    chunk_1 = 7 * y**3 - 36 * y**2 + 45 * y - 16
    chunk_2 = (18 * y - 12) * np.log(y)
    return  y / 36 * ( (chunk_1 + chunk_2) / ((y - 1)**4) )
print('E_H(mass)',E_H(charHm_100))
def w7_yy(mass):#NLO
    y = NLOyy(mass)
    
    chunk_1 = (8 * y**3 - 37 * y**2 + 18 * y) / ((y - 1)**4)
    chunk_2 = sp.spence( 1.0 / y)
    chunk_3 = (3 * y**3 + 23 * y**2 - 14 * y) * (np.log(y))**2 / ((y - 1)**5)
    chunk_4 = (21 * y**4 - 192 * y**3 - 174 * y**2 + 251 * y - 50) \
    * np.log(y) / (9 * (y - 1)**5)
    chunk_5 = (- 1202 * y**3 + 7569 * y**2 - 5436 * y + 797) / (108 * (y - 1)**4)
    return 2 * y / 9 * (chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5) - \
4 / 9 * E_H(mass)
def w8_yy(mass):#NLO
    y = NLOyy(mass)
    chunk_1 = (13 * y**3 - 17 * y**2 + 30 * y) / ((y - 1)**4)
    chunk_2 = sp.spence(1.0 / y)
    chunk_3 = - (17 * y**2 + 31 * y) / ((y - 1)**5) * (np.log(y))**2
    chunk_4 = (42 * y**4 + 318 * y**3 + 1353 * y**2 + 817 * y - 226) * \
np.log(y) / (36 * (y - 1)**5)
    chunk_5 = (- 4451 * y**3 + 7650 * y**2 - 18153 * y + 1130) / (216 * (y - 1)**4)
    return 1 / 6 * y * (chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5) - \
1 / 6 * E_H(mass)
def m7_yy(mass):#NLO
    y = NLOyy(mass)
    chunk_1 = - 14 * y**4 + 149 * y**3 - 153 * y**2 - 13 * y + 31
    chunk_2 = - (18 * y**3 + 138 * y**2 - 84 * y) * np.log(y)
    chunk_3 = (y - 1)**5
    return 1 / 27 * y * ( (chunk_1 + chunk_2) / chunk_3 )
def m8_yy(mass):#NLO
    y = NLOyy(mass)
    chunk_1 = - 7 * y**4 + 25 * y**3 - 279 * y**2 + 223 * y + 38 
    chunk_2 = (102 * y**2 + 186 * y) * np.log(y)
    chunk_3 = (y - 1)**5
    return 1 / 36 * y * ( (chunk_1 + chunk_2) / chunk_3)
def t7_yy(mass):#NLO
    y = NLOyy(mass)
    chunk_1 = 47 * y**3 - 63 * y**2 + 9 * y + 7 
    chunk_2 = - (18 * y**3 + 30 * y**2 - 24 * y) * np.log(y)
    return 1 / 3 * (y / 3 * ( (chunk_1 + chunk_2) / (y - 1)**5 ) )
def t8_yy(mass):#NLO
    y = NLOyy(mass)
    chunk_1 = - y**3 - 9 * y**2 + 9 * y + 1 + (6 * y**2 + 6 * y) * np.log(y)
    chunk_2 = (y - 1)**5
    return 1 / 3 * ( 2 * y * (chunk_1 / chunk_2) )

#Coupling (XY*)
def w7_xy(mass):#NLO
    y = NLOyy(mass)
    chunk_1 = (8 * y**2 - 28 * y + 12) / (3 * (y - 1)**3)
    chunk_2 = sp.spence(1.0 / y)
    chunk_3 = (3 * y**2 + 14 * y - 8) * (np.log(y))**2 / (3 * (y - 1)**4)
    chunk_4 = (4 * y**3 - 24 * y**2 + 2 * y + 6 ) * np.log(y) / (3 * (y - 1)**4)
    chunk_5 = (- 2 * y**2 + 13 * y  - 7) / ((y - 1)**3)
    return 4 / 3 * y * (chunk_1 * chunk_2 + chunk_3 + chunk_4 + chunk_5)
def w8_xy(mass):#NLO
    y = NLOyy(mass)
    c1 = (17 * y**2 - 25 * y + 36) / (2 * (y - 1)**3)
    c2 = sp.spence(1.0 / y)
    c3 = - (17 * y + 19) * (np.log(y))**2 / ((y - 1)**4)
    c4 = (14 * y**3 - 12 * y**2 + 187 * y + 3) * np.log(y) / (4 * (y - 1)**4)
    c5 = - (3 * (29 * y**2 - 44 * y + 143)) / (8 * (y - 1)**3)
    return 1 / 3 * y * (c1 * c2 + c3 + c4 + c5)
def m7_xy(mass):#NLO
    y = NLOyy(mass)
    c1 = - 8 * y**3 + 55 * y**2 - 68 * y + 21 
    c2 = - (6 * y**2 + 28 * y - 16) * np.log(y)
    c3 = (y - 1)**4
    return  2 / 9 * y * ( (c1 + c2) / c3 )
def m8_xy(mass):#NLO
    y = NLOyy(mass)
    c1 = - 7 * y**3 + 23 * y**2 - 97 * y  + 81
    c2 = (34 * y + 38) * np.log(y)
    c3 = (y - 1)**4
    return 1 / 6 * y * ( (c1 + c2) / c3 )
def t7_xy(mass):#NLO
    y = NLOyy(mass)
    c1 = 13 * y**2 - 20 * y + 7 - (6 * y**2 + 4 * y - 4) * np.log(y)
    c2 = (y - 1)**4
    return 2 / 3 * y * (c1 / c2)
def t8_xy(mass):#NLO
    y = NLOyy(mass)
    c1 = - y**2 - 4 * y + 5 + (4 * y + 2) * np.log(y)
    c2 = (y - 1)**4
    return 2 * y * (c1 / c2)
#####################################################
##################################################################
###################################################
#EQUATION 18
def c1_mu_effective(s,mass1,mass2,i1,j1,i2,j2):
    ratio1 = np.log(mt**2 / s**2)
    ratio = np.log(s**2 / mw**2)
    ratio2 = np.log(s**2/(mass1)**2)
    if np.any(mass2) == 0.0:
        ratio22 = 0.0
    else:
        ratio22 = np.log(s**2/(mass2)**2)
    ratio_muoverpi = NLOalpha_s(mw) / (4 * PI) 
#    print('yy()',yy())
#    print('w7_sm',w7_sm(),t7_sm() * (ratio1  - 4 / 3))
#    print('===',w8_yy(),w8_xy(),t8_xy(),t8_yy(),m8_xy())
    listsm7 =  w7_sm() + m7_sm() * ratio + \
        t7_sm() * (ratio1  - 4 / 3) #7
    listsm8 =  w8_sm() + m8_sm() * ratio + \
        t8_sm() * (ratio1  - 4 / 3) #8
    listyy7 =  w7_yy(mass1) + m7_yy(mass1) * ratio2 + \
        t7_yy(mass1) * (ratio1  - 4 / 3) #7
    listyy8 =  w8_yy(mass1) + m8_yy(mass1) * ratio2 + \
        t8_yy(mass1) * (ratio1  - 4 / 3) #8
    listxy7 =  w7_xy(mass1) + m7_xy(mass1) * ratio2 + \
        t7_xy(mass1) * (ratio1  - 4 / 3) #7
    listxy8 =  w8_xy(mass1) + m8_xy(mass1) * ratio2 + \
        t8_xy(mass1) * (ratio1  - 4 / 3) #8
    listyy72 =  w7_yy(mass2) + m7_yy(mass2) * ratio22 + \
        t7_yy(mass2) * (ratio1  - 4 / 3) #7
    listyy82 =  w8_yy(mass2) + m8_yy(mass2) * ratio22 + \
        t8_yy(mass2) * (ratio1  - 4 / 3) #8
    listxy72 =  w7_xy(mass2) + m7_xy(mass2) * ratio22 + \
        t7_xy(mass2) * (ratio1  - 4 / 3) #7
    listxy82 =  w8_xy(mass2) + m8_xy(mass2) * ratio22 + \
        t8_xy(mass2) * (ratio1  - 4 / 3) #8
#    print('LO_7sm()',c0_7sm())
#    print('LO_7xy()',c0_7xy())
#    print('LO_7yy()',c0_7yy())
#    print('NLO7_sm',listsm7 * ratio_muoverpi + c0_7sm())
#    print('NLO7_xy',listxy7 * ratio_muoverpi + c0_7xy())
#    print('NLO7_yy',listyy7 * ratio_muoverpi + c0_7yy())
#    print('w8_xy()', w8_xy(),m8_xy(),t8_xy(),ratio2)
#    print('listxy8',listxy8,listsm8,c0_8yy())
    
    result_NLO7 = listsm7 + np.abs(i1)**2 * listyy7 + np.abs(i2)**2 * listyy72 + \
                  np.array(j1) * listxy7 + np.array(j2) * listxy72
    result_LO7 = c0_7sm() + np.abs(i1)**2 * c0_7yy(mass1) + np.array(j1) * c0_7xy(mass1) +\
                  np.abs(i2)**2 * c0_7yy(mass2) + np.array(j2) * c0_7xy(mass2)
    result_LO2 = np.ones(len(result_LO7))
    result_NLO2 = np.zeros(len(result_NLO7))
    result_NLO8 = listsm8 + np.abs(i1)**2 * listyy8 + np.abs(i2)**2 * listyy82 + \
                  np.array(j1) * listxy8 + np.array(j2) * listxy82
    result_LO8 = c0_8sm() + np.abs(i1)**2 * c0_8yy(mass1) + np.array(j1) * c0_8xy(mass1) +\
                  np.abs(i2)**2 * c0_8yy(mass2) + np.array(j2) * c0_8xy(mass2)
#    print('result_LO8',result_LO8,result_NLO8,listsm8,listyy8,listxy8)
    result2 = result_LO2 + ratio_muoverpi * result_NLO2 
    result7 = result_LO7 + ratio_muoverpi * result_NLO7
    result8 = result_LO8 + ratio_muoverpi * result_NLO8
    return np.array([result_LO2,result_NLO2,result2,result_LO7,result_NLO7,result7,\
                     result_LO8,result_NLO8,result8])
#################################################################
# NLO Wislon coefficients at mu_b
#### Wilson coefficient at low scale(mu_b)
##############################LO
def C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2): #c0_7_eff(mu_b) LO
    eta = NLOalpha_s(s1) / NLOalpha_s(s2) # alpha_s (mu_w) / alpha_s(mu_b)
    step1 = eta**(16 / 23) * c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[3]
    step2 = (eta**(14 / 23) - eta**(16 / 23) ) * c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[6]
    result1 = 0.0
    for n in np.arange(0,8):
#        print('ai*hi*1.0',a_i[n] * h_i[n] * c0_eff(s1,i,j)[1])
        result1 += (eta)**(a_i[n]) * h_i[n] * c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[0]
    return step1 + 8 /3 * (step2) + result1
###############################NLO
def C1_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2): #c1_7_eff(mu_b) NLO
    eta = NLOalpha_s(s1) / NLOalpha_s(s2) # alpha_s (mu_w) / alpha_s(mu_b) NLO
    step1 = eta **(39 / 23) * c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[4]
#    print('c1_eff(s1,i,j)[6,m]',c1_eff(NLOalpha_s(mw),ll2,ll1)[6])
    step2 = 8 / 3 * (eta**(37 / 23) - eta**(39 / 23) ) *\
            c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[7]
    step3 = ((297664 / 14283) * eta**(16 / 23) - 7164416 / 357075 * eta**(14 / 23) +\
256868/14283 * (eta**(37 / 23)) - 6698884 / 357075 * eta**(39 / 23)) * \
             c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[6]
    step4 = 37208/4761 * (eta**(39 / 23) - eta**(16 / 23)) * \
            c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[3]
    result = 0.0
    c1_4effmw = E_0() + 2 / 3 * np.log(s1**2 / mw**2) + np.abs(i1)**2 * E_H(mass1) + \
                 np.abs(i2)**2 * E_H(mass2)
#    print('c1_eff(NLOalpha_s(mw),ll2,ll1)[3]',c1_eff(NLOalpha_s(mw),ll2,ll1)[3])
    for n in np.arange(0,8):
        step5 = e_i[n] * eta * c1_4effmw
        step6 = (f_i[n] + k_i[n] * eta) * c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[0]
        step7 = l_i[n] * eta * (15 + 6 * np.log(s1**2 / mw**2))
        result += (eta**(a_i[n])) * (step5 + step6 + step7)
    return step1 + step2 + step3 + step4 + result
def C0_1_eff(s2,s1,mass1,mass2,i1,j1,i2,j2):
    eta = NLOalpha_s(s1) / NLOalpha_s(s2) # alpha_s (mu_w) / alpha_s(mu_b)
    step1 =  eta**(6 / 23) - eta**(- 12 / 23) 
    return step1 * c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[0]
def C0_2_eff(s2,s1,mass1,mass2,i1,j1,i2,j2):
    eta = NLOalpha_s(s1) / NLOalpha_s(s2) # alpha_s (mu_w) / alpha_s(mu_b)
    step1 = 2 /3 * eta**(6 / 23) + 1 /3 * eta**(- 12 / 23) 
    return step1 * c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[0]
def C0_8_eff(s2,s1,mass1,mass2,i1,j1,i2,j2):
    eta = NLOalpha_s(s1) / NLOalpha_s(s2) # alpha_s (mu_w) / alpha_s(mu_b)
    step1 = eta**(14 / 23) * c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[6]
    result = 0.0
    for n in np.arange(0,5):
        result += h2_i[n] * eta**(a2_i[n]) * c1_mu_effective(s1,mass1,mass2,i1,j1,i2,j2)[0]
    return step1 + result
#####################################################################
####################################################################
####################################################################
def D_bar(s2,s1,mass1,mass2,i1,j1,i2,j2):#mu_b scale Reduced Amplitude LO
    # Riemann Zeta func- tion zeta_3
    zeta_3 = 1.2021
    L = np.log(zz)
    r2 = complex(2/243 * (-833 + 144 * PI**2 * zz**(3/2) + (1728 - 180 * PI**2 - 1296 * zeta_3 + \
       (1296 - 324 * PI**2) * L + 108 * L**2 + 36 * L**3 ) * zz + \
       (648 + 72 * PI**2 + (432 - 216 * PI**2) * L + 36 * L**3 ) * zz**2 + \
       (- 54 - 84 * PI**2 + 1092 * L - 756 * L**2) * zz**3) ,\
        16 * PI / 81 * (- 5 + (45 - 3 * PI**2 + 9 * L + 9 * L**2) * zz + \
        (- 3 * PI**2 + 9 * L**2) * zz**2 + (28 - 12 * L) * zz**3))
    r1 = - 1 / 6 * r2
    r7 = 32 / 9 - 8 / 9 * PI**2
    r8 = - 4 / 27 * complex(- 33 + 2 * PI**2, - 6 * PI )
    ans1 = C0_1_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) * \
    (r1 + 1 / 2.0 * gaeff.gamma0eff()[0][6] * np.log(mb**2 / s2**2 ))
    ans2 = C0_2_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) * \
    (r2 + 1 / 2.0 * gaeff.gamma0eff()[1][6] * np.log(mb**2 / s2**2 ))
    ans7 = C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) * \
    (r7 + 1 / 2.0 * gaeff.gamma0eff()[6][6] * np.log(mb**2 / s2**2 ))
    ans8 = C0_8_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) * \
    (r8 + 1 / 2.0 * gaeff.gamma0eff()[7][6] * np.log(mb**2 / s2**2 ))
    v_ub = ans1 + ans2 + ans7 + ans8 - 16 / 3.0 * \
    C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2)
    return C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) + \
NLOalpha_s(s2) / (4 * PI) *  ((C1_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) +  v_ub) )
####################################################################
def NLOD_bar(s2,s1,mass1,mass2,i1,j1,i2,j2):#mu_b scale Reduced Amplitude NLO
    zeta_3 = 1.2021
    L = np.log(zz)
    r2 = complex(2/243 * (-833 + 144 * PI**2 * zz**(3/2) + (1728 - 180 * PI**2 - 1296 * zeta_3 + \
       (1296 - 324 * PI**2) * L + 108 * L**2 + 36 * L**3 ) * zz + \
       (648 + 72 * PI**2 + (432 - 216 * PI**2) * L + 36 * L**3 ) * zz**2 + \
       (- 54 - 84 * PI**2 + 1092 * L - 756 * L**2) * zz**3) ,\
        16 * PI / 81 * (- 5 + (45 - 3 * PI**2 + 9 * L + 9 * L**2) * zz + \
        (- 3 * PI**2 + 9 * L**2) * zz**2 + (28 - 12 * L) * zz**3))
    r1 = - 1 / 6 * r2
    r7 = 32 / 9 - 8 / 9 * PI**2
    r8 = - 4 / 27 * complex(- 33 + 2 * PI**2, - 6 * PI )
    ans1 = C0_1_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) * \
    (r1 + 1 / 2.0 * gaeff.gamma0eff()[0][6] * np.log(mb**2 / s2**2 ))
    ans2 = C0_2_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) * \
    (r2 + 1 / 2.0 * gaeff.gamma0eff()[1][6] * np.log(mb**2 / s2**2 ))
    ans7 = C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) * \
    (r7 + 1 / 2.0 * gaeff.gamma0eff()[6][6] * np.log(mb**2 / s2**2 ))
    ans8 = C0_8_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) * \
    (r8 + 1 / 2.0 * gaeff.gamma0eff()[7][6] * np.log(mb**2 / s2**2 ))
    v_ub = ans1 + ans2 + ans7 + ans8 - 16 / 3.0 * \
    C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2)
    chunk1 = (C1_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) + v_ub )
    return NLOalpha_s(s2) * chunk1 / (4 * PI * C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2))
####################################################################
###########Decay_width of b > s gamma 
def decay_bsp(s2,s1,mass1,mass2,i1,j1,i2,j2):
    return gf**2 / (32 * PI**4) * (vts * vtb)**2 * \
    alpha_electroweak * mb**5 * np.abs(D_bar(s2,s1,mass1,mass2,i1,j1,i2,j2))**2
####################################################################
def Amp(s2,s1,mass1,mass2,i1,j1,i2,j2):# A for Decay_width of b > s gamma gluon
        c0_1 = C0_1_eff(s2,s1,mass1,mass2,i1,j1,i2,j2)
        c0_2 = C0_2_eff(s2,s1,mass1,mass2,i1,j1,i2,j2)
        c0_7 = C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2)
        c0_8 = C0_8_eff(s2,s1,mass1,mass2,i1,j1,i2,j2)
        f_22 = 16 * zz / 27 * quad(grand1, 0, 1.0 / zz  )[0] 
        print('f_22',f_22)
        f_27 = - 8 * zz**2 / 9 * complex(quad(grand2, 0, 1.0 / zz )[0], quad(grand3, 0, 1.0 / zz )[0] ) 
        print('f_27',f_27)
        f_11 = 1 /36 * f_22
        f_12 = - 1 /3 * f_22
        f_17 = - 1 /6 * f_27
        f_28 = - 1 /3 * f_27
        f_18 = - 1 /6 * f_28
        f_78 = 8 /9 * (25 /12 - PI**2 / 6)
        f_88 = 1 /27 * (16 /3 - 4 * PI**2 / 3 + 4 * np.log(mb / s2 ))
        summ11 = c0_1 * np.conjugate(c0_1) * f_11
        summ12 = c0_1 * np.conjugate(c0_2) * f_12
        summ17 = c0_1 * np.conjugate(c0_7) * f_17
        summ18 = c0_1 * np.conjugate(c0_8) * f_18
        summ22 = c0_2 * np.conjugate(c0_2) * f_22
        summ27 = c0_2 * np.conjugate(c0_7) * f_27
        summ28 = c0_2 * np.conjugate(c0_8) * f_28
        summ78 = c0_7 * np.conjugate(c0_8)* f_78
        summ88 = c0_8 * np.conjugate(c0_8) * f_88
        summ_all = summ11 + summ12 + summ17 + summ22 + summ27 + summ28\
    + summ78 + summ88 + summ18
        return NLOalpha_s(s2) / PI * summ_all.real
#print('Amp',Amp(ll2,ll1))
###########Decay_width of b > s gamma gluon
def decay_bspg(s2,s1,mass1,mass2,i1,j1,i2,j2):
        a1 = gf**2 / (32 * PI**4) * (vts * vtb)**2 * \
        alpha_electroweak * mb**5 
        return a1 * Amp(s2,s1,mass1,mass2,i1,j1,i2,j2)  
###########Decay_width of semileptonic 
def decay_SL():
    part1 = gf**2 /(192 * PI**3) * np.abs(vcb)**2 * mb**5 * g_z
    part2 = 1 - 2 * NLOalpha_s(run_quark_bar(mb)) * f_z / (3 * PI) \
    + delta_NP_SL / mb**2
    return part1 * part2
#print('Partial width of semileptonic decay', decay_SL() )
#################################################################
#Measured Semi- leptonic branching ratio B_SL
B_SL = 0.1049 # Phys. Rev. Lett. 76, 1570 – Published 4 March 1996 =  0.1049
#################################################################
################################################################
#################### Partial width of B_bar > X_s + gamma
def decay_B_bar_Xsg(s2,s1,mass1,mass2,i1,j1,i2,j2):
    a1 = gf**2 / (32 * PI**4) * np.abs(vts * vtb)**2 * alpha_electroweak * mb**5 
    chunk1 = np.abs(D_bar(s2,s1,mass1,mass2,i1,j1,i2,j2))**2 + \
             Amp(s2,s1,mass1,mass2,i1,j1,i2,j2) + delta_NP_ga / (mb**2) * \
             np.abs(C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2))**2 + \
             (delta_NP_c / mc**2) * \
    (C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2).conjugate() *\
     (C0_2_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) - \
      C0_1_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) * 1 / 6)).real
    return a1 * chunk1
def BR_B_Xs_gamma(s2,s1,mass1,mass2,i1,j1,i2,j2):
    return decay_B_bar_Xsg(s2,s1,mass1,mass2,i1,j1,i2,j2) \
          / decay_SL() * B_SL 
############
#  print('ll2,ll1',ll2,ll1,len(ll2))
#  print('xyfun-list',xyfun_list[m],yfun_list[m],len(xyfun_list[m]))
#  print('c0_7eff(i,j)',c0_7eff(yfun_list[m],xyfun_list[m]),len(c0_7eff(yfun_list[m],xyfun_list[m])))
#  print('c0_8eff(i,j)',c0_8eff(yfun_list[m],xyfun_list[m]),len(c0_8eff(yfun_list[m],xyfun_list[m])))
#  print('c0_eff(LOalpha_s(mw),i,j)',c0_eff(LOalpha_s(mw),yfun_list[m],xyfun_list[m]))
#  print('c1_eff(mu_w,i,j)',c1_eff(NLOalpha_s(mw),ll2,ll1), type( c1_eff(NLOalpha_s(mw),ll2,ll1)    ))
#  print('c1_eff(s1,i,j)[3]',c1_eff(NLOalpha_s(mw),ll2,ll1)[3])
#print('c_i_eff_muw(mu_w,i,j)',c_i_eff_muw(ll2,ll1 ))
#################################################################
##################################################################
#  print('C0_7_eff(mu_b)',c0_7_eff(mb,mw,yfun_list[m],xyfun_list[m])) #
#  print('C1_7_eff(mu_b)',c1_7_eff(mb,mw,yfun_list[m],xyfun_list[m]))
#  print('C0_2_eff(mu_b)',c0_2_eff(mb,mw,yfun_list[m],xyfun_list[m])) #
#  print('C0_8_eff(mu_b)',c0_8_eff(LOalpha_s(mb),LOalpha_s(mw),ll2,ll1)) #
#print('D_bar',D_bar(0.8,0.2), type(D_bar(0.8,0.2)))
#print('delta_D_bar',delta_D_bar(0.8,complex(0.2, 0.1)))
print('---------------------------------------------------------')
##################################################################
##################################################################
######### A_CP CP asymmetry  expression 
######### DOI: 10.1103/PhysRevD.58.094012 
def v_cp(i):
    return (5 + np.log(i) + np.log(i)**2 - PI**2 / 3) + \
 (np.log(i)**2 - PI**2 /3 ) * i + (28/ 9 - 4 /3 * np.log(i)) * i**2
def g(i,y): # i will = zz, y will = delta_cp
    stepfunction = np.heaviside(4 * i / y,1.0)
    part1 = (y**2 - 4 * y * i + 6 * i**2) * np.log(np.sqrt(y / (4 * i)) + \
            np.sqrt(y /(4 * i) + 1 ) ) 
    part2 = 3 * y * (y - 2 * i) / 4 * np.sqrt(stepfunction) 
    return  part1 - part2 
def b_cp(i,x):#delta_cp fraction of Energy cut
    return g(i,1) - g(i, 1 - x)
def A_cp(s2,s1,mass1,mass2,i1,j1,i2,j2): # CP asymmetry 
    epsilon_s = lanmda_ckm**2 * complex(- rho_ckm ,eta_ckm)#e_s = V*_usV_ub / V*_tsV_tb
    c2 = C0_2_eff(s2,s1,mass1,mass2,i1,j1,i2,j2)
    c7 = C0_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2) + NLOalpha_s(s2) / \
         (4 * PI) * C1_7_eff(s2,s1,mass1,mass2,i1,j1,i2,j2)
    c8 = C0_8_eff(s2,s1,mass1,mass2,i1,j1,i2,j2)
    part1 = NLOalpha_s(s2) / (np.abs(c7) )**2
    part2 = 40 / 81 * (c2 * np.conjugate(c7)).imag 
    part3 = 8 * zz / 9 * (v_cp(zz) + b_cp(zz, delta_cp)) *\
        ( (1 + epsilon_s) * (c2 * np.conjugate(c7) ) ).imag
    part4 = 4 / 9 * (c8 * np.conjugate(c7) ).imag
    part5 = 8 * zz / 27 * b_cp(zz, delta_cp) *\
        ( (1 + epsilon_s) * (c2 * np.conjugate(c8) ) ).imag
    return part1 * (part2 - part3  - part4 + part5)
print(yy(charHm_100))
print('---------------------------------------------------------')
print('BR(X_bar>Xs+gamma)',\
          BR_B_Xs_gamma(mb,100,100,100 + 20,\
                        [1.0],1j * np.array(np.arange(-5,6,1.0)),[0],[0]) )
print('A_CP',np.sort(A_cp(mb,mw,charHm_100,charHm_100 + 20 ,[1.0],[1.0],[0],[0])))
for m in np.arange(0,len(charHm_100)):
    print(m,'---------------------')
    print('charHm_100[m]',charHm_100[m],charHm_100[m]+ 20,m)
    print('Ti_sm,xy,yy()',t7_sm(),'xy',t7_xy(charHm_100[m]),'yy',t7_yy(charHm_100[m]))#
    print('wi_sm,xy,yy()',w7_sm(),'xy',w7_xy(charHm_100[m]),'yy',w7_yy(charHm_100[m]))#
    print('c1_mu_effective(scale,mhch),sm,xy,yy',\
          c1_mu_effective(mw,charHm_100[m],charHm_large[m],[1.0],[-1.0],[0],[0] ))
    print('C0_7_eff(s2,s1,i,j)',\
          C0_7_eff(mb,mw,charHm_100[m],charHm_large[m],[1.0],[-1.0],[0],[0]))
    print('C0_8_eff(s2,s1,i,j)',\
          C0_8_eff(mb,mw,charHm_100[m],charHm_large[m],[1.0],[-1.0],[0],[0]))
def Plot_3():#Figure 3 DOI: 10.1142/S0217751X17501457
    x_axis = np.array([ i for i in np.arange(100,1020,20)] )
    print(x_axis,type(x_axis),x_axis.shape)
    y11_axis = np.array(BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,[1.0],[-1.0],[0],[0]) )
    y12_axis = np.array(BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,[1.0/2.0],[-1.0/4.0],[0],[0]) )
    y130_axis = np.array(BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,[1.0/30.0],[-1.0/900],[0],[0]) )
    y21_axis = np.array(BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,[1.0],[1.0],[0],[0]) )
    y22_axis = np.array(BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,[1.0/2.0],[1.0],[0],[0]) )
    y230_axis = np.array(BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,[1.0/30.0],[1.0],[0],[0]) )
    plt.axis([100.0, 1000.0, 1.0, 7.0])
    plt.plot(x_axis,y11_axis / (1e-4))
    plt.plot(x_axis,y12_axis / (1e-4))
    plt.plot(x_axis,y130_axis / (1e-4))
    plt.plot(x_axis,y21_axis / (1e-4))
    plt.plot(x_axis,y22_axis / (1e-4) )
    plt.plot(x_axis,y230_axis / (1e-4))
    print('type I tanbeta = 1',y11_axis / (1e-4))
    plt.xlabel('$M_{H^{\pm}}$')
    plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.legend(('Type I tan$\\beta =$ 1', 'Type I tan$\\beta =$ 2', 'Type I tan$\\beta =$ 30',\
            'Type II tan$\\beta =$ 1', 'Type II tan$\\beta =$ 2', 'Type II tan$\\beta =$ 30'),
           loc='upper right', shadow=True,prop={'size': 7.8})
    plt.show()
    plt.close
###################################
def Plot_8_9():
    xim_axis = np.arange(-5.0,5.2,0.1)# figure 8
    XYimx_axis = [np.complex(-2,i)*1.0 for i in xim_axis]
    rangephi = np.arange(0,190,10) # figure 9
    print('rangephi', rangephi,len(rangephi))
    XYexpim_axis = [2 * complex(np.cos(i),np.sin(i)) * 0.5 for i in rangephi] 
    print('REALX,IMX:',[np.complex(-2,i)*1.0 for i in xim_axis])
    print('X = 2exp(i phi)',XYexpim_axis,len(XYexpim_axis))
   #mlx = MultipleLocator(1)
   #mly = MultipleLocator(0.25)
    y48imx_axis = BR_B_Xs_gamma(4.8,mhch,mhch,100 ,\
                        [1.0],XYimx_axis,[0.0],[0.0])
    y24imx_axis = BR_B_Xs_gamma(2.4,mhch,mhch,100,\
                        [1.0],XYimx_axis,[0.0],[0.0])
    y96imx_axis = BR_B_Xs_gamma(9.6,mhch,mhch,100,\
                        [1.0],XYimx_axis,[0.0],[0.0])
    y48phi_axis = BR_B_Xs_gamma(4.8,mhch,mhch,100 ,\
                        [0.5],XYexpim_axis,[0.0],[0.0])
    y24phi_axis = BR_B_Xs_gamma(2.4,mhch,mhch,100,\
                        [0.5],XYexpim_axis,[0.0],[0.0])
    y96phi_axis = BR_B_Xs_gamma(9.6,mhch,mhch,100,\
                        [0.5],XYexpim_axis,[0.0],[0.0])

    plt.xlim(-7, 7)
    plt.ylim(-2, 6.5)
#plt.axes().xaxis.set_minor_locator(mlx)
#plt.axes().yaxis.set_minor_locator(mly)
    plt.plot(xim_axis,y48imx_axis / (1e-4))
    plt.plot(xim_axis,y24imx_axis / (1e-4))
    plt.plot(xim_axis,y96imx_axis / (1e-4))
#plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.xlabel('Im(X)')
    plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.legend(('$\mu = 4.8$ GeV', '$\mu = 2.4$ GeV', '$\mu = 9.6$ GeV '),
           loc='upper right', shadow=True,prop={'size': 8})
    plt.show()
    plt.close
    plt.plot(rangephi,y48phi_axis / (1e-4))
    plt.plot(rangephi,y24phi_axis / (1e-4))
    plt.plot(rangephi,y96phi_axis / (1e-4))
    plt.xlabel('$\\phi$')
    plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.legend(('$\mu = 4.8$ GeV', '$\mu = 2.4$ GeV', '$\mu = 9.6$ GeV '),
           loc='upper right', shadow=True,prop={'size': 8})
    plt.show()
    plt.close
########################################################
Plot_3()
Plot_8_9()
#M = np.ones((3, 2))
#a = np.arange(3)
#print(M,a)
#print(a[:,np.newaxis])
#print(M * a[:,np.newaxis])

