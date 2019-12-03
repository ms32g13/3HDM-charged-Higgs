#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 13:20:18 2019
EDM of neutron 
@author: muyuansong0419
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import special as sp
from invariant import *
from scipy.integrate import quad
################################
#Charged Higgs contribution to NEUTRON EDM from paper :
#JHEP04(2014)076 
################################
mhadron = 1 # 1 GeV hadronic scale
#########################  
def beta_0(i): # beta_0 function below EQ2.33
    nc = 3
    nf = i
    return (11.0 * nc - 2.0 * nf)/3.0
def LOalpha_s(s,i): #alpha_s(mu) at LO
#    print(i,zz,alpmz,mz)
    beta0 = beta_0(i)
    v = 1 - beta0 * (alpmz / (2.0 * PI)) * np.log(mz / s)
    return alpmz / v 
################ NLO for strong coupling constant
def NLOalpha_s(s,i): #alpha_s(mu) at NLO
    beta0 = beta_0(i)
    beta1 = 116 / 3
    v = 1 - beta0 * (alpmz / (2.0 * PI)) * np.log(mz / s)
    ratio1 = np.log(v) / v
    return alpmz / v * (1 - (beta1 / beta0) * (alpmz / (4 * PI)) * ratio1)
def run_quark_bar(s,m,i):#  NLO
    c1 = np.log(m**2 / s**2)
    c2 = NLOalpha_s(s,i) / PI
    return m * (1 + c2 * c1 - 4 / 3 * c2 )
#########################
# Strong coupling constant of specific scale Q   One-loop
# from paper: PHYSICAL REVIEW D 81, 075016 (2010) Appendix A(4)
def gs(Q,i): 
    b0 = beta_0(i)
    denominator = 1.0 + b0 * NLOalpha_s(mz,i) / (2 * PI) *\
                  np.log(Q / mz)
    return NLOalpha_s(mz,i) / denominator
print('Strong coupling constant at scale '+ str(173) +':', gs(173,5))
print('running mass of b quark at 100 GeV',run_quark_bar(100,4.18,5))
print('mu_w scale at LO and NLO:',mw,LOalpha_s(mw,5),NLOalpha_s(mw,5))
print('mu_z scale at LO and NLO:',mz,LOalpha_s(mz,5),NLOalpha_s(mz,5))
print('mu_b scale at LO and NLO:',mb,LOalpha_s(mb,5),NLOalpha_s(mb,5),run_quark_bar(4.18,4.89,5))
print('mu_t scale at LO and NLO:',mt,LOalpha_s(mt,6),NLOalpha_s(mt,6))
print('168 scale at LO and NLO:',168,LOalpha_s(168,5),NLOalpha_s(168,5))
print('mu_c scale at LO and NLO:',mc,LOalpha_s(mc,4),NLOalpha_s(mc,4))
##########################    
def kc_p(e):# # EQ 2.33 gamma_{C_photon} / (2 * beta0) e: electric charge
    return  8 * e * 4 /3 
def kw(i):# # EQ 2.33 gamma_w / (2 * beta0) i: number of active color
    nf = i
    nc = 3
    return (nc + 2 * nf)/ (2.0 * beta_0(i) )
def kc(i):# # EQ 2.33 gamma_c / (2 * beta0)  i: number of active color
    nc = 3
    return (10 * 4 / 3 - 4 * nc)/ (2.0 * beta_0(i) )
def kgamma(i):# # EQ 2.33 gamma_gamma / (2 * beta0)
    return 2 * 4/3 / (2.0 * beta_0(i) )
def x_tH(mass):
    return  mt**2 / (mass)**2
def eta(x,y,i,j,k1,k2):# at NLO for EQ 2.35,2.36,2.37 
    # i number of active quark at high scale; j at  low scale 
    # x is high scale, y is low scale
    # k1: for higher scale ; k2: for lower scale
    return LOalpha_s(x,i)**k1 / LOalpha_s(y,j)**k2
#################################
def h_m_M(m,M): # EQ 4.14 m = quark mass, M = neutral Higgs (Additional Higgs) mass:
    def integrand(m,M,x,u):# EQ 4.14  intergrand intergrate respect with u
        numerator = u**3 * x**3 * (1 - x)
        part1 = m**2 * x * (1 - u * x)
        part2 = M**2 * (1 - u) * (1 - x)
        return numerator / (part1 + part2)**2
    def integ1(m,M,x): # EQ 4.14 intergrate respect with x
        return quad(integrand, 0.0, 1.0, args=(m,M,x))[0]
    return quad(lambda x: integ1(m,M,x),0.0,1.0)[0] * m**4 / 4.0
def h_m_Mlarge(m,M): #EQ 4.15 M >> m function
    part = M**2 / m**2
    return 1.0 / (4.0 * part ) * (np.log(part) - 3.0 / 2.0 )
#print('m = quark, M = neutral Higgs',h_m_M(mt,100) )
#print('M>>m',h_m_Mlarge(mc,150) )
#print(x_tH(np.arange(100,500,100)))
def toverH(x):# 2nd contibution in EQ 4.11
    return  (np.log(x) / (x - 1)**3 + \
                (x - 3) / (2.0 * (x - 1)**2)  )
def dC_btH(mass1,mass2,j1,j2): # d^{C}_{b}(mu_tH) EQ 4.11
    part1 = - gf  / ( np.sqrt(2) * 16 * PI**2) * np.abs(vtb)**2 * run_quark_bar(mt,mb,5)
#    print(part1,toverH(x_tH(mass1)), toverH(x_tH(mass2)) )
    part2 = np.array(j1).imag * x_tH(mass1) * toverH(x_tH(mass1)) \
          + np.array(j2).imag * x_tH(mass2) * toverH(x_tH(mass2))
    return part1 * part2
#print('dC_btH', dC_btH(mb,np.arange(100,1000,10),100,[complex(1.0,0.1)],[0.0]) )
def C_wmutH():       # EQ 4.13 neutral  Higgs contribution
    return 0.0       # keep no neutral contribution

#######################################
def zx(x):#EQ 4.23
        return x**2 / mt**2
def tt(x):#EQ 4.23
        return (1 - 3 * x ) / (x**2) * PI**2 /6 + (1/x - 5/2) * np.log(x) -\
        1/x - (2 - 1/x) * (1 - 1/x) * sp.spence(x ) 
def tb(x):#EQ 4.23
        return (2 * x - 1) / (x**2) * PI**2 / 6 + (3/2 - 1/x) * np.log(x) +\
        1/x - (1/x) * (2 - 1/x) *  sp.spence(x)  
def Ft_(mass):#EQ4.23 for q = t
        tth = tt(zx(mass))
        ttw = tt(zx(mw))
        return (tth - ttw)/(zx(mass) - zx(mw))
def Fb_(mass):#EQ4.23 for q = b
        tbh = tb(zx(mass))
        tbw = tb(zx(mw))
        return (tbh - tbw)/(zx(mass) - zx(mw))
def dgmma_dBZ(mass,j):# EQ 4.24 charged Higgs contribution
    part1 = - 12 * gf**2 * mw**2 / (4 * PI)**4 * mb
    part2 = np.abs(vtb)**2 * np.abs(vud)**2 
    part3 = np.array(j).imag * (Ft_(mass) * 2/3 - 1/3 * Fb_(mass) )
    return part1 * part2 * part3 
#print('C_wmuh(,mass1,mass2,j1,j2)',\
#      C_wmuh(mb,np.arange(100,500,50),600,[complex(0.0,0.1)],[complex(0.0,0.1)]))

#(C) EDM
def dn_CEDM(mass1,mass2,j1,j2): # 2.3
    f_pi = 0.13 # PION decay constant 
    # using PI-0 meson , so use mass of PI-0 
    m_pi = 0.135
    ga_c = 10 * 4 /3 - 4 * 3 # gamma_c
    ga_cga = 8 * 4 /3 * (-1/3) # gamma_{C_gamma} for bottom quark (-1/3)
    ga_ga = 2 * 4/3 # gamma_{gamma}
    part2 = - f_pi**2 * m_pi**2 / (mb + 0)
    c1 = eta(mt,mhadron,6,3,kgamma(6),kgamma(3)) *  (dgmma_dBZ(mass1,j1) + dgmma_dBZ(mass2,j2) )
    dgammad_muh  = c1 + ga_cga / (ga_ga - ga_c) * \
            ( eta(mt,mhadron,6,3,kgamma(6),kgamma(3)) - eta(mt,mhadron,6,3,kc(6),kc(3))) *  \
            dC_btH(mass1,mass2,j1,j2)
    
    dcq_muh = eta(mc,mhadron,4,3,kc(4),kc(3)) * eta(mb,mc,5,4,kc(5),kc(4)) * \
            eta(mt,mb,6,5,kc(6),kc(5)) * dC_btH(mass1,mass2,j1,j2)
    
    return (1.4 * dgammad_muh  + 1.1 * dcq_muh ) * 1.0 * (part2 / (0.225)**3)
    
#print(dn_CEDM(80,400,complex(100,-1),complex(100,0.1))) 
# Weinberg Operator 
def CW(mass1,mass2,j1,j2):   # EQ 2.37 and 2.5
    part1 = eta(mc,mhadron,4,2,kw(4),kw(2)) * eta(mb,mc,5,4,kw(5),kw(4))
    part2 = eta(mt,mb,6,5,kw(6),kw(5)) * C_wmutH()
    part3 = eta(mt,mb,6,5,kc(6),kc(5)) * gs(mb,5)**3 / (8 * PI**2 * mb) \
          * dC_btH(mass1,mass2,j1,j2)      
    return part1 * (part2 + part3) * 1 * 0.02  # 20 MeV 
# Total Neutron EDM contribution from charged Higgs in 3HDM
def dn(mass1,mass2,j1,j2):
    return  CW(mass1,mass2,j1,j2) #+ dn_CEDM(mass1,mass2,j1,j2)
print(CW(100,500,complex(1,0.35),0) )