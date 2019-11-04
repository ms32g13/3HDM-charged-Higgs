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
PI = np.pi
#########################  
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
def beta_0(nc,nf): # beta_0 function below EQ2.33
    return (11.0 * nc - 2.0 * nf)/3.0
#########################
# Strong coupling constant of specific scale Q   One-loop
# from paper: PHYSICAL REVIEW D 81, 075016 (2010) Appendix A(4)
def gs(Q): 
    nf = 5
    b0 = beta_0(3.0,nf)
    denominator = 1.0 + b0 * LOalpha_s(mz) / (2 * PI) *\
                  np.log(Q / mz)
    return LOalpha_s(mz) / denominator
#print('Strong coupling constant at scale '+ str(173) +':', gs(173))
##########################    
def eta(x,y):
    return LOalpha_s(x) / LOalpha_s(y)
def kw(i):# i = nf
    return (9 + 2 * i)/ (2.0 * beta_0(3.0,i) )
def kc(i):# i = nf
    return (10 * 4 / 3 - 4 * i)/ (2.0 * beta_0(3.0,i) )
def x_tH(mass):
    return  mt**2 / (mass)**2
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
def dC_btH(s2,mass1,mass2,j1,j2): # d^{C}_{b}(mu_tH) EQ 4.11
    def toverH(x):
        return x * (np.log(x) / (x - 1)**3 + \
                (x - 3) / (2.0 * (x - 1)**2)  )
    part1 = - gf * 2.0 / ( np.sqrt(2) * 16 * PI**2) * np.abs(vtb)**2 * mb
#    print(part1,toverH(x_tH(mass1)), toverH(x_tH(mass2)) )
    part2 = np.array(j1).imag * toverH(x_tH(mass1)) \
            + np.array(j2).imag * toverH(x_tH(mass2))
    return part1 * part2
#print('dC_btH', dC_btH(mb,np.arange(100,1000,10),100,[complex(1.0,0.1)],[0.0]) )
def C_wmutH():       # EQ 4.13 neutral  Higgs contribution
    return 0.0       # keep no neutral contribution
def C_wmuh(s2,mass1,mass2,j1,j2):   # EQ 2.37 
    part1 = eta(mc,mh)**kw(3) * eta(mb,mc)**kw(3) 
    part2 = eta(mt,mh)**kw(3) * C_wmutH()
    part3 = eta(mt,mb)**kc(3) * gs(mb)**3 / (8 * PI**2 * mb) \
          * dC_btH(s2,mass1,mass2,j1,j2) /2.0              
    return part1 * (part2 + part3)
#######################################
def zx(x):
        return x**2 / mt**2
def tt(x):
        return (1 - 3 * x ) / (x**2) * PI**2 /6 + (1/x - 5/2) * np.log(x) -\
        1/x - (2 - 1/x) * (1 - 1/x) * sp.spence(1.0 - x) 
def tb(x):
        return (2 * x - 1) / (x**2) * PI**2 / 6 + (3/2 - 1/x) * np.log(x) +\
        1/x - (1/x) * (2 - 1/x) *  sp.spence(1.0 - x) 
def dgmma_dBZ(mass1,mass2,j1,j2):# EQ 4.24 charged Higgs contribution
    part1 = 12 * gf**2 * mw**2 / (4 * PI)**4
    part2 = np.abs(vtb)**2 * np.abs(vud)**2 
    part3 = np.array(j1).imag * (2 /3 * (tt(zx(mass1)) - tt(zx(mw)))/ (zx(mass1) - zx(mw)) - \
              1/3 * 2 /3 * (tb(zx(mass1)) - tb(zx(mw)))/ (zx(mass1) - zx(mw)))
    part4 = np.array(j2).imag * (2 /3 * (tt(zx(mass2)) - tt(zx(mw)))/ (zx(mass2) - zx(mw)) - \
              1/3 * 2 /3 * (tb(zx(mass2)) - tb(zx(mw)))/ (zx(mass2) - zx(mw)))
    return part1 * part2 * (part3 + part4)
#print('C_wmuh(s2,mass1,mass2,j1,j2)',\
#      C_wmuh(mb,np.arange(100,500,50),600,[complex(0.0,0.1)],[complex(0.0,0.1)]))


