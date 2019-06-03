#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:13:53 2019

@author: muyuansong0419

3HDM Functions
"""
from invariant import *
import math
import numpy as np
import CEPCeehh as cepc
import eehh as fl
import matplotlib.pyplot as plt
#4*pi*SQRT(2); factor appears in partial width of 2 fermions
fac = 4.0 * np.sqrt(2.0) * PI 
#############################
#the = random.uniform(- math.pi/2, math.pi/2)#theta (i in loop)
#tbe = random.uniform(1.0,61.0) #tanbeta (j in loop)
#tga = random.uniform(1.0,61.0) #tangamma (k in loop)
#delta = random.uniform(0.0,2 * math.pi)#delta (l in loop)  #delta fixed#
the = np.arange(- PI / 2, PI/50 ,PI / 10)#theta (i in loop)
tbe = np.arange(1.0,61.0,1.0) #tanbeta (j in loop)
tga = np.arange(1.0,61.0,1.0) #tangamma (k in loop)
delta = np.arange(0.0,2.1 * PI , PI /6)#delta (l in loop)  #delta fixed#
A = []
B = []
read1 = str('')
read2 = str('')
i = - PI / 2.1 #theta
j = 40.0    #tangentbeta
k = 10.0 #tangamma
l = 0.0 # set delta (phase shift) to 0
x = np.arange(0.0,40.2,0.2) # x range
y = np.arange(0.0,0.62,0.02) # y range
z = np.arange(0.0,5.02,0.02) # z range
# set contour line value
linecs = np.arange(0.0,1.2,0.2) # percentage limit 
linecs1 = (0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01)   
linecs2 = (0.001,0.005,0.01,0.02,0.05,0.10) 
linecs3 = np.arange(-1.1,1.1,0.4) # REALXY Limit
linecs4 = np.arange(-0.1,0.11,0.01)#IMXY Limit
linecs5 = (0.0001,0.0005,0.001 ,0.002 ,0.005 ,0.01 )
print('--------------------------------------------------')   
# charged H - tau, neutrino(tau-type)
def lamtn (Z):
      couptn = Z**2
      return gf * mhch * (mtau**2) * couptn / fac
# Charged H - charm, strange quark
def lamcs (X,Y):   
      coupc = Y**2
      coups = X**2 
      return (3 * gf * mhch * (vcs**2)) * (mc**2 * coupc + ms**2 * coups) * \
           (1.0 + 17.0 * alpmz/(3.0 * PI)) / fac

# Charged H - charm, bottom quark
def lamcb (X,Y):
      coupc = Y**2
      coupb = X**2
      return (3 * gf * mhch * (vcb**2)) * (mc**2 * coupc + mb**2 * coupb) * \
           (1.0 + 17.0 * alpmz/(3.0 * PI)) / fac
####################################################
# BR : charged Higgs - tau , neutrino(tau)
def brtn (X,Y,Z):
       return lamtn(Z)/(lamtn(Z)+lamcs(X,Y)+lamcb(X,Y)+lamaw)
#       return 0.0
# BR : charged Higgs - charm , strange quark
def brcs (X,Y,Z):
       return lamcs(X,Y)/(lamtn(Z)+lamcs(X,Y)+lamcb(X,Y)+lamaw)
# BR : charged Higgs - charm, bottom quark
def brcb (X,Y,Z):
#       return lamcb(X,Y)/(lamtn(Z)+lamcs(X,Y)+lamcb(X,Y)+lamaw)
    return 0.05
# BR : charged Higgs - A ,W
#braw=LamAW/(LamTN+LamCS+LamCB+LamAW) 
############################################################################
#############################################################
###  Here evaluate the Partial  widths and Branching Ratio of t-> H+,b and t-> w,b
# Partial width for t-> H+,b
def lamtHb(X,Y):
       return (0.5/fac) * gf * mt * (mt**2 * Y**2 + mb**2 * X**2) * \
           (1.0 - mhch**2/mt**2)**2
# Partial width for t-> w,b
lamtWb = (0.5/fac) * gf * mt * (mt**2 + 2.0 * mw**2) * \
       (1.0 - mw**2/mt**2)**2
# Branching ratio for t-> H+,b
def brtHb(X,Y):
       return lamtHb(X,Y) /(lamtHb(X,Y) +lamtWb)
# Branching ratio for t-> w,b
def brtWb(X,Y):
       return lamtWb/(lamtHb(X,Y)+lamtWb)

###############################################################
### Calculate the c,b events ;tau,nv events; c,s events created from H+ in t-> H+,b 
#c,b events
def cbevents(X,Y,Z):
       return brtHb(X,Y) * brcb(X,Y,Z)
#c,s events
def csevents(X,Y,Z):
       return brtHb(X,Y) * brcs(X,Y,Z)
#tau,nv events
def tnevents(X,Y,Z):
       return brtHb(X,Y) * brtn(X,Y,Z)
#c,b events + c,s events
def totalcbcs(X,Y,Z): 
       return cbevents(X,Y,Z) + csevents(X,Y,Z)
################################################
############################################################
# |X| value for theta, tangentbeta, tangamma, delta (only care about first H+ this time)
# Condition d = 1 u = 2 e =3
def X_function(i,j,k,l):
        sga = k / math.sqrt(k**2 + 1.0) # sinegamma
        cga = 1.0 / math.sqrt(1.0 + k**2)# cosgamma
        sthe = math.sin(i)#sinetheta
        cthe = math.cos(i)#costheta
        cde = math.cos(l)#cosdelta
        sde = math.sin(l)#sindelta
        return complex((- cthe * j * cde - sthe * cga)/sga , (- cthe * j * sde/sga))
#############################################################################
# |Y| value for theta, tangentbeta, tangamma, delta (only care about first H+ this time)
# Condition d = 1 u = 2 e =3
def Y_function(i,j,k,l):
        sga = k / math.sqrt(k**2 + 1.0) # sinegamma
        cga = 1.0 / math.sqrt(1.0 + k**2)# cosgamma
        sthe = math.sin(i)#sinetheta
        cthe = math.cos(i)#costheta
        cde = math.cos(l)#cosdelta
        sde = math.sin(l)#sindelta
        return complex((- cthe * cde / j + sthe * cga) / sga , - cthe * sde /(j * sga))
#############################################################################   
# |Z| value for theta, tangentbeta, tangamma, delta (only care about first H+ this time)
# Condition d = 1 u = 2 e =3
def Z_function(i,k):
        sthe = math.sin(i)#sinetheta
        return sthe * k
#############################################################################
#unitary matrix (diagonalization matrix)
#### theta, tangentbeta, tangamma, delta for U^dagger
def U(i,j,k,l):
        sga = k / math.sqrt(k**2 + 1.0) # sinegamma
        cga = 1.0 / math.sqrt(1.0 + k**2)# cosgamma
        sthe = math.sin(i)#sinetheta
        cthe = math.cos(i)#costheta
        cde = math.cos(l)#cosdelta
        sde = math.sin(l)#sindelta
        cbe = 1.0 / math.sqrt(1.0 + j**2)# cosbeta
        sbe = j / math.sqrt(j**2 + 1.0) # sinebeta
        ud1 = sga * cbe
        ud2 = complex(- cthe * sbe * cde  - sthe * cga * cbe, cthe * sbe *sde)
        ud3 = complex(sthe * sbe * cde - cthe * cga * cbe, sthe * sbe * sde)
        uu1 = sga * sbe
        uu2 = complex(cthe * cbe * cde - sthe * cga * sbe, cthe * cbe * sde)
        uu3 = complex(- sthe * cbe * cde - cthe * cga * sbe, - sthe * cbe * sde)
        ul1 = cga
        ul2 = sthe * sga
        ul3 = cthe * sga
        return [[ud1,ud2,ud3],[uu1,uu2,uu3],[ul1,ul2,ul3]] 
###################################################################
def start():#choose what 2 parameters in total 4 parameters (theta,tanbeta,tangamma,delta)
        global read1,A,B,read2
        while True:
            read1 = input('x-axis:number 0 to 3 (0 for theta, 1 for tanbeta,2 for tangamma,3 for delta):')
            if read1 == str(0) :
                 A = the 
                 print("{}".format("theta"))
                 print('A',A)
                 break
            elif read1 == str(1):
                 A = tbe
                 print("{}".format("tanbeta"))
                 print('A',A)
                 break
            elif read1 == str(2):
                 A = tga                 
                 print("{}".format("tangamma"))
                 print('A',A)
                 break
            elif read1 == str(3):
                 A = delta                 
                 print("{}".format("delta")) 
                 print('A',A)
                 break                 
            else :
                 print("Not correct variable")
    
        while True:
            read2 = input('y-axis:number 0 to 3 (0 for theta, 1 for tanbeta,2 for tangamma,3 for delta):')
            if read2 != read1:
              if read2 == str(0):
                  B = the
                  print("{}".format("theta"))
                  print('B',B)
                  break
              elif read2 == str(1):
                  B = tbe
                  print("{}".format("tanbeta"))
                  print('B',B)
                  break
              elif read2 == str(2):
                  B = tga
                  print("{}".format("tangamma"))
                  print('B',B)
                  break
              elif read2 == str(3):
                  B = delta
                  print("{}".format("delta")) 
                  print('B',B)
                  break
              else:
                  print("Not correct variable")
        else :
              print("this is already used")      
        return

#############################################################################
def start1():# choose model
        global X2,Y2,Z2
        while True:
            read0 = input('Choose type of3HDM (1 for I, 2 for II ,\
3 for Leptonic-specific, 4 for flipped,5 for Democratic):')
            if read0 == str(1) :
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Z2
                 print('Model:',read0)
                 break
            elif read0 == str(2):
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[1][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #Z2
                 print('Model:',read0)
                 break
            elif read0 == str(3):
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #Z2
                 print('Model:',read0)
                 break
            elif read0 == str(4):
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Z2
                 print('Model:',read0)
                 break   
            elif read0 == str(5):
                 def X2(i,j,k,l):
                     return U(i,j,k,l)[0][1] / U(i,j,k,l)[0][0] #X2
                 def Y2(i,j,k,l):
                     return - U(i,j,k,l)[1][1] / U(i,j,k,l)[1][0] #Y2
                 def Z2(i,j,k,l) :
                     return U(i,j,k,l)[2][1] / U(i,j,k,l)[2][0] #Z2
                 print('Model:',read0)
                 break 
              
            else :
                 print("Not correct type 3HDM")
            print('Model:',read0)
        return 
#############################################################################
def complexyfunction(i,j,k,l):# (XY^*) function 
        return (X2(i,j,k,l) * (np.conj(Y2(i,j,k,l))))
###############################################################################
def start3():
#        global X2,Y2,Z2 #
#        print('x2',X2)
        reference_array = [i,j,k,l]
        read1_int = int(read1)
        read2_int = int(read2)
        for var_b in B:
            for var_a in A:
                my_tuple = ()
                for counter in range(4):
                    if read1_int == counter:
                        my_tuple += (var_a,)
                    elif read2_int == counter:
                        my_tuple += (var_b,)
                    else:
                        my_tuple += (reference_array[counter],)
                absolutexlist.append(abs(X2(*my_tuple))) # unpacked bracket of my_tuple by * sign
                absoluteylist.append(abs(Y2(*my_tuple)))
                absolutezlist.append(abs(Z2(*my_tuple)))
                xyfun.append(complexyfunction(*my_tuple))
                BRCBfinal.append(brcb(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple))))
#                BRCB2final.append(brcb(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple)))**2)
                BRCSfinal.append(brcs(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple))))
                BRTNfinal.append(brtn(abs(X2(*my_tuple)),abs(Y2(*my_tuple)),abs(Z2(*my_tuple))))
        print('lenbrcb',len(BRCBfinal))
        BRCBfinal_list.append(BRCBfinal)
        BRCSfinal_list.append(BRCSfinal)
        BRTNfinal_list.append(BRTNfinal)
def max_valueposition(xxx):# x has to be np.array ; The max value postion of long array
    yyy = np.where(xxx == max(xxx))[0]
    return int(yyy)
#######################################################################
#######################################################################
#####################################################
#plot labels 
start()
start1()
readlist = ("$\\theta$","tan$\\beta$","tan$\\gamma$","$\\delta$") 
for n in np.arange(0,len(cepc.mhch_list_CEPC)):
    mhch = cepc.mhch_list_CEPC[n]
    print(mhch)
    alpmhch = alpmz /(1.0 + (alpmz/(12.0 * PI)) * (33.0 - 2 * nf) * \
                  math.log((mhch**2)/(mz**2))) 

    alpmb= alpmz /(1.0 + (alpmz/(12.0 * PI))*(33.0 - 2 * nf) * \
               math.log(4.89**2/(mz**2)))
    print('alpmb :',alpmb, ' alpmhch:',alpmhch)

# alpha_s at sale 2 GeV.Can read this off from Rev.Part.Properties
    alp2gev = 0.117 /(1.0 + (0.117 /(12.0 * PI)) * (33.0 - 2 * nf) *\
                  math.log(2.0**2 / mz**2))
#####################################################
    alp2gev = 0.310 # Take this value will have same figure as Andrew and Stefano' paper
#    [ PHYSICAL REVIEW D 85, 115002 (2012) ]
#####################################################
# for running of ms at 2 GeV to ms at 4.88 GeV (mb pole mass)                                                                                                                                                                                                                          
    fc2gev=1.0 + 1.101 * alp2gev/PI + 1.39 * (alp2gev/PI)**2 + 1.09 * \
     (alp2gev/PI)**3
    fcmb=1.0 + 1.101 * alpmb/PI + 1.39 * (alpmb/PI)**2 + 1.09 * \
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
    print('ms_value',ms,mb)
#mb = 2.95# mass bottom quark
#ms =0.055
    print('MC MB MS',mc , mb, ms)
###################################################################################                                                                                                                                                                                                  
    if mhch > ma:
      lamaw=9.0 * gf**2 * mw**4 * mhch * gaw /(16.0 * PI**3)
    else:
      lamaw=0.00
###################################################################################
# Produce brancing ratio results
    absolutexlist = []
    absoluteylist = []
    absolutezlist = []
    xyfun = []
    BRCBfinal = []
    BRCSfinal = []
    BRTNfinal = []
    BRCBfinal_list = []
    BRCSfinal_list = []
    BRTNfinal_list = []
    start3()
    BRCB_TN = np.array(BRCBfinal) * np.array(BRTNfinal)# product of BRCB and BRTN
    BRCB_CS = np.array(BRCBfinal) * np.array(BRCSfinal)# product of BRCB and BRCS
    BRCBPLUSCS = np.array(BRCBfinal) + np.array(BRCSfinal)# Sum of BRCB and BRCS
#### 2 jet tagged plot array
    twotagsig = cepc.CEPCevent[n] * (fl.eeHHcbtn_1bsignal(BRCBfinal,BRTNfinal) +\
                fl.eeHHcbtn_0bsignal(BRCBfinal,BRTNfinal) + \
                fl.eeHHcstn_1bsignal(BRCSfinal,BRTNfinal) + fl.eeHHcstn_0bsignal(BRCSfinal,BRTNfinal) \
                ) * fl.selection_2j 
#                              / np.sqrt(fl.backgroundtagging2())
#### 4 jet 2b-tagged plot array
    fourtag2bsig = cepc.CEPCevent[n] * (fl.eeHHcbcb_2bsignal(BRCBfinal,BRCBfinal) +\
                   fl.eeHHcbcb_1bsignal(BRCBfinal,BRCBfinal) + \
                   fl.eeHHcbcs_2bsignal(BRCBfinal,BRCSfinal) + fl.eeHHcbcs_1bsignal(BRCBfinal,BRCSfinal) + \
                   fl.eeHHcbcb_0bsignal(BRCBfinal,BRCBfinal) + fl.eeHHcbcs_0bsignal(BRCBfinal,BRCSfinal) + \
                   fl.eeHHcscs_0bsignal(BRCSfinal,BRCSfinal)) 
#                           \   * fl.epsilon / np.sqrt(fl.backgroundtagging())
#### 4 jet 1b-tagged plot array
    fourtag1bsig = cepc.CEPCevent[n] * (fl.real_b_cbcb(BRCBfinal,BRCBfinal) +\
                   fl.fake_b_cbcb(BRCBfinal,BRCBfinal) + fl.real_b_cbcs(BRCBfinal,BRCSfinal) + \
                   fl.fake_b_cbcs(BRCBfinal,BRCSfinal) + fl.real_b_cscs(BRCSfinal,BRCSfinal) + \
                   fl.fake_b_cscs(BRCSfinal,BRCSfinal)) 
#                             \ * fl.epsilon / np.sqrt(fl.backgroundtagging())
    ###############################################
    #PLOT OF CONTOUR WITH RANGE OS 4 PARAMETERS INTO |X|,|Y|,|Z = 0.1| WITH CS,CB,TAUNV
# (4 parameters):A,B, BRCS contour plot [in the :reshape(y,x) not reshape(x,y)]
    plt.figure(1)
    plt.subplot(221)
    Contourbrcs = plt.contourf(A,B, \
          np.resize(BRCSfinal,len(BRCSfinal)).reshape(len(B),len(A)),\
          colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
    plt.axis([min(A), max(A),min(B), max(B)])
#plt.clabel(Contourbrcs, inline= 0.01, fontsize=10)# contour level show
    plt.colorbar(Contourbrcs)
    plt.title('BR($H^{\pm} \longrightarrow $ cs)')#,$M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,cs.png')
#    plt.show()
#    plt.close()
#    print('--------------------------------------------------')
# (4 parameters):A,B, BRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
    plt.subplot(222)
    Contourbrcb = plt.contourf(A,B, \
           np.resize(BRCBfinal,len(BRCBfinal)).reshape(len(B),len(A)),\
           colors = ['black','royalblue','purple','darkgreen','brown','red','black'])
#    plt.clabel(Contourbrcb, inline= 0.02, fontsize= 9)# contour level show
    plt.colorbar(Contourbrcb)
    plt.title('BR($H^{\pm} \longrightarrow $ cb)')#,$M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,cb.png')
#    plt.show()
#    plt.close()
#    print('--------------------------------------------------')
#(4 parameters): A,B, BRCB contour plot [in the :reshape(y,x) not reshape(x,y)]
    plt.subplot(223)
    Contourbrtn = plt.contourf(A,B, \
           np.resize(BRTNfinal,len(BRTNfinal)).reshape(len(B),len(A)),\
           colors = ['black','royalblue','purple','darkgreen','brown','red','black'],levels = linecs)
#plt.clabel(Contourbrtn)# contour level show
    plt.colorbar(Contourbrtn)
    plt.title('BR($H^{\pm} \longrightarrow $ $\\tau \\nu_\\tau $)')#,$M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#   plt.savefig('M{H^{\pm}}= '+ str(mhch) +' GeV,TN.png')
#    plt.show()
#    plt.close()
#    print('--------------------------------------------------')
#(4 parameters): A,B, BRCB + BRCS contour plot 
    plt.subplot(224)
    ContourBRCBCS = plt.contourf(A,B, \
           np.resize(BRCBPLUSCS ,len(BRCBPLUSCS )).reshape(len(B),len(A)),\
           colors = ['black','royalblue','purple','yellow','brown','red','gray','green'])#, levels = np.arange(0.06, 0.16,0.02))
#    plt.clabel(ContourBRCBCS, inline= 0.02, fontsize= 9)# contour level show
    plt.colorbar(ContourBRCBCS)
    plt.title(' BR($H^{\pm} \longrightarrow $ cb) + BR($H^{\pm} \longrightarrow $ cs) ')#, $M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('CharHCBplusCS_CEPC'+ str(mhch) +'.png')
    plt.subplots_adjust(top=0.99, bottom=0.06, left=0.1, right=0.99, hspace= 0.6,
                    wspace=0.3)
    plt.show()
    plt.close()
# (4 parameters):A,B, 4jet tagged plots
    plt.figure()
    signal4jet = plt.contourf(A,B, \
        np.resize( fourtag2bsig ,\
              len( fourtag2bsig)).\
        reshape(len(B),len(A)),\
        colors = ['black','royalblue','purple','yellow','brown','red','gray','green'])# ,levels = np.arange(1.0,6.0,1.0))
    plt.title('Signal of $H^{\pm}$ 2btagged 4jet '+\
             ', $M_{H^{\pm}}$= '+ str(mhch) +' GeV')
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.colorbar(signal4jet)
    plt.savefig('sig_4jet_CEPC'+ str(mhch) +'.png')
    plt.show()
    plt.close()
###############################################
# (4 parameters):A,B, 4jet tagged plots 1-b
    plt.figure()
    signaloneb4jet = plt.contourf(A,B, \
        np.resize( fourtag1bsig ,\
              len( fourtag1bsig)).\
        reshape(len(B),len(A)),\
        colors = ['black','royalblue','purple','yellow','brown','red','gray','green'])# ,levels = np.arange(1.0,6.0,1.0))
    plt.title('Signal of $H^{\pm}$ 1btagged 4jet '+\
             ', $M_{H^{\pm}}$= '+ str(mhch) +' GeV')
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.colorbar(signaloneb4jet)
    plt.savefig('sig_1b4jet_CEPC'+ str(mhch) +'.png')
    plt.show()
    plt.close()
################################
#(4 parameters):A,B,4jet notagging plot (BR($H^{\pm} \longrightarrow $ cb + cs))
    plt.figure()
    signalnotag4jet = plt.contourf(A,B, \
        np.resize(cepc.CEPCevent[n] * (BRCBPLUSCS**2),\
              len(cepc.CEPCevent[n] * (BRCBPLUSCS**2))).\
        reshape(len(B),len(A)),colors = ['black','royalblue','purple','yellow','brown','red','gray','green'])
    plt.colorbar(signalnotag4jet)
    plt.title('Signal 4jet notag,$M_{H^{\pm}}$= '+ str(mhch) +' GeV')#plot title
    plt.xlabel(readlist[int(read1)])# x-axis label
    plt.ylabel(readlist[int(read2)])# y-axis label
#    plt.xscale('log')
#    plt.yscale('log')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.savefig('sig_4jetnotag_CEPC'+ str(mhch) +'.png')
    plt.show()
    plt.close()