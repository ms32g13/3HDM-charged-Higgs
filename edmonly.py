#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:27:11 2019

@author: muyuansong0419
b_sgamma.py plot_function
"""
from invariant import *
import bsgamma as bsg
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as axe
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import exercise as exe
from neutronedm import dC_btH, dn_CEDM,CW,dn,de
import pandas as pd
import random
#from exe import xyfun,xyfun3,yfun,yfun3,U,X2,X3,Y2,Y3,Z2,Z3,\
#        complexyfunction,complexyfunction3,read1,read2,A,B,readlist
#print('A,B',exe.A,exe.B,exe.read1,exe.read2,len(exe.B),type(exe.read1))
#Plot Functions
def ABarray4(): # using [A,B] to plot BR(B_bar > X_s +gamma)
    i = exe.i # theta
    j = exe.j  # tanbeta
    k = exe.k #tangamma
    l = exe.l #delta
#    print('i,j,k,l',i,j,k,l)
    longlist = []
    reference_array = [i,j,k,l]
    for var_b in exe.B:
            for var_a in exe.A:
                my_tuple = ()
                for counter in range(4):
                    if int(exe.read1) == counter:
                        my_tuple += (var_a,)
                    elif int(exe.read2) == counter:
                        my_tuple += (var_b,)
                    else:
                        my_tuple += (reference_array[counter],)
                longlist.append(my_tuple)
    return longlist
##################################################################################
def plt_A_B_bsgnedm(i,j):# Bsgamma-result and N-EDM in {A,B} plane
    mass_axis1,mass_axis2 = i,j
    print(ABarray4()[0],mass_axis1,len(ABarray4()))
    print(ABarray4()[1],mass_axis2,len(ABarray4()))
    resultb = []
    resultn = []
    resulte = []
#B>Xs+gamma SECTION    
    for n in np.arange(0,len(ABarray4()) ):
        y3hdm= bsg.BR_B_Xs_gamma(mb,mw,mass_axis1,mass_axis2,\
                        exe.Y2(*ABarray4()[n] ), exe.complexyfunction(*ABarray4()[n] ),\
                        exe.Y3(*ABarray4()[n] ), exe.complexyfunction3(*ABarray4()[n] )) 
        resultb.append(y3hdm / (1e-4) )
#Nedm SECTION    
        nedm3hdm = abs(dn(mass_axis1,mass_axis2, exe.complexyfunction(*ABarray4()[n]),\
                    exe.complexyfunction3(*ABarray4()[n]) ) / (5.06e13)  )\
                        / 1e-26
        resultn.append( nedm3hdm )
#eedm SECTION
        eedm3hdm = abs(de(mass_axis1,mass_axis2,exe.yconjz2(*ABarray4()[n]),\
                    exe.yconjz3(*ABarray4()[n]) ) /1e-29  )
        resulte.append( eedm3hdm )
#########
    ned = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultn).flatten()  ,len(np.array(resultn).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,\
          levels = np.array([0.0,1.8]),colors = ['red']  )
#########
    bsgamm = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultb).flatten()  ,len(np.array(resultb).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,\
          levels = np.array([2.99,3.55]),colors = ['green'] )
#########
#    eed = plt.contourf(exe.A, exe.B, \
#           np.resize(np.array(resulte).flatten()  ,len(np.array(resulte).flatten() ) ).\
#          reshape(len(exe.B),len(exe.A)) ,\
#          levels = np.array([0.0,1.1]),colors = ['blue']  )
    
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) and NEDM in '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2) )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
#    plt.axis([0,60,-1.6,0]) #{tanbeta/tangamma,theta} plane
#    plt.axis([0,2 * PI ,-1.6,0]) #{theta,delta} plane
    plt.axis([0,60,0,60])# {tanbeta,tangamma} plane
    plt.savefig( str("%02d" % mass_axis1) + str("%02d"% mass_axis2) +'bsg.png')
    plt.show()
    plt.close()
def plt_A_b_eedm(i,j):#E-EDM in {A,B} plane
    mass_axis1,mass_axis2 = i,j
    resulte = []
    for n in np.arange(0,len(ABarray4()) ):
#        print('1,',exe.yconjz2(*ABarray4()[n]))
        eedm3hdm = abs(de(mass_axis1,mass_axis2,exe.yconjz2(*ABarray4()[n]),exe.yconjz3(*ABarray4()[n]) ) )
        resulte.append( eedm3hdm )
    eed = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resulte).flatten()  ,len(np.array(resulte).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,\
          levels = np.array([0.0,1.1e-29]),colors = ['blue']  )
#    plt.colorbar(eed)
    plt.title('E-EDM in '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2) )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
#    plt.axis([0,2 * PI ,-1.6,0]) #{theta,delta} plane
    plt.axis([0,60,0,60])# {tanbeta,tangamma} plane
#    plt.savefig( str("%02d" % mass_axis1) + str("%02d"% mass_axis2) +'eedm.png')
    plt.show()
    plt.close()
###############################################################################
###############################################################################
###############################################################################
def NEUTRONEDMtext():#save Neutron edm result in txt file
    
    f = open('Neutron_EDM.txt','w')
    f.write('%s %s %s %s %s %s\n' % ("   N","   MH+1", "    MH+2",\
                                 "  X_1Y_1*", "   X_2Y_2*","  Result") )
    for n in range(len(ABarray4()) ):
#    print('dn',n, dn(80,170,exe.complexyfunction(*ABarray4()[n]).imag, \
#                       exe.complexyfunction3(*ABarray4()[n]).imag) )
        f.write( "%5.0f %5.1f %5.1f %10.10e %10.10e %10.15e\n" % (n,80,170,\
            - exe.complexyfunction(*ABarray4()[n]),\
            - exe.complexyfunction3(*ABarray4()[n]), \
            dn(80,170,- exe.complexyfunction(*ABarray4()[n]), \
            - exe.complexyfunction3(*ABarray4()[n])) ) )
    f.close()
    return
def NEDMfigure4_plot():#NEDM figure 4
    ################### First charged Higgs
    masss = 120 # first charged Higgs
    ################### Second charged Higgs
    masss2 = 300
    #################### mass different
    massdiffer = 20
    massdiffer_list = np.array([ i for i in np.arange(20,80,20)] )
    massdiffer_list2 = np.array([ i for i in np.arange(0,355,50)] )
    m_axis = np.array([ i for i in np.arange(0,550,50)] )
    m1_axis = np.array([ i for i in np.arange(0,1021,1)] )
    xy_axis = np.array([ i for i in np.arange(0,1.1,0.1)] )
    xy1_axis = np.array([ i for i in np.arange(-5,5.1,0.1)] )
    longarray = np.array([i for i in np.arange(-1,1.2,0.2)])
    arrayxy1xy2 =  np.array([42,12,20,40])# Maximum of Im(XY^*) can get
    empty = []
    emptye =[]
    for j in xy1_axis:
        for i in m1_axis:
            # 1GeV =  5.06e13 cm^{-1}
            nedm = abs(dn(i,300,[complex(0,j)],[complex(0,0)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
            eedm = abs(de(i,300,[complex(0,j)],[complex(0,0)]))
            empty.append(nedm )
            emptye.append(eedm)
#            print('nedm',i,j,nedm )
    result = plt.contourf(m1_axis, xy1_axis, \
           np.resize(np.array(empty),len(np.array(empty) )).\
           reshape(len(xy1_axis),len(m1_axis)), \
           levels = np.array([0.0,1.8e-26])      
           )
    plt.colorbar(result)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('IM($X_2Y_2^*$)')
    plt.title('NEDM in 2HDM')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([100,500, 0, 1.0])
    plt.show()
    plt.close()
    resultedm = plt.contourf(m1_axis, xy1_axis, \
           np.resize(np.array(emptye),len(np.array(emptye) )).\
           reshape(len(xy1_axis),len(m1_axis)), \
           levels = np.array([0.0,1.1e-29])      
           )
    plt.colorbar(resultedm)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('IM($X_2^*Z_2$)')
    plt.title('e-EDM in 2HDM')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([100,500, 0, 0.005])
    plt.show()
    plt.close()
#################
    empty3hdm = []
    empty3hdme = []
    for j in xy1_axis:
        for i in m1_axis:
            # 1GeV =  5.06e13 cm^{-1}
                nedm = abs(dn(i,masss2,[complex(0,j)],[ complex(0,- j)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
                eedm2 = abs(de(i,masss2,[complex(0,j)],[complex(0,-j)]))
                empty3hdm.append(nedm )
                empty3hdme.append(eedm2)
    result8 = plt.contourf(m1_axis, xy1_axis, \
                          np.resize(np.array(empty3hdm),len(np.array(empty3hdm) )).\
                          reshape(len(xy1_axis),len(m1_axis)), \
                          levels = np.array([0.0,1.8e-26])      
                          )
    plt.colorbar(result8)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('IM($X_2Y_2^*$)')
    plt.title('NEDM in 3HDM with $M_{H^{\pm}_2} = $' + str(masss2))
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([0,500, 0, 1])
    plt.show()
#        plt.savefig('m1' + str(masss2) + 'imxy2'+ str("%.2f" % m) + '.png' )
    plt.close()
    result8e = plt.contourf(m1_axis, xy1_axis, \
                          np.resize(np.array(empty3hdme),len(np.array(empty3hdme) )).\
                          reshape(len(xy1_axis),len(m1_axis)), \
                          levels = np.array([0.0,1.1e-29])      
                          )
    plt.colorbar(result8e)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('IM($X_2^*Z_2$)')
    plt.title('e-EDM in 3HDM with $M_{H^{\pm}_2} = $' + str(masss2))
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([0,500, 0, 0.1])
    plt.show()
#    plt.savefig('m1' + str(masss2) + 'imxy2'+ str("%.2f" % m) + '.png' )
    plt.close()
##################
    
    for m in arrayxy1xy2:
        threehdm = []
        threehdme = []
        for j in massdiffer_list2:
            for i in m1_axis:
             # 1GeV =  5.06e13 cm^{-1}
                nedm00 = abs(dn(i,i + j,[complex(0,m)],[ complex(0,- m)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
                eedm00 = abs(de(i,i + j,[complex(0,m)],[complex(0, - m)]))
                threehdm.append(nedm00 )
                threehdme.append(eedm00)
        result9 = plt.contourf(m1_axis, massdiffer_list2, \
                          np.resize(np.array(threehdm),len(np.array(threehdm) )).\
                          reshape(len(massdiffer_list2),len(m1_axis)), \
                          levels = np.array([0.0,1.8e-26])      
                          )
        plt.colorbar(result9)
        plt.xlabel('$M_{H^{\pm}_{1}}$')
        plt.ylabel('Mass difference')
#        plt.title('n-EDM with IM($X_2Y_2^*$) = $\sqrt{2}$, IM($X_3Y_3^*$) = - $\sqrt{2}$' )
        plt.title('n-EDM with IM($X_2Y_2^*$) = '+ str('%.2g'% m) + \
              ', IM($X_3Y_3^*$) = - ' + str('%.2g'% m) )
        plt.axis([0,500,0, 20])
        plt.show()
#        plt.savefig('m1masdifer'+ 'nedm.png' )
        plt.close()
        result9e = plt.contourf(m1_axis, massdiffer_list2, \
                          np.resize(np.array(threehdme),len(np.array(threehdme) )).\
                          reshape(len(massdiffer_list2),len(m1_axis)), \
                          levels = np.array([0.0,1.1e-29])      
                          )

        
        plt.colorbar(result9e)
        plt.xlabel('$M_{H^{\pm}_{1}}$')
        plt.ylabel('Mass difference')
        plt.title('e-EDM with IM($X_2^*Z_2$) = '+ str('%.2g'% m) + \
              ', IM($X_3^*Z_3$) =  ' + str('%.2g'% - m) )
        plt.savefig('m1masdiferimxz2'+ str('%.2g'% m) + 'imxz3'+ str('%.2g'% m) + 'eedm.png' )
        plt.axis([0,500,0, 1])
        plt.show()
        plt.close()
#####################################
   # for m in m1_axis:
    empty20 = []
    empty40 = []
    empty60 = []
    empty80 = []
    masssdif = [20,40,60,80]
    xx_ary = np.array([0.35,0.45])
    yy_ary = np.array([0.80,0.65,0.50,0.35])
    for j in longarray:
        for i in longarray:
                    # 1GeV =  5.06e13 cm^{-1}
                    nedm20 = abs(dn(masss,masss +  masssdif [0],[complex(0,i)],\
                                [complex(0,j)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
                    nedm40 = abs(dn(masss,masss + masssdif [1],[complex(0,i)],\
                                [complex(0,j)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
                    nedm60 = abs(dn(masss,masss + masssdif [2],[complex(0,i)],\
                                [complex(0,j)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
                    nedm80 = abs(dn(masss,masss + masssdif [3],[complex(0,i)],\
                                [complex(0,j)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
                    empty20.append(nedm20 )
                    empty40.append(nedm40 )
                    empty60.append(nedm60 )
                    empty80.append(nedm80 )
    result20 = plt.contour(longarray, longarray, \
                                   np.resize(np.array(empty20),len(np.array(empty20) )).\
                                   reshape(len(longarray),len(longarray)), \
                                   levels = np.array([0.0,1.8e-26]),      
                                   colors='green',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 20')
    result40 = plt.contour(longarray, longarray, \
                                   np.resize(np.array(empty40),len(np.array(empty40) )).\
                                   reshape(len(longarray),len(longarray)), \
                                   levels = np.array([0.0,1.8e-26]),      
                                   colors='black',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 40')
    result60 = plt.contour(longarray, longarray, \
                                   np.resize(np.array(empty60),len(np.array(empty60) )).\
                                   reshape(len(longarray),len(longarray)), \
                                   levels = np.array([0.0,1.8e-26]),      
                                   colors='blue',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 60')
    result80 = plt.contour(longarray, longarray, \
                                   np.resize(np.array(empty80),len(np.array(empty80) )).\
                                   reshape(len(longarray),len(longarray)), \
                                   levels = np.array([0.0,1.8e-26]),      
                                   colors='red',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 80')
    
    plt.xlabel('IM($X_2Y_2^*$)')
    plt.ylabel('IM($X_3Y_3^*$)')
    plt.plot(xx_ary, np.array([yy_ary[0]] * len(xx_ary)),color='green')
    plt.plot(xx_ary, np.array([yy_ary[1]] * len(xx_ary)),color='black')
    plt.plot(xx_ary, np.array([yy_ary[2]] * len(xx_ary)),color='blue')
    plt.plot(xx_ary, np.array([yy_ary[3]] * len(xx_ary)),color='red')
    plt.text(0.5, yy_ary[0] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [0]) )
    plt.text(0.5, yy_ary[1] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [1]) )
    plt.text(0.5, yy_ary[2] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [2]) )
    plt.text(0.5, yy_ary[3] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [3]) )
    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) )  
#    plt.title('NEDM with $M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +'+ str("%3d" % mn)+ 'GeV')
#   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([-1.0,1.0,-1.0, 1.0])
#    plt.savefig(str("%3d" % m) + '_' + str("%3d" % (mn + m) ) \
#                        + '_imxy2_' + '_imxy3_'+ '.png' )
    plt.show()
    plt.close()
    
    empty3 = []
    mass90 = 90
    for j in xy_axis:
        for i in xy_axis:
            # 1GeV =  5.06e13 cm^{-1}
            nedm3 = abs(dn(mass90,mass90 + massdiffer,[complex(1,i)],[complex(1,j)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
            
            empty3.append(nedm3 )
#            print('nedm',i,j,nedm )
    result3 = plt.contourf(xy_axis, xy_axis, \
           np.resize(np.array(empty3),len(np.array(empty3) )).\
           reshape(len(xy_axis),len(xy_axis)), \
           levels = np.array([0,1.8e-26])      
           )
    plt.colorbar(result3)
    plt.xlabel('IM($X_2Y_2^*$)')
    plt.ylabel('IM($X_3Y_3^*$)')
    plt.title('NEDM in 3HDM' )
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([-1.0, 1.0,-1.0, 1.0])
#    plt.show()
    plt.close()
    
    empty4 = []
    for j in m1_axis:
        for i in m1_axis:
            # 1GeV =  5.06e13 cm^{-1}
            nedm4 = abs(dn(i,j,[complex(0,42)],[complex(0,- 42)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
            
            empty4.append(nedm4 )
#            print('nedm',i,j,nedm )
    result4 = plt.contourf(m1_axis, m1_axis, \
           np.resize(np.array(empty4),len(np.array(empty4) )).\
           reshape(len(m1_axis),len(m1_axis)), \
           levels = np.array([0.0,1.8e-26])      
           )
    plt.colorbar(result4)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('NEDM with IM($X_2Y_2^*$) = 42, IM($X_3Y_3^*$) = - 42')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([0,1000,0,1000])
    plt.show()
    plt.close()
    return

def nedm3hdm_plot():#NEDM figure 4 3hdm
    masss = 100 # first charged Higgs
    masss2 = 200 # first charged Higgs
    for i in np.arange(20,100,20):
        resultxy2imag = []
        resultxy3imag = []
        for n in np.arange(0,len(ABarray4()) ):
            resultxy2imag.append( (  exe.complnedm2(*ABarray4()[n] ))  )
            resultxy3imag.append( (  exe.complnedm3(*ABarray4()[n] ))  )
        nedm = abs(dn(masss,i+ masss,resultxy2imag,resultxy3imag) / (5.06e13) )/ 1e-26# dn [GeV^{-1}] > # cm
        nedm2 = abs(dn(masss2,i + masss2,resultxy2imag,resultxy3imag) / (5.06e13) )/ 1e-26# dn [GeV^{-1}] > # cm
#    print('empty',len(nedm),np.array(nedm).flatten())
        result = plt.contourf(exe.A, exe.B, \
                         np.resize(np.array(nedm).flatten(),len(np.array(nedm).flatten())).\
                                  reshape(len(exe.B),len(exe.A)), \
                                 levels = np.array([0.0,1.8])        
                                  )
        plt.colorbar(result)
        plt.xlabel(exe.readlist[int(exe.read1)])
        plt.ylabel(exe.readlist[int(exe.read2)])
        plt.title('NEDM in 3HDM,'+ str(masss)+','+ str(i+ masss))
        plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
        plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
        plt.show()
        plt.close
        result2 = plt.contourf(exe.A, exe.B, \
                         np.resize(np.array(nedm2).flatten(),len(np.array(nedm2).flatten())).\
                                  reshape(len(exe.B),len(exe.A)), \
                                 levels = np.array([0.0,1.8])        
                                  )
        plt.colorbar(result2)
        plt.xlabel(exe.readlist[int(exe.read1)])
        plt.ylabel(exe.readlist[int(exe.read2)])
        plt.title('NEDM in 3HDM,'+ str(masss2) +','+ str(i+ masss2))
        plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
        plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([80,500, 0.0, 0.5])
        plt.show()
        plt.close
    return
#######################################
def nedm3hdm_xy1_xy2():# NEDM figure 4 3hdm
    ################### First charged Higgs
    masss = 200 # first charged Higgs
    masssdif = [20,40,60,80]#mass differ
    xx_ary = np.array([0.1,0.15])
    yy_ary = np.array([0.48,0.40,0.32,0.24])
    resultxy2imag = []
    resultxy3imag = []
    empty20 = []
    empty40 = []
    empty60 = []
    empty80 = []
    for n in np.arange(0,len(ABarray4()) ):
        resultxy2imag.append( (  exe.complnedm2(*ABarray4()[n] ))  )
        resultxy3imag.append( (  exe.complnedm3(*ABarray4()[n] ))  )
#    print(np.array(resultxy2imag).flatten())
    xy1 = np.array(resultxy2imag).flatten()
    xy2 = np.array(resultxy3imag).flatten()
    for i in xy1:
        for j in xy2:
            nedm20 = abs(dn(masss,masss +  masssdif [0],[i],\
                                [j]) / (5.06e13) ) / 1e-26# dn [GeV^{-1}] > # cm
            nedm40 = abs(dn(masss,masss +  masssdif [1], [i],\
                                [j]) / (5.06e13) ) / 1e-26# dn [GeV^{-1}] > # cm
            nedm60 = abs(dn(masss,masss +  masssdif [2], [i],\
                                [j]) / (5.06e13) ) / 1e-26# dn [GeV^{-1}] > # cm
            nedm80 = abs(dn(masss,masss +  masssdif [3], [i],\
                                [j]) / (5.06e13) ) / 1e-26# dn [GeV^{-1}] > # cm
            empty20.append(nedm20 )
            empty40.append(nedm40 )
            empty60.append(nedm60 )
            empty80.append(nedm80 )
    result20 = plt.contour(np.imag(xy1),np.imag(xy2), \
               np.resize(empty20,len(empty20 )).\
               reshape(len(np.imag(xy2)),len(np.imag(xy1))), \
               levels = np.array([0.0,1.8]),      
               colors='green',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 20')
#    plt.colorbar(result20)
#    plt.text(4, yy_ary[0] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [0]) )
#    plt.text(0.5, yy_ary[1] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [1]) )
#    plt.text(0.5, yy_ary[2] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [2]) )
#    plt.text(0.5, yy_ary[3] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [3]) )
    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss))
#    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) \
#              + ', $M_{H^{\pm}_{2}}$ = '+ str("%3d" % (masss +  masssdif [0])) )  
#   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([-0.3,0.3,-0.3, 0.3])
#    plt.xlabel('IM($X_2Y_2^*$)')
#    plt.ylabel('IM($X_3Y_3^*$)')
#    plt.show()
#    plt.close()
    plt.plot(xx_ary, np.array([yy_ary[0] - 0.05] * len(xx_ary)),color='green')
    plt.plot(xx_ary, np.array([yy_ary[1] - 0.05] * len(xx_ary)),color='black')
    plt.plot(xx_ary, np.array([yy_ary[2] - 0.05] * len(xx_ary)),color='blue')
    plt.plot(xx_ary, np.array([yy_ary[3] - 0.05] * len(xx_ary)),color='red')
    plt.text(xx_ary[1]+0.05, yy_ary[0] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [0]) )
    plt.text(xx_ary[1]+0.05, yy_ary[1] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [1]) )
    plt.text(xx_ary[1]+0.05, yy_ary[2] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [2]) )
    plt.text(xx_ary[1]+0.05, yy_ary[3] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [3]) )
#####################
    result40 = plt.contour(np.imag(xy1),np.imag(xy2), \
               np.resize(np.array(empty40).flatten(),len(np.array(empty40).flatten() )).\
               reshape(len(np.imag(xy2)),len(np.imag(xy1))), \
               levels = np.array([0.0,1.8]),      
               colors='black',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 40')
#    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) \
#              + ', $M_{H^{\pm}_{2}}$ = '+ str("%3d" % (masss +  masssdif [1])) )
#    plt.xlabel('IM($X_2Y_2^*$)')
#    plt.ylabel('IM($X_3Y_3^*$)')
#    plt.show()
#    plt.close()
#################
    result60 = plt.contour(np.imag(xy1),np.imag(xy2), \
               np.resize(np.array(empty60).flatten(),len(np.array(empty60).flatten() )).\
               reshape(len(np.imag(xy2)),len(np.imag(xy1))), \
               levels = np.array([0.0,1.8]),      
               colors='blue',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 60')
#    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) \
#              + ', $M_{H^{\pm}_{2}}$ = '+ str("%3d" % (masss +  masssdif [2])) )
#    plt.xlabel('IM($X_2Y_2^*$)')
#    plt.ylabel('IM($X_3Y_3^*$)')
#    plt.show()
#    plt.close()
    result80 = plt.contour(np.imag(xy1),np.imag(xy2), \
               np.resize(np.array(empty80).flatten(),len(np.array(empty80).flatten() )).\
               reshape(len(np.imag(xy2)),len(np.imag(xy1))), \
               levels = np.array([0.0,1.8]),      
               colors='red',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 80')
#    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) \
#              + ', $M_{H^{\pm}_{2}}$ = '+ str("%3d" % (masss +  masssdif [3])) )
#    plt.axis([-2,2,-2, 1])
    plt.xlabel('IM($X_2Y_2^*$)')
    plt.ylabel('IM($X_3Y_3^*$)')
    plt.show()
    plt.close()
    return
def nedm_mhp1_massdiff():#[mhch1,mass_differ] in 3hdm nedm result
    mass_diff = np.arange(0,400,50)
    mhp1 = np.arange(0,520,20)
    x2y2 = 42.4 #np.sqrt(2)# 
    x3y3 = -42.4 #- np.sqrt(2)
    list_1 = []
    for j in mass_diff:
        for i in mhp1:
            nedm_0 = abs(dn(i,i + j,[complex(0, x2y2)],\
                                    [complex(0, x3y3)]) / (5.06e13)/ 1e-26 )# dn [GeV^{-1}] > # cm
            list_1.append(nedm_0)
    list_1plot = plt.contourf(mhp1,mass_diff, \
               np.resize(list_1,len(list_1 )).\
               reshape(len(mass_diff),len(mhp1)), \
               levels = np.array([0.0,1.8]))
    plt.colorbar(list_1plot)
    plt.axis([0,500,0, 20])
    plt.title('NEDM with IM($X_2Y_2^*$) =  42.4'+ \
              ', IM($X_3Y_3^*$) =  - 42.4')
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('Mass difference')
    plt.show()
    plt.close()
    return
def nedm_mhp1_mhp2():#[mhch1,mhch2] in 3hdm nedm result
    mhp1 = np.arange(0,1050,50)
    mhp2 = np.arange(0,1050,50)
    x2y2 = 42.4 # np.sqrt(2)
    x3y3 = - 42.4 #- np.sqrt(2)
    list_1 = []
    for j in mhp2:
        for i in mhp1:
            nedm_1 = abs(dn(i,j,[complex(0,x2y2)],\
                                    [complex(0,x3y3)]) / (5.06e13)/ 1e-26 )# dn [GeV^{-1}] > # cm
            list_1.append(nedm_1)
    list_2plot = plt.contourf(mhp1,mhp2, \
               np.resize(list_1,len(list_1 )).\
               reshape(len(mhp2),len(mhp1)), \
               levels = np.array([0.0,1.8]))
    plt.colorbar(list_2plot)
    plt.axis([0,1000,0, 1000])
    plt.title('NEDM with IM($X_2Y_2^*$) =  42.4'+ \
              ', IM($X_3Y_3^*$) = - 42.4')
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.show()
    plt.close()
    return
def nedm_3dparameter():#3D scan for 3 mixing parameters against Neutron EDM limit
    masss = 100 # first charged Higgs
    masss2 = 200 # first charged Higgs
    masssdif = [20,40,60,80]#mass differ
    list_a =[]
    list_b =[]
    list_c =[]
    n2 = []
    n3 = []
    result = ([],[],[],[])
    result2 = ([],[],[],[])
    print('result',result[0],type(result[0]),type(n2))
    for i in exe.the:
        for j in exe.tbe:
#            for k in exe.tga:
                for l in exe.delta:
                       my_tuple = (i,j,exe.k,l)  
                       two = exe.complnedm2(*my_tuple)
                       three = exe.complnedm3(*my_tuple)
                       list_a.append(i)
                       list_b.append(j)
                       list_c.append(l)
                      
                       n2.append(two)
                       n3.append(three)
                       
                       result[0].append(abs(dn(masss,masss +  masssdif [0],two,\
                                three) / (5.06e13) ) / 1e-26)
                       result[1].append(abs(dn(masss,masss +  masssdif [1],two,\
                                three) / (5.06e13) ) / 1e-26)
                       result[2].append(abs(dn(masss,masss +  masssdif [2],two,\
                                three) / (5.06e13) ) / 1e-26)
                       result[3].append(abs(dn(masss,masss +  masssdif [3],two,\
                                three) / (5.06e13) ) / 1e-26)
                       result2[0].append(abs(dn(masss2,masss2 +  masssdif [0],two,\
                                three) / (5.06e13) ) / 1e-26)
                       result2[1].append(abs(dn(masss2,masss2 +  masssdif [1],two,\
                                three) / (5.06e13) ) / 1e-26)
                       result2[2].append(abs(dn(masss2,masss2 +  masssdif [2],two,\
                                three) / (5.06e13) ) / 1e-26)
                       result2[3].append(abs(dn(masss2,masss2 +  masssdif [3],two,\
                                three) / (5.06e13) ) / 1e-26)
#    print(*np.array(result[0]) )
    print('imn2min,max',min(abs(np.array(n2).imag)), max(abs(np.array(n2).imag)) )
    print('imn2min,max',min(abs(np.array(n3).imag)), max(abs(np.array(n3).imag)) )
    print('ren2min,max',min(np.array(n2).real), max(np.array(n2).real))
    index = [i for i, j in enumerate(np.array(n2).real) if j == max(np.array(n2).real)]#position of max
    print('ren3min,max',min(np.array(n3).real), max(np.array(n3).real) )
#    print(len(exe.the),len(exe.tbe),len(exe.delta),len(result[0]),type(list_a),type(result[0]))
    d1 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'n2': [*n2]})
    d2 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'n3': [*n3]})
    d10 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'nedm': [*result[0]]})
    d11 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'nedm': [*result[1]]})
    d12 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'nedm': [*result[2]]})
    d13 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'nedm': [*result[3]]})
    d20 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'nedm': [*result2[0]]})
    d21 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'nedm': [*result2[1]]})
    d22 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'nedm': [*result2[2]]})
    d23 = pd.DataFrame({'theta': [*list_a],'tanbeta': [*list_b],'delta': [*list_c],
                           'nedm': [*result2[3]]})
    d1.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/nxy2'+str(exe.k) +'.csv',\
                   index = None, header=True)
    d2.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/nxy3'+str(exe.k) +'.csv',\
                   index = None, header=True)
    d10.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/n10'+str(exe.k) +'.csv',\
                   index = None, header=True)
    d11.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/n11'+str(exe.k) +'.csv',\
                   index = None, header=True)
    d12.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/n12'+str(exe.k) +'.csv',\
                   index = None, header=True)
    d13.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/n13'+str(exe.k) +'.csv',\
                   index = None, header=True)
    d20.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/n20'+str(exe.k) +'.csv',\
                   index = None, header=True)
    d21.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/n21'+str(exe.k) +'.csv',\
                   index = None, header=True)
    d22.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/n22'+str(exe.k) +'.csv',\
                   index = None, header=True)
    d23.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/n23'+str(exe.k) +'.csv',\
                   index = None, header=True)
    return
def mp1_imxy_mp2_edm():#3D scan for [MHp1,imxy,MHp2]against NEDM
    mp1 = np.arange(0,550,50)
    mp2 = np.arange(0,550,50)
    imxy = np.arange(0,42,2)
    nedmresult = []
    lista = []
    listb = []
    listc = []
    for i in mp1:
        for j in mp2:
            for m in imxy:
                lista.append(i)
                listb.append(j)
                listc.append(m)
                nedmresult.append(abs(dn(i,j,complex(0,m),\
                                complex(0,- m)) / (5.06e13)  / 1e-26 ) )
    dlist = pd.DataFrame({'$M_{H^{\pm}_1}}$': [*lista],'$M_{H^{\pm}_2}}$': [*listb],'delta': [*listc],
                           'nedm': [*nedmresult]})
    dlist.to_csv \
(r'/Users/muyuansong0419/Desktop/Allen/PHDmeeting/3HDM codes/mp1mp2imxynedm.csv',\
                   index = None, header=True)
    return
def plt_A_B_nedm(i,j):#[A,B] plane with MHP1 = i, MHp2 = j Nedm
    resultn = []
#nedm result
    for n in np.arange(0,len(ABarray4()) ):
        nedm3hdm = abs(dn(i,j, exe.complexyfunction(*ABarray4()[n]),\
                               exe.complexyfunction3(*ABarray4()[n]) ) / (5.06e13)  )
        resultn.append( nedm3hdm )
#########
    ned = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultn).flatten()  ,len(np.array(resultn).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,\
          levels = np.array([0.0,1.8e-26])  )
#           levels = np.arange(3.0,4.2,0.2),\
#          colors = ['black','royalblue','purple','darkgreen',\
#                    'brown','red','gray','orange','pink'])
    plt.colorbar(ned)
    plt.title('Neutron EDM with '\
                    + str("%02d" % i) +', ' + str("%02d"% j) )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    return
##############################################################################
def rexy_imxy_bsg_nedm(mm1,mm2,i22,i33):
#imn2min,max 0.0 42.9827952239115
#imn2min,max 0.0 42.9827952239115
#ren2min,max 0.0008955247801336608 43.58968705776075
#ren3min,max 0.0013616841467565938 43.32971016574854
    m1,m2 = mm1,mm2
    i2,i3 = i22,i33
    re2 = np.arange(- 43.6, 43.3)
    re3 = np.arange(- 43.3, 43.6)
    im2 = np.arange(-43,44)
    im3 = np.arange(-43,44)
    j2 = re2 + 1j * im2
    j3 = re3 + 1j * im3 
#    print(j2.real,j3)
    resultb = []
    resultn = []
    for j in j2.imag:
        for i in j2.real:
#B>Xs+gamma SECTION    
            y3hdm= bsg.BR_B_Xs_gamma(mb,mw,m1,m2,\
                        i2, complex(i,j) ,\
                        i3, complex(10,-j) ) 
            resultb.append(y3hdm / (1e-4) )
#Nedm SECTION    
            nedm3hdm = abs(dn(m1,m2, complex(i,j),\
                               complex(i,-j) ) / (5.06e13) / 1e-26  )
            resultn.append( nedm3hdm )     
#    print(resultn)
    #########
    ned = plt.contourf(j2.real, j2.imag, \
           np.resize(np.array(resultn).flatten()  ,len(np.array(resultn).flatten() ) ).\
          reshape(len(j2.imag),len(j2.real)) ,\
          levels = np.array([0.0,1.8]),colors = ['red']  )
#########
    bsgamm = plt.contourf(j2.real, j2.imag, \
           np.resize(np.array(resultb).flatten()  ,len(np.array(resultb).flatten() ) ).\
          reshape(len(j2.imag),len(j2.real)) ,\
          levels = np.array([2.99,3.55]),colors = ['green'] )
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) and NEDM in '\
                    + str("%02d" % m1) +', ' + str("%02d"% m2) )
    plt.xlabel('Re$(X_2Y_2^*)$')
    plt.ylabel('Im$(X_2Y_2^*)$')
#    plt.axis([-40,-15, -10,15])
#    plt.savefig('rexy2imxy2'+ str("%02d" % m1) + str("%02d"% m2) +'.png')
    plt.show()
    plt.close()
    return
###########################
def monte_carlo_rexy_imxy(mm1,mm2,i22,i33,samples):
    result = 0
    m1 ,m2 = mm1,mm2
    i2,i3 = i22,i33
    xbsg_list = []
    ybsg_list = []
    xnedm_list = []
    ynedm_list = []
    for _ in range(samples):
        x = random.uniform(-40,40)
        y = random.uniform(-40,40)
        if float(bsg.BR_B_Xs_gamma(mb,mw,m1,m2,\
                        i2, complex(x,y) ,\
                        i3, complex(0.1,- y) ) / (1e-4)) >= 2.99 and float(bsg.BR_B_Xs_gamma(mb,mw,m1,m2,\
                        i2, complex(x,y) ,\
                        i3, complex(0.1,- y) ) / (1e-4)) <= 3.55 :
            xbsg_list.append(x)
            ybsg_list.append(y)
        if float(abs(dn(m1,m2, complex(i,j),\
           complex(i,-j) ) / (5.06e13) / 1e-26  )) >=0 and float(abs(dn(m1,m2, complex(i,j),\
                               complex(i,-j) ) / (5.06e13) / 1e-26  )) <=1.1:
            xnedm_list.append(x)
            ynedm_list.append(y)
            result += 1
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.scatter(xbsg_list, ybsg_list, s=20, c='b')
    ax.scatter(xnedm_list, ynedm_list, s=20, c='r')
    plt.show()
    plt.close()
    return
        
##############################################################################
##############################################################################
##########################Neutron EDM PLOT SECTION
#plt_A_B_nedm(100,105)#[A,B] plane with MHP1 = i, MHp2 = j Nedm
#NEDMfigure4_plot()#Neutron EDM plot figure 4
#nedm3hdm_plot()#Neutron EDM in 3hdm-plot with exercise file
#nedm3hdm_xy1_xy2()# NEDM figure 4 3hdm
#nedm_3dparameter()#  NEDM to plot 3D based on 4 mixing parameters
#nedm_mhp1_massdiff()#[mhch1,mass_differ] in 3hdm nedm result
#nedm_mhp1_mhp2()#[mhch1,mhch2] in 3hdm nedm result
#mp1_imxy_mp2_edm()#3D scan for [MHp1,imxy,MHp2]against NEDM
################################################################################
################################################################################
for i in np.array([85,120,200,500]): 
    rexy_imxy_bsg_nedm(85,i,1,0.5)
    plt_A_B_bsgnedm(85,i)
#    plt_A_b_eedm(85,i)
#monte_carlo_rexy_imxy(85,200,10,1,100000)