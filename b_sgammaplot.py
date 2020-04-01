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
from neutronedm import dC_btH, dn_CEDM,CW,dn
import pandas as pd
#from exe import xyfun,xyfun3,yfun,yfun3,U,X2,X3,Y2,Y3,Z2,Z3,\
#        complexyfunction,complexyfunction3,read1,read2,A,B,readlist
########################
axs = np.arange(1, 61, 1)
#print('A,B',exe.A,exe.B,exe.read1,exe.read2,len(exe.B),type(exe.read1))
#Plot Functions
def Plot_3():#Figure 3 DOI: 10.1142/S0217751X17501457
    x_axis = np.array([ i for i in np.arange(100,1020,20)] )
    print(x_axis,type(x_axis),x_axis.shape)
    y11_axis = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                                      [1.0],[-1.0],[0],[0]) )
    y12_axis = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                                      [1.0/2.0],[-1.0/4.0],[0],[0]) )
    y130_axis = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                                      [1.0/30.0],[-1.0/900],[0],[0]) )
    y21_axis = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                                      [1.0],[1.0],[0],[0]) )
    y22_axis = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                                      [1.0/2.0],[1.0],[0],[0]) )
    y230_axis = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                                      [1.0/30.0],[1.0],[0],[0]) )
    plt.axis([100.0, 1000.0, 1.0, 7.0])
    plt.plot(x_axis,y11_axis / (1e-4))
    plt.plot(x_axis,y12_axis / (1e-4))
    plt.plot(x_axis,y130_axis / (1e-4))
    plt.plot(x_axis,y21_axis / (1e-4))
    plt.plot(x_axis,y22_axis / (1e-4) )
    plt.plot(x_axis,y230_axis / (1e-4))
    print('type I tanbeta = 1',y11_axis / (1e-4))
    print('type I tanbeta = 30',y130_axis / (1e-4))
    plt.xlabel('$M_{H^{\pm}}$')
    plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.title('Figure 3: BR($\\bar{B} \\to X_{s} \gamma$) VS. $M_{H^{\pm}_{1}}$ ' )
    plt.legend(('Type I tan$\\beta =$ 1', 'Type I tan$\\beta =$ 2', 'Type I tan$\\beta =$ 30',\
            'Type II tan$\\beta =$ 1', 'Type II tan$\\beta =$ 2', 'Type II tan$\\beta =$ 30'),
           loc='upper right', shadow=True,prop={'size': 7.8})
    plt.show()
    plt.close
###################################
def plotfigure1():
    mhch = 100
    
    x_axis = [ i for i in np.arange(5,91,1)[::-1]] 
    c0eff7sm = []
    c0eff7yy = []
    c0eff7xy = []
    ceff7sm = []
    ceff7yy = []
    ceff7xy = []
    for i in x_axis:
        c0eff7sm.append(bsg.onlyfigure1(i,80,mhch)[0])
        c0eff7yy.append(bsg.onlyfigure1(i,80,mhch)[1])
        c0eff7xy.append(bsg.onlyfigure1(i,80,mhch)[2])
        ceff7sm.append(bsg.onlyfigure1(i,80,mhch)[3])
        ceff7yy.append(bsg.onlyfigure1(i,80,mhch)[4])
        ceff7xy.append(bsg.onlyfigure1(i,80,mhch)[5])
#        print(i,bsg.onlyfigure1(i,mhch,500)[3])
    print('c0eff7xy',c0eff7xy)
    print('ceff7xy',ceff7xy)
    print('ceff7sm',ceff7sm)
    print('c0eff7sm',c0eff7sm)
    print('c0eff7yy',c0eff7yy)
    print('ceff7yy',ceff7yy)
    plt.xlim(80, 0)
    plt.ylim(-0.4, 0)
    plt.plot(x_axis,c0eff7sm)
    plt.plot(x_axis,ceff7sm)
    plt.grid(axis='x', linestyle='-', color='0.4') # show x-axis grid line
    plt.grid(axis='y', linestyle='-', color='0.4') # show x-axis grid line
    plt.legend(('$C^{0,eff}_{7,SM}$', '$C^{eff}_{7,SM}$'),
           loc='upper left', shadow=True,prop={'size': 10})
    plt.show()
    plt.close()
    
    plt.xlim(80, 0)
    plt.ylim(-0.4, 0)
    plt.plot(x_axis,c0eff7yy)
    plt.plot(x_axis,c0eff7xy)
    plt.plot(x_axis,ceff7yy)
    plt.plot(x_axis,ceff7xy)
    plt.grid(axis='x', linestyle='-', color='0.4') # show x-axis grid line
    plt.grid(axis='y', linestyle='-', color='0.4') # show x-axis grid line
    plt.legend(('$C^{0,eff}_{7,YY}$', '$C^{0,eff}_{7,XY}$','$C^{eff}_{7,YY}$',\
                '$C^{eff}_{7,XY}$'),
           loc='lower right', shadow=True,prop={'size': 9})
    plt.show()
    plt.close()
    return 
###################################
def Plot4_8_9():#Boltzmati's Paper
    
    x_axis = np.arange(-10.0,10.25,0.25)# figure 8
    XYimx_axis = [complex(-2,i)*1.0 for i in x_axis]
    rangephi = np.arange(0,np.pi,0.01) # figure 9
#    print('rangephi', rangephi,len(rangephi))
    XYexpim_axis = [  complex(np.cos(j),np.sin(j)) for j in rangephi] 
#    print('REALX,IMX:',[np.complex(-2,i)*1.0 for i in xim_axis],XYimx_axis)
#    print('X = 2exp(i phi)',XYexpim_axis,len(XYexpim_axis))
    mhch = 100
    y48_axis = bsg.BR_B_Xs_gamma(4.8,mhch,mhch,mhch + bsg.mass_differ ,\
                        [1.0],x_axis ,[0.0],[0.0])
    y24_axis = bsg.BR_B_Xs_gamma(2.4,mhch,mhch,mhch + bsg.mass_differ ,\
                        [1.0],x_axis ,[0.0],[0.0])
    y96_axis = bsg.BR_B_Xs_gamma(9.6,mhch,mhch,mhch + bsg.mass_differ ,\
                        [1.0],x_axis ,[0.0],[0.0])
    y48imx_axis = bsg.BR_B_Xs_gamma(4.8,mhch,mhch,mhch + bsg.mass_differ ,\
                        [1.0],XYimx_axis,[0.0],[0.0])
    y24imx_axis = bsg.BR_B_Xs_gamma(2.4,mhch,mhch,mhch + bsg.mass_differ,\
                        [1.0],XYimx_axis,[0.0],[0.0])
    y96imx_axis = bsg.BR_B_Xs_gamma(9.6,mhch,mhch,mhch + bsg.mass_differ,\
                        [1.0],XYimx_axis,[0.0],[0.0])
    y48phi_axis = bsg.BR_B_Xs_gamma(4.8,mhch,mhch,mhch + bsg.mass_differ ,\
                        [0.5],XYexpim_axis,[0.0],[0.0])
    y24phi_axis = bsg.BR_B_Xs_gamma(2.4,mhch,mhch,mhch + bsg.mass_differ,\
                        [0.5],XYexpim_axis,[0.0],[0.0])
    y96phi_axis = bsg.BR_B_Xs_gamma(9.6,mhch,mhch,mhch + bsg.mass_differ,\
                        [0.5],XYexpim_axis,[0.0],[0.0])
#    print('----',XYimx_axis)
#    print(x_axis)
#    print('48',y48_axis * (1e4))
#    print('24',y24_axis * (1e4))
#    print('96',y96_axis * (1e4))
#    print('48im',y48imx_axis * (1e4) )
#    print('24im',y24imx_axis * (1e4))
#    print('96im',y96imx_axis * (1e4))
    plt.xlim(-10, 2)
    plt.ylim(-5, 10)
    plt.plot(x_axis,y48_axis * (1e4))
    plt.plot(x_axis,y24_axis * (1e4))
    plt.plot(x_axis,y96_axis * (1e4))
    plt.grid(axis='x', linestyle='-', color='0.4') # show x-axis grid line
    plt.grid(axis='y', linestyle='-', color='0.4') # show x-axis grid line
    plt.xlabel('X')
    plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.title('Figure4' )
    plt.legend(('$\mu = 4.8$ GeV', '$\mu = 2.4$ GeV', '$\mu = 9.6$ GeV '),
           loc='lower left', shadow=True,prop={'size': 8})
    plt.show()

    plt.xlim(-7, 7)
    plt.ylim(-2, 7)
    plt.plot(x_axis,y48imx_axis * (1e4))
    plt.plot(x_axis,y24imx_axis * (1e4))
    plt.plot(x_axis,y96imx_axis * (1e4))
#plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.grid(axis='y', linestyle='-', color='0.75') # show x-axis grid line
    plt.xlabel('Im(X)')
    plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.title('Figure8' )
    plt.legend(('$\mu = 4.8$ GeV', '$\mu = 2.4$ GeV', '$\mu = 9.6$ GeV '),
           loc='lower right', shadow=True,prop={'size': 8})
    plt.show()
    plt.close
    plt.plot(rangephi,y48phi_axis / (1e-4))
    plt.plot(rangephi,y24phi_axis / (1e-4))
    plt.plot(rangephi,y96phi_axis / (1e-4))
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.grid(axis='y', linestyle='-', color='0.75') # show x-axis grid line
    plt.xlabel('$\\phi$')
    plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.title('Figure9' )
    plt.legend(('$\mu = 4.8$ GeV', '$\mu = 2.4$ GeV', '$\mu = 9.6$ GeV '),
           loc='upper right', shadow=True,prop={'size': 8})
    plt.xlim(0, np.pi)
    plt.ylim(0, 8.0)
    plt.show()
    plt.close
def Plot_4() :#Andrew and Stefano Paper
    x_axis = np.array([ i for i in np.arange(100,1020,20)] )
#    print('II 2HDM tanbeta = 2',BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,\
#                                              [1.0/2.0],[1.0],[0],[0])\
#            )
    theta_c = bsg.PI/4.0
    tan_beta = 2.0
    Y1_array = - 1.0/tan_beta * np.cos(theta_c)
    X1_array =  1.0/ tan_beta  * np.cos(theta_c) 
    X2_array = - 1.0/tan_beta * np.sin(theta_c)
    Y2_array = 1.0/tan_beta * np.sin(theta_c)
    y12_2hdm = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                                      [1.0/2.0],[-1.0/4.0],[0],[0]) )
    y12_3hdm = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                     [Y1_array],[np.array(X1_array) * np.conjugate(Y1_array)],\
                     [Y2_array],[np.array(X2_array) * np.conjugate(Y2_array)]) )
    plt.plot(x_axis,y12_2hdm / (1e-4))
    plt.plot(x_axis,y12_3hdm / (1e-4))
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.title('Figure 4:$\\theta = \\pi$/4,$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' \
                   + str(bsg.mass_differ)+ ' GeV' )
    plt.legend(('Type I tan$\\beta =$ 2 2HDM', 'Type I tan$\\beta =$ 2 3HDM '))
#    for n in np.arange(0,len(axs)):
        
#        y22_axis3hdm= BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,\
#                        Y2(*array4()[n]),complexyfunction(*array4()[n]),\
#                        Y3(*array4()[n]),complexyfunction3(*array4()[n])) 
#        plt.plot(x_axis,y22_axis3hdm / (1e-4))
    plt.axis([100.0, 1000.0, 1.0, 7.0])
    plt.show()
    plt.close()
    return
def Plot_5() :#Andrew and Stefano Paper
    x_axis = np.array([ i for i in np.arange(100,1020,20)] )
#    print('II 2HDM tanbeta = 2',BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,\
#                                              [1.0/2.0],[1.0],[0],[0])\
#            )
#    cpphase = np.exp(complex(0,PI))
    y22_2hdm = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                                      [1.0/2.0],[1.0],[0],[0]) )
    plt.plot(x_axis,y22_2hdm / (1e-4))
    theta_c = - bsg.PI/4.0
    tan_beta = 2
    cos_beta = 1 / np.sqrt(1 + tan_beta**2)
#    tan_gamma = n
    for n in np.array([2,4,7,10,20]):
        
        X1_array =  - tan_beta  * np.cos(theta_c) - n / (cos_beta) * np.sin(theta_c) #\
                       # * cpphase
        Y1_array = - 1.0/tan_beta * np.cos(theta_c)
        X2_array =  tan_beta  * np.sin(theta_c) - n / (cos_beta) * np.cos(theta_c) #\
                       # * cpphase
        Y2_array = 1.0/tan_beta * np.sin(theta_c)
        
        y12_3hdm = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                     [Y1_array],[np.array(X1_array) * np.conjugate(Y1_array)],\
                     [Y2_array],[np.array(X2_array) * np.conjugate(Y2_array)]) )
#        print(y12_3hdm/ (1e-4))
        plt.plot(x_axis,y12_3hdm / (1e-4))
        plt.xlabel('$M_{H^{\pm}_{1}}$')
        plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
        plt.title('Figure 5:$\\theta = -\\pi$/4, $M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' \
                   + str(bsg.mass_differ)+ ' GeV' )
    plt.legend(('Type II tan$\\beta =$ 2 2HDM', 'Type II tan$\\gamma =$ 2 3HDM ',\
                    'Type II tan$\\gamma =$ 4 3HDM ','Type II tan$\\gamma =$ 7 3HDM ',\
                    'Type II tan$\\gamma =$ 10 3HDM ','Type II tan$\\gamma =$ 20 3HDM '),\
               loc='upper right', shadow=True,prop={'size': 7.5} )    
#    for n in np.arange(0,len(axs)):
        
#        y22_axis3hdm= BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,\
#                        Y2(*array4()[n]),complexyfunction(*array4()[n]),\
#                        Y3(*array4()[n]),complexyfunction3(*array4()[n])) 
#        plt.plot(x_axis,y22_axis3hdm / (1e-4))
    plt.axis([100.0, 1000.0, 1.0, 6.0])
    plt.show()
    plt.close()
    return
######################################################################
def plot_under_Heatherbasis(i,j,k,l):
    m1_axis = np.array([ i for i in np.arange(50,550,50)] )
    m2_axis = np.array([ i for i in np.arange(50,1050,50)] )
    m2 = m2_axis[0]
    m1 = m1_axis[0]
    print('i,j,k,l',i,j,k,l)
#    xx, yy = np.meshgrid(m1_axis, m2_axis)
    empty = []
    for m2 in m2_axis:
        for m1 in m1_axis:
            threehdm = bsg.BR_B_Xs_gamma(mb,mw,m1,m2,\
                        [exe.Y2(i,j,k,l)],[- exe.complexyfunction(i,j,k,l)],\
                        [exe.Y3(i,j,k,l)],[- exe.complexyfunction3(i,j,k,l)])
            empty.append(threehdm)
    result = plt.contourf(m1_axis, m2_axis, \
           np.resize(np.array(empty) / (1e-4),len(np.array(empty) / (1e-4))).\
           reshape(len(m2_axis),len(m1_axis)), \
           colors = ['black','royalblue','purple','darkgreen','brown','red','gray','orange'],\
           levels = np.array([3.11,3.5,3.87])      
           )
    plt.colorbar(result)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([50,200, 120.0, 180.0])
    plt.show()
    plt.close()
    return
######################################################################
def plot_under_deltascan(i,j,k,l):
    m1_axis = np.array([ i for i in np.arange(50,550,50)] )
    m2_axis = np.array([ i for i in np.arange(50,1050,50)] )
    m2 = m2_axis[0]
    m1 = m1_axis[0]
#    print('i,j,k,l',i,j,k,l)
#    xx, yy = np.meshgrid(m1_axis, m2_axis)
    print('i,j,k,l',i,j,k,l)
    empty = []
    for m2 in m2_axis:
        for m1 in m1_axis:
            threehdm = bsg.BR_B_Xs_gamma(mb,mw,m1,m2,\
                        [exe.Y2(i,j,k,l)],[- exe.X2(i,j,k,l) * np.conjugate(exe.Y2(i,j,k,l) )],\
                        [exe.Y3(i,j,k,l)],[- exe.X3(i,j,k,l) * np.conjugate(exe.Y3(i,j,k,l) )])
            empty.append(threehdm)
    result = plt.contourf(m1_axis, m2_axis, \
           np.resize(np.array(empty) / (1e-4),len(np.array(empty) / (1e-4))).\
           reshape(len(m2_axis),len(m1_axis)), \
           colors = ['black','royalblue','purple','darkgreen','brown','red','gray','orange'],\
           levels = np.array([3.00,3.55])      
           )
    plt.colorbar(result)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([50,200, 50.0, 1000.0])
    plt.axis([0,500, 0.0, 1000.0])
    plt.show()
    plt.close()
def plot_Hp1_Hp2():# [mH+1,mH+2] for fixed B_bar > X_s + gamma
#    fixedarray = (- bsg.PI/2.1,10,60,0.0) #mixing matrix parameters
    tangamma = 3.0#tangamma
    theta = - bsg.PI/ 4
    tanbeta = 2
    X1_array =  - tanbeta  * np.cos(theta) - tangamma / (1.0 / np.sqrt(1.0 + tanbeta**2)) * np.sin(theta)
    Y1_array = - 1.0/tanbeta * np.cos(theta)
    X2_array =  tanbeta  * np.sin(theta) - tangamma / (1.0 / np.sqrt(1.0 + tanbeta**2))  * np.cos(theta)
    Y2_array = 1.0/tanbeta * np.sin(theta)
#    print('X1',X1_array,Y1_array,X2_array,Y2_array)
    m1_axis = np.array([ i for i in np.arange(50,550,20)] )
    m2_axis = np.array([ i for i in np.arange(50,1050,20)] )
    empty =[]
    m2 = m2_axis[0]
    m1 = m1_axis[0]
#    xx, yy = np.meshgrid(m1_axis, m2_axis)
    for m2 in m2_axis:
        for m1 in m1_axis:
            threehdm = bsg.BR_B_Xs_gamma(mb,mw,m1,m2,\
                        [Y1_array],[ X1_array * np.conjugate(Y1_array)],\
                        [Y2_array],[ X2_array * np.conjugate(Y2_array)]) 
            empty.append(threehdm)
    result = plt.contourf(m1_axis, m2_axis, \
           np.resize(np.array(empty) / (1e-4),len(np.array(empty) / (1e-4))).\
           reshape(len(m2_axis),len(m1_axis)), \
           colors = ['black','royalblue','purple','darkgreen','brown','red','gray','orange'],\
           levels = np.array([2.99,3.87])         
         )
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.colorbar(result)
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([90,200, 75.0, 300.0])
    plt.show()
    plt.close()
    return
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
def plt_A_B_xy():# A_B XY*
    resultxy2real = []
    resultxy3real = []
    resultxy2imag = []
    resultxy3imag = []
    for n in np.arange(0,len(ABarray4()) ):
        resultxy2real.append( (- exe.complexyfunction(*ABarray4()[n] )).real  )
        resultxy3real.append( (- exe.complexyfunction3(*ABarray4()[n] )).real )
        resultxy2imag.append( (- exe.complexyfunction(*ABarray4()[n] )).imag  )
        resultxy3imag.append( (- exe.complexyfunction3(*ABarray4()[n] )).imag )
    #########
    xy2real = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultxy2real).flatten()  ,len(np.array(resultxy2real).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,levels = np.arange(0,7,1)   )
    plt.colorbar(xy2real)
    plt.title('Re$X_{2} Y_{2}^{*} $')
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
#    plt.axis([0,30,0,30])
    plt.savefig( 'xy2real.png')
    plt.show()
    plt.close()
    
    xy3real = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultxy3real).flatten()  ,len(np.array(resultxy3real).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,levels = np.arange(-4,2,0.8)    )
    plt.colorbar(xy3real)
    plt.title('Re$X_{3} Y_{3}^{*} $')
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
#    plt.axis([0,30,0,30])
    plt.savefig('xy3real.png')
    plt.show()
    plt.close()
    #########
    xy2imag = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultxy2imag).flatten()  ,len(np.array(resultxy2imag).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) )#,levels = np.array([0,0.1,0.2,0.3,0.4,0.5])     )
    plt.colorbar(xy2imag)
    plt.title('Im$X_{2} Y_{2}^{*} $')
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
#    plt.axis([0,30,0,30])
    plt.savefig( 'xy2imag.png')
    plt.show()
    plt.close()
    
    xy3imag = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultxy3imag).flatten()  ,len(np.array(resultxy3imag).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) )# ,levels = np.array([-0.5,-0.4,-0.3,-0.2,-0.1,0])    )
    plt.colorbar(xy3imag)
    plt.title('Im$X_{3} Y_{3}^{*} $')
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
#    plt.axis([0,30,0,30])
    plt.savefig('xy3imag.png')
    plt.show()
    plt.close()
##################################################################################
def plt_A_B_bsg(i,j):# Bsgamma-result in {A,B} plane
    mass_axis1,mass_axis2 = i,j
    print(ABarray4()[0],mass_axis1,len(ABarray4()))
    print(ABarray4()[1],mass_axis2,len(ABarray4()))
    resultb = []
    resultn = []
#B>Xs+gamma SECTION    
    for n in np.arange(0,len(ABarray4()) ):
        y3hdm= bsg.BR_B_Xs_gamma(mb,mw,mass_axis1,mass_axis2,\
                        exe.Y2(*ABarray4()[n] ), exe.complexyfunction(*ABarray4()[n] ),\
                        exe.Y3(*ABarray4()[n] ), exe.complexyfunction3(*ABarray4()[n] )) 
        resultb.append(y3hdm / (1e-4) )
#Nedm SECTION    
        nedm3hdm = abs(dn(mass_axis1,mass_axis2, exe.complexyfunction(*ABarray4()[n]),\
                               exe.complexyfunction3(*ABarray4()[n]) ) / (5.06e13)  )
        resultn.append( nedm3hdm )
#########
    ned = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultn).flatten()  ,len(np.array(resultn).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,\
          levels = np.array([0.0,3.3e-26]),colors = ['red']  )
#########
    bsgamm = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultb).flatten()  ,len(np.array(resultb).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,\
          levels = np.array([3.00,3.55]),colors = ['green'] )
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) and NEDM in '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2) )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    plt.savefig( str("%02d" % mass_axis1) + str("%02d"% mass_axis2) +'bsg.png')
    plt.show()
    plt.close()
####################################################################################
def plt_A_B_cp(i,j):#CP-asymmetry in {A,B} plane
    result_cp = []
    mass_axis1,mass_axis2 = i,j
    for n in np.arange(0,len(ABarray4()) ):
        cpasymetry = bsg.newa_cp(mb,mw,mass_axis1,mass_axis2,\
                        exe.Y2(*ABarray4()[n] ), exe.complexyfunction(*ABarray4()[n] ),\
#                        [0.0],[0.0])
                        exe.Y3(*ABarray4()[n] ), exe.complexyfunction3(*ABarray4()[n] )) 
        result_cp.append(cpasymetry)
    cpresult = plt.contourf(exe.A, exe.B, \
                            np.resize(np.array(result_cp).flatten()  ,\
                            len(np.array(result_cp).flatten() ) ).\
                            reshape(len(exe.B),len(exe.A)), \
         )# levels = np.array([-12,-10,-8,-6,-4,-2,0,2,4]) )
    plt.colorbar(cpresult)
    plt.title('CP-asymmetry with charged Higgs: '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2)+' GeV.' )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    plt.savefig( 'cp'+ str("%02d" % mass_axis1) + str("%02d"% mass_axis2) +'.png')
    plt.show()
    plt.close()
def plt_A_B_cpd(i,j):#CPD-asymmetry in {A,B} plane
    result_cpd = []
    mass_axis1,mass_axis2 = i,j
    for n in np.arange(0,len(ABarray4()) ):
        cpdasymetry = bsg.newa_cpd(mb,mw,mass_axis1,mass_axis2,\
                        exe.Y2(*ABarray4()[n] ), exe.complexyfunction(*ABarray4()[n] ),\
#                        [0.0],[0.0])
                        exe.Y3(*ABarray4()[n] ), exe.complexyfunction3(*ABarray4()[n] )) 
        result_cpd.append(cpdasymetry)
    cpresult = plt.contourf(exe.A, exe.B, \
                            np.resize(np.array(result_cpd).flatten()  ,\
                            len(np.array(result_cpd).flatten() ) ).\
                            reshape(len(exe.B),len(exe.A)), \
          )#levels = np.arange(-20,-8,2) )
    plt.colorbar(cpresult)
    plt.title('CPD-asymmetry with charged Higgs: '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2)+' GeV.' )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    plt.savefig( 'cpd'+ str("%02d" % mass_axis1) + str("%02d"% mass_axis2) +'.png')
    plt.show()
    plt.close()
def cpijkl(m1,m2,i,j,k,l):#CP-asymmetry as theta, tanbeta, tangamma, delta
    if m2 == 0:
        cpresult = bsg.newa_cp(mb,mw,m1,500,\
                        exe.Y2(i,j,k,l), exe.complexyfunction(i,j,k,l ),\
                        [0.0],[0.0])
    else:
        cpresult = bsg.newa_cp(mb,mw,m1,m2,\
                        exe.Y2(i,j,k,l), exe.complexyfunction(i,j,k,l ),\
#                        [0.0],[0.0])
                        exe.Y3(i,j,k,l), exe.complexyfunction3(i,j,k,l )) 
    print('cpresult',cpresult)
    return
######################################################################################
    return
def numerical():
    mass_axis = (80.0,250.0)
    result = []
    for n in np.arange(0,len(ABarray4()) ):
        y3hdm= bsg.BR_B_Xs_gamma(mb,mw,mass_axis[0],mass_axis[1],\
                        exe.Y2(*ABarray4()[n] ),- exe.complexyfunction(*ABarray4()[n] ),\
                        exe.Y3(*ABarray4()[n] ),- exe.complexyfunction3(*ABarray4()[n] )) 
#        print(y3hdm / (1e-4),n)
        result.append(y3hdm / (1e-4) )
    return np.concatenate(result).ravel()
def scanscanplot():#scan each parameter plots
    for m in exe.tbe:
        plot_under_deltascan(exe.i,m,exe.k, exe.l)
        plt.savefig(str(m) + '.png')
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
    m_axis = np.array([ i for i in np.arange(100,550,50)] )
    m1_axis = np.array([ i for i in np.arange(0,1020,20)] )
    xy_axis = np.array([ i for i in np.arange(0,1.1,0.1)] )
    xy1_axis = np.array([ i for i in np.arange(-5,5.1,0.1)] )
    longarray = np.array([i for i in np.arange(-1,1.2,0.2)])
    arrayxy1xy2 =  np.array([4,12,24,40])# Maximum of Im(XY^*) can get
    empty = []
    for j in xy1_axis:
        for i in m1_axis:
            # 1GeV =  5.06e13 cm^{-1}
            nedm = abs(dn(i,300,[complex(0,j)],[complex(0,0)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
            
            empty.append(nedm )
#            print('nedm',i,j,nedm )
    result = plt.contourf(m1_axis, xy1_axis, \
           np.resize(np.array(empty),len(np.array(empty) )).\
           reshape(len(xy1_axis),len(m1_axis)), \
           levels = np.array([0.0,1.5e-26,3.3e-26])      
           )
    plt.colorbar(result)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('IM($X_2Y_2^*$)')
    plt.title('NEDM in 2HDM')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([80,500, 0, 1.0])
    plt.show()
    plt.close()
#################
    empty3hdm = []
    for j in xy1_axis:
        for i in m1_axis:
            # 1GeV =  5.06e13 cm^{-1}
                nedm = abs(dn(i,masss2,[complex(0,j)],[ complex(0,- j)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
            
                empty3hdm.append(nedm )
#            print('nedm',i,j,nedm )
    result8 = plt.contourf(m1_axis, xy1_axis, \
                          np.resize(np.array(empty3hdm),len(np.array(empty3hdm) )).\
                          reshape(len(xy1_axis),len(m1_axis)), \
                          levels = np.array([0.0,1.1e-26,2.2e-26,3.3e-26])      
                          )
    plt.colorbar(result8)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('IM($X_2Y_2^*$)')
    plt.title('NEDM in 3HDM with $M_{H^{\pm}_2} = $' + str(masss2))
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([50,500, 0, 1.0])
    plt.show()
#        plt.savefig('m1' + str(masss2) + 'imxy2'+ str("%.2f" % m) + '.png' )
    plt.close()
##################
    
    for m in arrayxy1xy2:
        threehdm = []
        for j in massdiffer_list2:
            for i in m1_axis:
             # 1GeV =  5.06e13 cm^{-1}
                nedm00 = abs(dn(i,i + j,[complex(0,m)],[ complex(0,- m)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
                threehdm.append(nedm00 )
        result9 = plt.contourf(m1_axis, massdiffer_list2, \
                          np.resize(np.array(threehdm),len(np.array(threehdm) )).\
                          reshape(len(massdiffer_list2),len(m1_axis)), \
                          levels = np.array([0.0,1.1e-26,2.2e-26,3.3e-26])      
                          )
        plt.colorbar(result9)
        plt.xlabel('$M_{H^{\pm}_{1}}$')
        plt.ylabel('Mass difference')
#        plt.title('NEDM with IM($X_2Y_2^*$) = $\sqrt{2}$, IM($X_3Y_3^*$) = - $\sqrt{2}$' )
        plt.title('NEDM with IM($X_2Y_2^*$) = '+ str('%.2g'% m) + \
              ', IM($X_3Y_3^*$) = - ' + str('%.2g'% m) )
        plt.axis([0,500,0, 80])
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
                                   levels = np.array([0.0,3.3e-26]),      
                                   colors='green',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 20')
    result40 = plt.contour(longarray, longarray, \
                                   np.resize(np.array(empty40),len(np.array(empty40) )).\
                                   reshape(len(longarray),len(longarray)), \
                                   levels = np.array([0.0,3.3e-26]),      
                                   colors='black',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 40')
    result60 = plt.contour(longarray, longarray, \
                                   np.resize(np.array(empty60),len(np.array(empty60) )).\
                                   reshape(len(longarray),len(longarray)), \
                                   levels = np.array([0.0,3.3e-26]),      
                                   colors='blue',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 60')
    result80 = plt.contour(longarray, longarray, \
                                   np.resize(np.array(empty80),len(np.array(empty80) )).\
                                   reshape(len(longarray),len(longarray)), \
                                   levels = np.array([0.0,3.3e-26]),      
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
           levels = np.array([0,3.3e-26])      
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
            nedm4 = abs(dn(i,j,[complex(0,2)],[complex(0,- 2)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
            
            empty4.append(nedm4 )
#            print('nedm',i,j,nedm )
    result4 = plt.contourf(m1_axis, m1_axis, \
           np.resize(np.array(empty4),len(np.array(empty4) )).\
           reshape(len(m1_axis),len(m1_axis)), \
           levels = np.array([0.0,1.1e-26,2.2e-26,3.3e-26])      
           )
    plt.colorbar(result4)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('NEDM with IM($X_2Y_2^*$) = 2, IM($X_3Y_3^*$) = - 2')
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
                                 levels = np.array([0.0,1.1,2.2,3.3,4.4])        
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
                                 levels = np.array([0.0,1.1,2.2,3.3])        
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
    masss = 100 # first charged Higgs
    masssdif = [20,40,60,80]#mass differ
    xx_ary = np.array([-0.1,0.0])
    yy_ary = np.array([0.1,0.0,-0.1,-0.2])
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
               levels = np.array([0.0,1.1,2.2,3.3]),      
               colors='green',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 20')
#    plt.colorbar(result20)
#    plt.text(4, yy_ary[0] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [0]) )
#    plt.text(0.5, yy_ary[1] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [1]) )
#    plt.text(0.5, yy_ary[2] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [2]) )
#    plt.text(0.5, yy_ary[3] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [3]) )
#    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss))
    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) \
              + ', $M_{H^{\pm}_{2}}$ = '+ str("%3d" % (masss +  masssdif [0])) )  
#   plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([-0.3,0.3,-0.3, 0.3])
    plt.xlabel('IM($X_2Y_2^*$)')
    plt.ylabel('IM($X_3Y_3^*$)')
    plt.show()
    plt.close()
#    plt.plot(xx_ary, np.array([yy_ary[0]] * len(xx_ary)),color='green')
#    plt.plot(xx_ary, np.array([yy_ary[1]] * len(xx_ary)),color='black')
#    plt.plot(xx_ary, np.array([yy_ary[2]] * len(xx_ary)),color='blue')
#    plt.plot(xx_ary, np.array([yy_ary[3]] * len(xx_ary)),color='red')
#    plt.text(xx_ary[1]+0.05, yy_ary[0] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [0]) )
#    plt.text(xx_ary[1]+0.05, yy_ary[1] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [1]) )
#    plt.text(xx_ary[1]+0.05, yy_ary[2] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [2]) )
#    plt.text(xx_ary[1]+0.05, yy_ary[3] - 0.05, r'$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ +' +  str(masssdif [3]),size = 10 )
#####################
    result40 = plt.contour(np.imag(xy1),np.imag(xy2), \
               np.resize(np.array(empty40).flatten(),len(np.array(empty40).flatten() )).\
               reshape(len(np.imag(xy2)),len(np.imag(xy1))), \
               levels = np.array([0.0,1.1,2.2,3.3]),      
               colors='black',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 40')
    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) \
              + ', $M_{H^{\pm}_{2}}$ = '+ str("%3d" % (masss +  masssdif [1])) )
    plt.xlabel('IM($X_2Y_2^*$)')
    plt.ylabel('IM($X_3Y_3^*$)')
    plt.show()
    plt.close()
#################
    result60 = plt.contour(np.imag(xy1),np.imag(xy2), \
               np.resize(np.array(empty60).flatten(),len(np.array(empty60).flatten() )).\
               reshape(len(np.imag(xy2)),len(np.imag(xy1))), \
               levels = np.array([0.0,1.1,2.2,3.3]),      
               colors='blue',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 60')
    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) \
              + ', $M_{H^{\pm}_{2}}$ = '+ str("%3d" % (masss +  masssdif [2])) )
    plt.xlabel('IM($X_2Y_2^*$)')
    plt.ylabel('IM($X_3Y_3^*$)')
    plt.show()
    plt.close()
    result80 = plt.contour(np.imag(xy1),np.imag(xy2), \
               np.resize(np.array(empty80).flatten(),len(np.array(empty80).flatten() )).\
               reshape(len(np.imag(xy2)),len(np.imag(xy1))), \
               levels = np.array([0.0,1.1,2.2,3.3]),      
               colors='red',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 80')
    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) \
              + ', $M_{H^{\pm}_{2}}$ = '+ str("%3d" % (masss +  masssdif [3])) )
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
               levels = np.array([0.0,1.1,2.2,3.3]))
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
               levels = np.array([0.0,1.1,2.2,3.3]))
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
                                complex(0,- m)) / (5.06e13) ) / 1e-26)
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
          levels = np.array([0.0,1.1e-26,2.2e-26,3.3e-26])  )
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
def rexy_imxy_bsg_nedm():
#imn2min,max 0.0 42.9827952239115
#imn2min,max 0.0 42.9827952239115
#ren2min,max 0.0008955247801336608 43.58968705776075
#ren3min,max 0.0013616841467565938 43.32971016574854
    m1,m2 = 85, 200
    i2,i3 = 5,5
    re2 = np.arange(- 43.6, 43.3)
    re3 = np.arange(- 43.3, 43.6)
    im2 = np.arange(-43,44)
    im3 = np.arange(-43,44)
    j2 = re2 + 1j * im2
    j3 = re3 + 1j * im3 
    resultb = []
    resultn = []
    for i in j2.real:
        for j in j2.imag:
#B>Xs+gamma SECTION    
            y3hdm= bsg.BR_B_Xs_gamma(mb,mw,m1,m2,\
                        i2, complex(i,j) ,\
                        i3, complex(i,-j) ) 
            resultb.append(y3hdm / (1e-4) )
#Nedm SECTION    
            nedm3hdm = abs(dn(m1,m2, complex(i,j),\
                               complex(i,-j) ) / (5.06e13)  )/ 1e-26
            resultn.append( nedm3hdm )       
    #########
    ned = plt.contourf(j2.real, j2.imag, \
           np.resize(np.array(resultn).flatten()  ,len(np.array(resultn).flatten() ) ).\
          reshape(len(j2.imag),len(j2.real)) ,\
          levels = np.array([0.0,3.3]),colors = ['red']  )
#########
    bsgamm = plt.contourf(j2.real, j2.imag, \
           np.resize(np.array(resultb).flatten()  ,len(np.array(resultb).flatten() ) ).\
          reshape(len(j2.imag),len(j2.real)) ,\
          levels = np.array([3.00,3.55]),colors = ['green'] )
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) and NEDM in '\
                    + str("%02d" % m1) +', ' + str("%02d"% m2) )
    plt.xlabel('Re$(X_2Y_2^*)$')
    plt.ylabel('Im$(X_2Y_2^*)$')
    plt.axis([-15,15, -20,20])
    plt.savefig('rexy2imxy2'+ str("%02d" % m1) + str("%02d"% m2) +'.png')
    plt.show()
    plt.close()
    
    resultb2 = []
    resultn2 = []
    for i in j2.real:
        for j in j3.imag:
#B>Xs+gamma SECTION    
            y3hdm= bsg.BR_B_Xs_gamma(mb,mw,m1,m2,\
                        i2, complex(i,-j) ,\
                        i3, complex(i,j) ) 
            resultb2.append(y3hdm / (1e-4) )
#Nedm SECTION    
            nedm3hdm = abs(dn(m1,m2, complex(i,-j),\
                               complex(i,j) ) / (5.06e13)  )/ 1e-26
            resultn2.append( nedm3hdm )       
    #########
    ned2 = plt.contourf(j2.real, j3.imag, \
           np.resize(np.array(resultn2).flatten()  ,len(np.array(resultn2).flatten() ) ).\
          reshape(len(j2.imag),len(j2.real)) ,\
          levels = np.array([0.0,3.3]),colors = ['red']  )
#########
    bsgamm2 = plt.contourf(j2.real, j3.imag, \
           np.resize(np.array(resultb2).flatten()  ,len(np.array(resultb2).flatten() ) ).\
          reshape(len(j2.imag),len(j2.real)) ,\
          levels = np.array([3.00,3.55]),colors = ['green'] )
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) and NEDM in '\
                    + str("%02d" % m1) +', ' + str("%02d"% m2) )
    plt.xlabel('Re$(X_2Y_2^*)$')
    plt.ylabel('Im$(X_3Y_3^*)$')
    plt.axis([-15, 15, -20,20])
    plt.savefig('rexy2imxy3'+ str("%02d" % m1) + str("%02d"% m2) +'.png')
    plt.show()
    plt.close()
    return
##############################################################################
##############################################################################
##################PLOT SECTION for B> Xs + Gamma
#numerical()
#Plot_3()
#Plot_4()
#Plot_5()
#Plot4_8_9()
#plotfigure1()
#plot_Hp1_Hp2()
#plot_under_Heatherbasis(exe.i,exe.j,exe.k,exe.l)
#cpijkl(85,165,- np.pi/4,5,1,0)#CP-asymmetry as theta, tanbeta, tangamma, delta
#############################################
#starttime = exe.time.time()
#scanscanplot()
#plot_under_deltascan(exe.i,exe.j,exe.k, exe.l)
#endtime = exe.time.time()
#print(endtime - starttime)
#plt_A_B_bsg(100,105)# Bsgamma-result in {A,B} plane
#plt_A_B_xy()#A_B_XY*
#for jjjj in np.arange(100,1100,200):
#        plt_A_B_bsg(85,jjjj )#A_B_bsgamma light light
#for iiii in np.array([100,200,400,600]):
#        plt_A_B_bsg(85,iiii )#A_B_bsgamma light heavy
#        plt_A_B_bsg(200,iiii)#A_B_bsgamma heavy heavy
#        plt_A_B_cp(85,iiii) #A_B_CP-asymmetry
#        plt_A_B_cpd(85,iiii)##A_B_CPD-asymmetry
################################################################################
################################################################################    
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
rexy_imxy_bsg_nedm()