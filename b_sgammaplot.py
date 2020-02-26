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
import random
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
           levels = np.array([3.11,3.5,3.87])      
           )
    plt.colorbar(result)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([50,200, 50.0, 1000.0])
    plt.axis([50,200, 100.0, 170.0])
def plot_Hp1_Hp2():# [mH+1,mH+2] for fixed B_bar > X_s + gamma
#    fixedarray = (- bsg.PI/2.1,10,60,0.0) #mixing matrix parameters
    tangamma = 5.0#tangamma
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
           levels = np.array([2.99,3.5,3.85])         
         )
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.colorbar(result)
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([90,200, 75.0, 300.0])
#    plt.show()
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
          reshape(len(exe.B),len(exe.A))     )
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
          reshape(len(exe.B),len(exe.A))     )
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
          reshape(len(exe.B),len(exe.A))     )
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
          reshape(len(exe.B),len(exe.A))     )
    plt.colorbar(xy3imag)
    plt.title('Im$X_{3} Y_{3}^{*} $')
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
#    plt.axis([0,30,0,30])
    plt.savefig('xy3imag.png')
    plt.show()
    plt.close()
def plt_A_B_bsg(i,j):#A_B bsgamma
    mass_axis1,mass_axis2 = i,j
    print(ABarray4()[0],mass_axis1,len(ABarray4()))
    print(ABarray4()[1],mass_axis2,len(ABarray4()))
    resultb = []
    resultn = []
#B>Xs+gamma SECTION    
    for n in np.arange(0,len(ABarray4()) ):
        y3hdm= bsg.BR_B_Xs_gamma(mb,mw,mass_axis1,mass_axis2,\
                        exe.Y2(*ABarray4()[n] ), exe.complexyfunction(*ABarray4()[n] ),\
#                        [0.0],[0.0])
                        exe.Y3(*ABarray4()[n] ), exe.complexyfunction3(*ABarray4()[n] )) 
        resultb.append(y3hdm / (1e-4) )
    
#########
    bsgamma = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultb).flatten()  ,len(np.array(resultb).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,\
          levels = np.array([3.11,3.5,3.87]) )
#           levels = np.arange(3.0,4.2,0.2),\
#          colors = ['black','royalblue','purple','darkgreen',\
#                    'brown','red','gray','orange','pink'])
    plt.colorbar(bsgamma)
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$ in '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2) )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([0,2.0 * bsg.PI, - 1.4, - 0.2])
#    plt.axis([0,30,0,30])
    plt.savefig( str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2) +'bsg.png')
    plt.show()
    plt.close()
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
######################################################################
############################################
##PLOT SECTION
#numerical()
#Plot_3()
#Plot_4()
#Plot_5()
Plot4_8_9()
#plotfigure1()
#plot_Hp1_Hp2()
#plot_under_Heatherbasis(exe.i,exe.j,exe.k,exe.l)
#############################################
#starttime = exe.time.time()
def scanscanplot():#scan each parameter plots
    for m in exe.tbe:
        plot_under_deltascan(exe.i,m,exe.k, exe.l)
        plt.savefig(str(m) + '.png')
        plt.show()
        plt.close()
#scanscanplot()
#endtime = exe.time.time()
#print(endtime - starttime)
#plt_A_B_xy()#A_B_XY*
#for jjjj in np.arange(80,180,20):
#    for iiii in np.arange(180,750,100):
#        plt_A_B_bsg(jjjj,iiii )#A_B_bsgamma
##############################################################
#print('Newcp-asymmetry',bsg.newa_cp(mb,mw,80,400,\
#                        exe.Y2(*ABarray4()[0] ),exe.complexyfunction(*ABarray4()[0] ),\
#                        exe.Y3(*ABarray4()[0] ),exe.complexyfunction3(*ABarray4()[0] )))
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
def figure4_plot():#NEDM figure 4
    ################### First charged Higgs
    masss = 200 # first charged Higgs
    ################### Second charged Higgs
    masss2 = 500
    #################### mass different
    massdiffer = 20
    massdiffer_list = np.array([ i for i in np.arange(20,80,20)] )
    m_axis = np.array([ i for i in np.arange(100,550,50)] )
    m1_axis = np.array([ i for i in np.arange(50,550,50)] )
    xy_axis = np.array([ i for i in np.arange(0,1.1,0.1)] )
    xy1_axis = np.array([ i for i in np.arange(-5,5.1,0.1)] )
    longarray = np.array([i for i in np.arange(-1,1.2,0.2)])
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
    for m in longarray:
        empty3hdm = []
        for j in xy1_axis:
            for i in m1_axis:
            # 1GeV =  5.06e13 cm^{-1}
                nedm = abs(dn(i,masss2,[complex(0,j)],[complex(0,m)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
            
                empty3hdm.append(nedm )
#            print('nedm',i,j,nedm )
        result8 = plt.contourf(m1_axis, xy1_axis, \
                          np.resize(np.array(empty3hdm),len(np.array(empty3hdm) )).\
                          reshape(len(xy1_axis),len(m1_axis)), \
                          levels = np.array([0.0,1.5e-26,3.3e-26])      
                          )
        plt.colorbar(result8)
        plt.xlabel('$M_{H^{\pm}_{1}}$')
        plt.ylabel('IM($X_2Y_2^*$)')
        plt.title('NEDM in 3HDM')
        plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
        plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
        plt.axis([50,500, 0, 1.0])
        plt.show()
        plt.savefig('m1' + str(masss2) + 'imxy2'+ str("%.2f" % m) + '.png' )
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
    
    plt.xlabel('IM($X_1Y_1^*$)')
    plt.ylabel('IM($X_2Y_2^*$)')
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
            nedm4 = abs(dn(i,j,[complex(1,0.3)],[complex(1,-0.4)]) / (5.06e13) )# dn [GeV^{-1}] > # cm
            
            empty4.append(nedm4 )
#            print('nedm',i,j,nedm )
    result4 = plt.contourf(m1_axis, m1_axis, \
           np.resize(np.array(empty4),len(np.array(empty4) )).\
           reshape(len(m1_axis),len(m1_axis)), \
           levels = np.array([0.0,1.5e-26,3.3e-26])      
           )
    plt.colorbar(result4)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('NEDM in 3HDM' )
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([50,500,50,500])
    plt.show()
    plt.close()
    return

def nedm3hdm_plot():#NEDM figure 4 3hdm
    resultxy2imag = []
    resultxy3imag = []
    for n in np.arange(0,len(ABarray4()) ):
        resultxy2imag.append( (  exe.complnedm2(*ABarray4()[n] ))  )
        resultxy3imag.append( (  exe.complnedm3(*ABarray4()[n] ))  )
#        print('xy2',(- exe.complexyfunction(*ABarray4()[n] )).imag)
#        print('xy3',(- exe.complexyfunction3(*ABarray4()[n] )).imag)
#    print('XY2',resultxy2imag)
#    print('xy3',resultxy3imag)
    nedm = abs(dn(100,120,resultxy2imag,resultxy3imag) / (5.06e13) )# dn [GeV^{-1}] > # cm
        
    print('empty',len(nedm),np.array(nedm).flatten())
    result = plt.contourf(exe.A, exe.B, \
                         np.resize(np.array(nedm).flatten(),len(np.array(nedm).flatten())).\
                                  reshape(len(exe.B),len(exe.A)), \
                                 levels = np.array([0.0,3.3e-26])        
                                  )
    plt.colorbar(result)
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    plt.title('NEDM in 3HDM')
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([80,500, 0.0, 0.5])
#    plt.show()
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
               levels = np.array([0,3.3]),      
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
               levels = np.array([0.0,3.3]),      
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
               levels = np.array([0,3.3]),      
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
               levels = np.array([0.0,3.3]),      
               colors='red',labels='$M_{H^{\pm}_{2}}$ = $M_{H^{\pm}_{1}}$ + 80')
    plt.title('NEDM with $M_{H^{\pm}_{1}}$ = '+ str("%3d" % masss) \
              + ', $M_{H^{\pm}_{2}}$ = '+ str("%3d" % (masss +  masssdif [3])) )
    plt.xlabel('IM($X_2Y_2^*$)')
    plt.ylabel('IM($X_3Y_3^*$)')
    plt.show()
    plt.close()
    return
###############################################################
#Neutron EDM PLOT SECTION
#    for m in np.arange(0,len(ABarray4()) ):
#        nedm3hdm= C_wmuh(mb,mass_axis[0],mass_axis[1],\
#                        exe.complexyfunction(*ABarray4()[m] ),\
#                        exe.complexyfunction3(*ABarray4()[m] )) 
#        result1.append(nedm3hdm )
#    nedm = plt.contourf(exe.A, exe.B, \
#           np.resize(np.array(result1).flatten()  ,len(np.array(result1).flatten() ) ).\
#          reshape(len(exe.B),len(exe.A)) ,\
#           levels = np.arange(2.0,7.0,1.0),\
#          colors = ['black','royalblue','purple','darkgreen',\
#                    'brown','red','gray','orange','pink'])
#    plt.colorbar(nedm)
#    plt.title('Neutron EDM in ['\
#                    + str("%02d" % mass_axis[0]) +',' + str("%02d"% mass_axis[1]) +']')
#    plt.xlabel(exe.readlist[int(exe.read1)])
#    plt.ylabel(exe.readlist[int(exe.read2)])
#    plt.show()
#    plt.close()    
#figure4_plot()#Neutron EDM plot
#nedm3hdm_plot()#Neutron EDM in 3hdm-plot with exercise file
#nedm3hdm_xy1_xy2()# NEDM figure 4 3hdm