#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 11:39:02 2020

@author: muyuansong0419
"""
#########B>Sgamma related plots only
from invariant import *
import bsgamma as bsg
import numpy as np
import matplotlib.pyplot as plt
import exercise as exe
#####################################################################
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
######################################################################
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
##########################################
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
           levels = np.array([2.99,3.55])      
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
#########################################
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
########################################
##################################################################################
def plt_A_B_bsg(i,j):# Bsgamma-result in {A,B} plane
    mass_axis1,mass_axis2 = i,j
    print(ABarray4()[0],mass_axis1,len(ABarray4()))
    print(ABarray4()[1],mass_axis2,len(ABarray4()))
    resultb = []
#B>Xs+gamma SECTION    
    for n in np.arange(0,len(ABarray4()) ):
        y3hdm= bsg.BR_B_Xs_gamma(mb,mw,mass_axis1,mass_axis2,\
                        exe.Y2(*ABarray4()[n] ), exe.complexyfunction(*ABarray4()[n] ),\
                        exe.Y3(*ABarray4()[n] ), exe.complexyfunction3(*ABarray4()[n] )) 
        resultb.append(y3hdm / (1e-4) )
#########
    bsgamm = plt.contourf(exe.A, exe.B, \
           np.resize(np.array(resultb).flatten()  ,len(np.array(resultb).flatten() ) ).\
          reshape(len(exe.B),len(exe.A)) ,\
          levels = np.array([2.99,3.55]),colors = ['green'] )
    plt.colorbar(bsgamm)
    plt.title('BR($\\bar{B} \\to X_{s} \gamma$) in '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2) )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    plt.axis([0,6.5,-1.6,0])
#    plt.axis([1,60,-1.6,0])
#    plt.axis([0,60,0,60])
    plt.savefig( str("%02d" % mass_axis1) + str("%02d"% mass_axis2) +'bsg.png')
    plt.show()
    plt.close()
####################################################################################

def plt_A_B_cps(i,j):#CP-asymmetry in {A,B} plane
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
        cmap = plt.cm.get_cmap('RdBu_r') )# levels = np.array([-12,-10,-8,-6,-4,-2,0,2,4]) )
    plt.colorbar(cpresult)
    plt.title('$A_{CP}(B \\to X_{s}\gamma)$ with charged Higgs: '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2)+' GeV.' )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    plt.axis([0,6.5,-1.6,0])
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
          cmap = plt.cm.get_cmap('RdBu_r'))#levels = np.arange(-20,-8,2) )
    plt.colorbar(cpresult)
    plt.title('$A_{CP}(B \\to X_{d}\gamma)$ with charged Higgs: '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2)+' GeV.' )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    plt.axis([0,6.5,-1.6,0])
    plt.savefig( 'cpd'+ str("%02d" % mass_axis1) + str("%02d"% mass_axis2) +'.png')
    plt.show()
    plt.close()
def plt_A_B_untag(i,j):#Untag-asymmetry in {A,B} plane
    result_untag = []
    mass_axis1,mass_axis2 = i,j
    for n in np.arange(0,len(ABarray4()) ):
        untagg = bsg.untag_cp(mb,mw,mass_axis1,mass_axis2,\
                        exe.Y2(*ABarray4()[n] ), exe.complexyfunction(*ABarray4()[n] ),\
#                        [0.0],[0.0])
                        exe.Y3(*ABarray4()[n] ), exe.complexyfunction3(*ABarray4()[n] )) 
        result_untag.append(untagg)
    result = plt.contourf(exe.A, exe.B, \
                            np.resize(np.array(result_untag).flatten()  ,\
                            len(np.array(result_untag).flatten() ) ).\
                            reshape(len(exe.B),len(exe.A)), \
          cmap = plt.cm.get_cmap('RdBu_r'))#levels = np.arange(-20,-8,2) )
    plt.colorbar(result)
    plt.title('$A_{CP} (B \\to X_{s + d} \gamma )$ with charged Higgs: '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2)+' GeV.' )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    plt.axis([0,6.5,-1.6,0])
    plt.savefig( 'untag'+ str("%02d" % mass_axis1) + str("%02d"% mass_axis2) +'.png')
    plt.show()
    plt.close()
def plt_A_B_cpsdiffer(i,j):#Delta-CPS-asymmetry in {A,B} plane
    result_deltas = []
    mass_axis1,mass_axis2 = i,j
    for n in np.arange(0,len(ABarray4()) ):
        cpsdif = bsg.newdifferacps(mb,mw,mass_axis1,mass_axis2,\
                        exe.Y2(*ABarray4()[n] ), exe.complexyfunction(*ABarray4()[n] ),\
#                        [0.0],[0.0])
                        exe.Y3(*ABarray4()[n] ), exe.complexyfunction3(*ABarray4()[n] )) 
        result_deltas.append(cpsdif)
    result = plt.contourf(exe.A, exe.B, \
                            np.resize(np.array(result_deltas).flatten()  ,\
                            len(np.array(result_deltas).flatten() ) ).\
                            reshape(len(exe.B),len(exe.A)), \
          cmap = plt.cm.get_cmap('RdBu_r'))#levels = np.arange(-20,-8,2) )
    plt.colorbar(result)
    plt.title('$\\Delta_{X_s\gamma}$ with charged Higgs: '\
                    + str("%02d" % mass_axis1) +', ' + str("%02d"% mass_axis2)+' GeV.' )
    plt.xlabel(exe.readlist[int(exe.read1)])
    plt.ylabel(exe.readlist[int(exe.read2)])
    plt.axis([0,6.5,-1.6,0])
    plt.savefig( 'cpsdiffer'+ str("%02d" % mass_axis1) + str("%02d"% mass_axis2) +'.png')
    plt.show()
    plt.close()
#################################################################################
def mp1_mp2_cps(i,j,k,l):
    m1_axis = np.array([ i for i in np.arange(10,550,50)] )
    m2_axis = np.array([ i for i in np.arange(10,1050,50)] )
    m2 = m2_axis[0]
    m1 = m1_axis[0]
    empty = []
    for m2 in m2_axis:
        for m1 in m1_axis:
            acpp= bsg.newa_cp(mb,mw,m1,m2,\
                        [exe.Y2(i,j,k,l)],[- exe.X2(i,j,k,l) * np.conjugate(exe.Y2(i,j,k,l) )],\
                        [exe.Y3(i,j,k,l)],[- exe.X3(i,j,k,l) * np.conjugate(exe.Y3(i,j,k,l) )])
            empty.append(acpp)
    resultcps = plt.contourf(m1_axis, m2_axis, \
           np.resize(np.array(empty),len(np.array(empty))).\
           reshape(len(m2_axis),len(m1_axis)), \
#           colors = ['black','royalblue','purple','darkgreen','brown','red','gray','orange'],\
#           levels = np.array([-0.08,0.2])      
           )
    plt.colorbar(resultcps)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('$A_{CP} (b \\to s \gamma)$ for 3HDM')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([50,200, 50.0, 1000.0])
    plt.axis([0,500, 0.0, 1000.0])
    plt.show()
    plt.close()
    return
def mp1_mp2_cpsdiffer(i,j,k,l):# Delta_CP(B>Xs + gamma)
    print(i,j,k,l)
    m1_axis = np.array([ i for i in np.arange(10,550,50)] )
    m2_axis = np.array([ i for i in np.arange(10,1050,50)] )
    m2 = m2_axis[0]
    m1 = m1_axis[0]
    emptycpsd = []
    for m2 in m2_axis:
        for m1 in m1_axis:
            acppd= bsg.newdifferacps(mb,mw,m1,m2,\
                        [exe.Y2(i,j,k,l)],[ - exe.X2(i,j,k,l) * np.conjugate(exe.Y2(i,j,k,l) )],\
                        [exe.Y3(i,j,k,l)],[ - exe.X3(i,j,k,l) * np.conjugate(exe.Y3(i,j,k,l) )])
#            print(acppd)
            emptycpsd.append(acppd)
    resultcpsd = plt.contourf(m1_axis, m2_axis, \
           np.resize(np.array(emptycpsd),len(np.array(emptycpsd))).\
           reshape(len(m2_axis),len(m1_axis)), \
#           colors = ['black','royalblue','purple','darkgreen','brown','red','gray','orange'],\
#           levels = np.array([0.0,0.5,1,1.5,2])      
           )
    plt.colorbar(resultcpsd)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('$\\Delta_{X_s\gamma}$ for 3HDM')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([50,200, 50.0, 1000.0])
    plt.axis([0,500, 0.0, 500.0])
    plt.show()
    plt.close()
    return
def mp1_A_cpsdiffer(i,j,k,l):
    m1_axis = np.array([ i for i in np.arange(1,551,25)] )
#    m2_axis = np.array([ i for i in np.arange(10,1050,50)] )
    cps_li = []
    for j in exe.tbe:
        print('j',j)
        for m1 in m1_axis:
            acpp= bsg.newa_cp(mb,mw,m1,300,\
                        [exe.Y2(i,j,k,l)],[- exe.X2(i,j,k,l) * np.conjugate(exe.Y2(i,j,k,l) )],\
                        [exe.Y3(i,j,k,l)],[- exe.X3(i,j,k,l) * np.conjugate(exe.Y3(i,j,k,l) )])
            cps_li.append(acpp)
    resultcps = plt.contourf(m1_axis, exe.tbe, \
           np.resize(np.array(cps_li),len(np.array(cps_li))).\
           reshape(len(exe.tbe),len(m1_axis)), \
#           colors = ['black','royalblue','purple','darkgreen','brown','red','gray','orange'],\
#           levels = np.array([0,0.5,1,1.5,2])      
           )    
    plt.colorbar(resultcps)
    
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel(exe.readlist[int(exe.read1)])
    plt.title('$\\Delta_{X_s\gamma}$ for 3HDM')  
    plt.show()
    plt.close()
    return
def mp1_mp2_cpuntag(i,j,k,l):# untag_CP(B>X(s+d) + gamma)
    print(i,j,k,l)
    m1_axis = np.array([ i for i in np.arange(10,570,20)] )
    m2_axis = np.array([ i for i in np.arange(10,1070,20)] )
    m2 = m2_axis[0]
    m1 = m1_axis[0]
    emptytag = []
    for m2 in m2_axis:
        for m1 in m1_axis:
            tag= bsg.untag_cp(mb,mw,m1,m2,\
                        [exe.Y2(i,j,k,l)],[ - exe.X2(i,j,k,l) * np.conjugate(exe.Y2(i,j,k,l) )],\
                        [exe.Y3(i,j,k,l)],[ - exe.X3(i,j,k,l) * np.conjugate(exe.Y3(i,j,k,l) )])
#            print(acppd)
            emptytag.append(tag)
#    print(emptytag)
    resulttag = plt.contourf(m1_axis, m2_axis, \
           np.resize(np.array(emptytag),len(np.array(emptytag))).\
           reshape(len(m2_axis),len(m1_axis)), \
#           colors = ['black','royalblue','purple','darkgreen','brown','red','gray','orange'],\
          levels = np.array([0.4,0.6,0.8,1.2,1.4,1.6])      
           )
    plt.colorbar(resulttag)
    plt.xlabel('$M_{H^{\pm}_{1}}$')
    plt.ylabel('$M_{H^{\pm}_{2}}$')
    plt.title('$A_{CP} (B \\to X_{s + d} \gamma )$ for 3HDM')
#    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
#    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
#    plt.axis([50,200, 50.0, 1000.0])
    plt.axis([0,500, 0.0, 500.0])
    plt.show()
    plt.close()
    return
###############################################################################
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
#mp1_mp2_cps(exe.i,exe.j,exe.k,exe.l)# [MHP1,MHP2] for CP-asymmetry
#for jjjj in np.array([85,200,500,700]):
#        plt_A_B_bsg(85,jjjj )#A_B_bsgamma light light
#for iiii in np.array([85,120,140,200,600]):
#        plt_A_B_bsg(85,iiii )#A_B_bsgamma light heavy
#        plt_A_B_bsg(200,iiii)#A_B_bsgamma heavy heavy
#        plt_A_B_cps(85,iiii) #A_B_CP-asymmetry
#        plt_A_B_cpd(85,iiii)##A_B_CPD-asymmetry
#        plt_A_B_cpsdiffer(85,iiii) #Delta-CPS-asymmetry in {A,B} plane
#        plt_A_B_untag(85,iiii)#Untag-asymmetry in {A,B} plane
mp1_mp2_cps(exe.i,exe.j,exe.k,exe.l)
mp1_mp2_cpsdiffer(exe.i,exe.j,exe.k,exe.l)
#mp1_A_cpsdiffer(exe.i,exe.j,exe.k,exe.l)
mp1_mp2_cpuntag(exe.i,exe.j,exe.k,exe.l)