#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:27:11 2019

@author: muyuansong0419
b_sgamma.py plot_function
"""
from invariant import *
import b_sgamma as bsg
import numpy as np
import matplotlib.pyplot as plt
from exercise import xyfun,xyfun3,yfun,yfun3,U,X2,X3,Y2,Y3,Z2,Z3,\
        complexyfunction,complexyfunction3,read1,read2,A,B
########################
axs = np.arange(1, 61, 1)
print('A,B',A,B,read1,read2,len(B))
def ABarray4(): # using [A,B] to plot BR(B_bar > X_s +gamma)
    i = - bsg.PI/2.1
    j = 40.0  # tanbeta
    k = 40.0 #tangamma
    l = 0.0 #delta
    longlist = []
    reference_array = [i,j,k,l]
    for var_b in B:
            for var_a in A:
                my_tuple = ()
                for counter in range(4):
                    if int(read1) == counter:
                        my_tuple += (var_a,)
                    elif int(read2) == counter:
                        my_tuple += (var_b,)
                    else:
                        my_tuple += (reference_array[counter],)
                longlist.append(my_tuple)
    return longlist
print(ABarray4()[0],len(ABarray4()))
def array4():
    tuple_array = []
    for i in axs:
        my_tuple = ( - bsg.PI/2.1,i,60,0.0)
#        print(my_tuple)
        tuple_array.append(my_tuple)
    return tuple_array
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
def Plot_8_9():
    xim_axis = np.arange(-5,5.2,0.1)# figure 8
    XYimx_axis = [complex(-2,i)*1.0 for i in xim_axis]
    rangephi = np.arange(0,3.14,0.01) # figure 9
#    print('rangephi', rangephi,len(rangephi))
    XYexpim_axis = [ np.exp(complex(0,j)) for j in rangephi] 
#    print('REALX,IMX:',[np.complex(-2,i)*1.0 for i in xim_axis])
#    print('X = 2exp(i phi)',XYexpim_axis,len(XYexpim_axis))
   #mlx = MultipleLocator(1)
   #mly = MultipleLocator(0.25)
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
    plt.title('Figure8' )
    plt.legend(('$\mu = 4.8$ GeV', '$\mu = 2.4$ GeV', '$\mu = 9.6$ GeV '),
           loc='upper right', shadow=True,prop={'size': 8})
    plt.show()
    plt.close
    plt.plot(rangephi,y48phi_axis / (1e-4))
    plt.plot(rangephi,y24phi_axis / (1e-4))
    plt.plot(rangephi,y96phi_axis / (1e-4))
    plt.xlabel('$\\phi$')
    plt.ylabel('BR($\\bar{B} \\to X_{s} \gamma$) $\\times 10^{4}$')
    plt.title('Figure9' )
    plt.legend(('$\mu = 4.8$ GeV', '$\mu = 2.4$ GeV', '$\mu = 9.6$ GeV '),
           loc='upper right', shadow=True,prop={'size': 8})
    plt.xlim(0, 3)
    plt.ylim(0, 8.0)
    plt.show()
    plt.close
def Plot_4() :
    x_axis = np.array([ i for i in np.arange(100,1020,20)] )
#    print('II 2HDM tanbeta = 2',BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,\
#                                              [1.0/2.0],[1.0],[0],[0])\
#            )
    Y1_array = - 1.0/2.0 * np.cos(bsg.PI/4.0)
    X1_array =  1.0/ 2.0  * np.cos(bsg.PI/4.0) 
    X2_array = - 1.0/2.0 * np.sin(bsg.PI/4.0)
    Y2_array = 1.0/2.0 * np.sin(bsg.PI/4.0)
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
    plt.axis([100.0, 1000.0, 2.0, 4.0])
    plt.show()
    plt.close()
    return
def Plot_5() :
    x_axis = np.array([ i for i in np.arange(100,1020,20)] )
#    print('II 2HDM tanbeta = 2',BR_B_Xs_gamma(mb,mw,x_axis,x_axis + 20,\
#                                              [1.0/2.0],[1.0],[0],[0])\
#            )
    cpphase = np.exp(complex(0,PI))
    y22_2hdm = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                                      [1.0/2.0],[1.0],[0],[0]) )
    plt.plot(x_axis,y22_2hdm / (1e-4))
    for n in np.array([2,4,7,10,20]):
        
        X1_array =  - 2.0  * np.cos(- bsg.PI/4.0) - n / (1/np.sqrt(5)) * np.sin(- bsg.PI/4.0) \
                        * cpphase
        Y1_array = - 1.0/2.0 * np.cos(- bsg.PI/4.0)
        X2_array =  2.0  * np.sin(- bsg.PI/4.0) - n / (1/np.sqrt(5)) * np.cos(- bsg.PI/4.0) \
                        * cpphase
        Y2_array = 1.0/2.0 * np.sin(- bsg.PI/4.0)
        
        y12_3hdm = np.array(bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                     [Y1_array],[np.array(X1_array) * np.conjugate(Y1_array)],\
                     [Y2_array],[np.array(X2_array) * np.conjugate(Y2_array)]) )
        
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
def plot_under_Heatherbasis():
    x_axis = np.array([ i for i in np.arange(80,220,20)] )
    for n in np.arange(0,len(axs)):
        
        y_axis3hdm= bsg.BR_B_Xs_gamma(mb,mw,x_axis,x_axis + bsg.mass_differ,\
                        Y2(*array4()[n]),complexyfunction(*array4()[n]),\
                        Y3(*array4()[n]),complexyfunction3(*array4()[n])) 
        plt.plot(x_axis,y_axis3hdm / (1e-4))
    plt.axis([80,200, 3.0, 10.0])
    return
def plot_Hp1_Hp2():# [mH+1,mH+2] for fixed B_bar > X_s + gamma
#    fixedarray = (- bsg.PI/2.1,10,60,0.0) #mixing matrix parameters
    n = 5 #tangamma
    X1_array =  - 2.0  * np.cos(- bsg.PI/4.0) - n / (1.0/np.sqrt(5)) * np.sin(- bsg.PI/4.0)
    Y1_array = - 1.0/2.0 * np.cos(- bsg.PI/4.0)
    X2_array =  2.0  * np.sin(- bsg.PI/4.0) - n / (1.0/np.sqrt(5)) * np.cos(- bsg.PI/4.0)
    Y2_array = 1.0/2.0 * np.sin(- bsg.PI/4.0)
    m1_axis = np.array([ i for i in np.arange(80,220,7.5)] )
    m2_axis = np.array([ i for i in np.arange(100,1050,50)] )
    empty =[]
    m2 = m2_axis[0]
    m1 = m1_axis[0]
#    xx, yy = np.meshgrid(m1_axis, m2_axis)
    for m2 in m2_axis:
        for m1 in m1_axis:
            threehdm = bsg.BR_B_Xs_gamma(mb,mw,m1,m2,\
                        [Y1_array],[X1_array * np.conjugate(Y1_array)],\
                        [Y2_array],[X2_array * np.conjugate(Y2_array)]) 
            empty.append(threehdm)
    result = plt.contourf(m1_axis, m2_axis, \
           np.resize(np.array(empty) / (1e-4),len(np.array(empty) / (1e-4))).\
           reshape(len(m2_axis),len(m1_axis)), \
           colors = ['black','royalblue','purple','darkgreen','brown','red','gray','orange'])#,\
#           levels = np.array([2.99,3.85])         )
    plt.colorbar(result)
    plt.grid(axis='y', linestyle='-', color='0.75') # show y-axis grid line
    plt.grid(axis='x', linestyle='-', color='0.75') # show x-axis grid line
    plt.axis([80,200, 100.0, 300.0])
    plt.show()
    plt.close
    return
def plt_A_B_bsg():
    mass_axis = np.array([80.0,1000.0])
    result = []
    for n in np.arange(0,len(ABarray4()) ):
        y3hdm= bsg.BR_B_Xs_gamma(mb,mw,mass_axis[0],mass_axis[1],\
                        Y2(*ABarray4()[n] ),complexyfunction(*ABarray4()[n] ),\
                        Y3(*ABarray4()[n] ),complexyfunction3(*ABarray4()[n] )) 
        result.append(y3hdm / (1e-4) )
    y = plt.contour(A, B, \
           np.resize(np.array(result).flatten()  ,len(np.array(result).flatten() ) ).\
           reshape(len(B),len(A)) ,
           levels = np.arange(1.0,11.0,1.0), \
           colors = ['black','royalblue','purple','darkgreen','brown','red','gray','orange'])
    plt.colorbar(y)
    plt.show()
    plt.close
    return
def numerical():
    mass_axis = (80.0,250.0)
    result = []
    for n in np.arange(0,len(ABarray4()) ):
        y3hdm= bsg.BR_B_Xs_gamma(mb,mw,mass_axis[0],mass_axis[1],\
                        Y2(*ABarray4()[n] ),complexyfunction(*ABarray4()[n] ),\
                        Y3(*ABarray4()[n] ),complexyfunction3(*ABarray4()[n] )) 
#        print(y3hdm / (1e-4),n)
        result.append(y3hdm / (1e-4) )
    return np.concatenate(result).ravel()
######################################################################
#numerical()
#Plot_3()
#Plot_4()
#Plot_5()
#Plot_8_9()
#plot_Hp1_Hp2()
#plot_under_Heatherbasis()
plt_A_B_bsg()
#print(numerical())