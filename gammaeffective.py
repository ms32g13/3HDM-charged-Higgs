#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 17:40:01 2019

@author: muyuansong0419
"""

import numpy as np
# gamma_matrices for anomalous dimension matrix (Wilson coefficients)
def gamma0eff():#gamma_0_eff_ji matrix values
    gamma = np.zeros((8,8))
    gamma[0][0] = -4
    gamma[0][1] = 8 / 3
    gamma[0][3] = - 2 / 9
    gamma[0][6] = - 208/243
    gamma[0][7] = 173/162
    gamma[1][0] = 12
    gamma[1][3] = 4/3
    gamma[1][6] = 416/81
    gamma[1][7] = 70/27
    gamma[2][3] = - 52 /3
    gamma[2][5] = 2
    gamma[2][6] = - 176/81
    gamma[2][7] = 14/27
    gamma[3][2] = - 40 /9
    gamma[3][3] = - 100 /9
    gamma[3][4] = 4 /9
    gamma[3][5] = 5/6
    gamma[3][6] = - 152 /243
    gamma[3][7] = - 587/162
    gamma[4][3] = - 256/3
    gamma[4][5] = 20
    gamma[4][6] = 4624/243
    gamma[4][7] = 4772/81
    gamma[5][2] = - 256 /9
    gamma[5][3] = 56/9
    gamma[5][4] = 40/9
    gamma[5][5] = - 2/3
    gamma[5][6] = 4624/243
    gamma[5][7] = 4772/81
    gamma[6][6] = 32/3
    gamma[7][6] = - 32/9
    gamma[7][7] = 28/3
    return gamma
def gamma1eff():#gamma_1_eff_ji matrix values
    gamma = np.zeros((8,8))
    gamma[0][0] = - 355/9
    gamma[0][1] = - 502/27
    gamma[0][2] = - 1412 / 243
    gamma[0][3] = - 1369/ 243
    gamma[0][4] = 134 /243
    gamma[0][5] = -35 / 162
    gamma[0][6] = - 818/243
    gamma[0][7] = 3779/324
    gamma[1][0] = - 35/3
    gamma[1][1] = - 28/3
    gamma[1][2] = - 416 /81
    gamma[1][3] = 1280/81
    gamma[1][4] = 56/81
    gamma[1][5] = 35/27
    gamma[1][6] = 508/81
    gamma[1][7] = 1841 / 108
    gamma[2][2] = - 4468 /81
    gamma[2][3] = - 31469/81
    gamma[2][4] = 400 /81
    gamma[2][5] = 3373 /108
    gamma[2][6] = 22348 /243
    gamma[2][7] = 10178/81
    gamma[3][2] = - 8158 /243
    gamma[3][3] = - 59399 / 243
    gamma[3][4] = 269/486
    gamma[3][5] = 12899 /648
    gamma[3][6] = - 17584 / 243
    gamma[3][7] = - 172471 / 648
    gamma[4][2] = - 251680 /81
    gamma[4][3] = - 128648 /81
    gamma[4][4] = 23836 / 81
    gamma[4][5] = 6106 /27
    gamma[4][6] = 1183696 / 729
    gamma[4][7] = 2901296 / 243
    gamma[5][2] = 58640 /243
    gamma[5][3] = - 26348 /243
    gamma[5][4] = - 14324 /243
    gamma[5][5] = - 2551/162
    gamma[5][6] = 2480344/2187
    gamma[5][7] = - 3296257 / 729
    gamma[6][6] = 4688 /27
    gamma[7][6] = - 2192 /81
    gamma[7][7] = 4063 /27
    return gamma