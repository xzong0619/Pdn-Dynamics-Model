# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 13:36:33 2018

@author: wangyf
"""

import numpy as np
import lattice_functions as lf
import pandas as pd
#%%
'''
Mother graph
2D coorindates
'''

# first layer
l1 = np.array([(0, 0), (1, 0), (1 / 2, 3**0.5 / 2), (3 / 2, 3**0.5 / 2), (0, 3**0.5), (1, 3**0.5),
               (1 / 2, -3**0.5 / 2), (3 / 2, -3**0.5 / 2), (0, -3**0.5), (1, -3**0.5),
               (-1, 0), (-1 / 2, -3**0.5 / 2), (-3 / 2, -3**0.5 / 2), (-1, 3**0.5),
               (-1 / 2, 3**0.5 / 2), (-3 / 2, 3**0.5 / 2), (-1, -3**0.5)])
# second layer
l2 = np.array([np.sum(l1[[13, 14, 15]], 0),
               np.sum(l1[[2, 4, 14]], 0),
               np.sum(l1[[2, 3, 5]], 0),
               np.sum(l1[[0, 10, 14]], 0),
               np.sum(l1[[0, 1, 2]], 0),
               np.sum(l1[[10, 11, 12]], 0),
               np.sum(l1[[0, 6, 11]], 0),
               np.sum(l1[[1, 6, 7]], 0),
               np.sum(l1[[11, 8, 16]], 0),
               np.sum(l1[[8, 9, 6]], 0),
               np.sum(l1[[0, 2, 14]], 0),
               np.sum(l1[[0, 10, 11]], 0),
               np.sum(l1[[0, 1, 6]], 0),          
               ]) / 3

# third layer
l3 = np.array([np.sum(l1[[0, 2, 14]], 0),
               np.sum(l1[[0, 10, 11]], 0),
               np.sum(l1[[0, 1, 6]], 0),
               np.sum(l1[[6, 8, 11]], 0),
               np.sum(l1[[0]],0) ])/ 3
# fourth layer
l4 = np.array([l1[0]])

'''
Add z coordinate to make Pd atoms in tetrahedron
'''
dz = 6**0.5/3 
l1d = lf.add_z(l1, dz)
l2d = lf.add_z(l2, 2 * dz)
l3d = lf.add_z(l3, 3 * dz)
l4d = lf.add_z(l4, 4 * dz)

mother = np.concatenate((l1d, l2d, l3d, l4d), axis=0)

'''
Configuragtiions
'''

config = [[0],                      #Pd1
          [0,1],                    #Pd2
          [0,1,2],                  #Pd3
          [0,1,2,6],                #Pd4
          [0,1,2,21],
          [0,1,2,6,11],             #Pd5
          [0,1,2,6,29],
          [0,1,2,6,14,11],          #Pd6
          [0,1,2,6,11,3],
          [0,1,2,6,11,7],
          [0,1,2,6,11,29],
          [0,1,2,6,21,29],
          [0,1,2,6,11,10,14],        #Pd7
          [0,1,2,6,11,3,8],
          [0,1,2,6,14,11,4],
          [0,1,2,6,14,11,3],
          [0,1,2,6,14,11,21],
          [0,1,2,6,11,3,21],
          [0,1,2,6,11,7,29],
          [0,1,2,6,11,21,23],
          [0,1,2,6,11,14,21,23],     #Pd8
          [0,1,2,6,11,21,23,24],
          [0,1,2,6,11,7,21,23,24],   #Pd9
          [0,1,2,6,11,10,14,21,23],
          [0,1,2,6,11,10,14,21,28],
          [0,1,2,6,11,10,14,29,28,27], #Pd10
          [0,1,2,6,11,10,14,8,21,23],
          [0,1,2,6,11,7,21,23,24,32],
          [0,1,2,6,11,10,14,20,21,23,24], #Pd11
          [0,1,2,6,11,10,14,20,21,23,7],
          [0,1,2,6,11,10,14,20,21,23,8],
          [0,1,2,6,11,10,14,20,21,23,34],
          [0,1,2,6,11,10,14,20,21,23,24,7],  #Pd12
          [0,1,2,6,11,10,14,20,21,23,24,22],
          [0,1,2,6,11,10,14,20,21,23,24,7,22], #Pd13
          [0,1,2,6,11,10,14,20,21,23,24,22,18],
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12], #14
          [0,1,2,6,11,10,14,20,21,23,24,7,22,18],
          [0,1,2,6,11,10,14,20,21,23,24,22,18,34],
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18], #15
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18,4], #16
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18,4, 30], #17
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18,4, 34], 
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18,4, 30,31],#18
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18,4, 30,31,32],#19
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18,4, 30,31,32,35],#20
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18,4, 30,31,32,35, 33],#21
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18,4, 30,31,32,35, 25],
          [0,1,2,6,11,10,14,20,21,23,24,7,22,12,18,4, 30,31,32,35, 16]]


'''
Configuration Energy
'''
'''
Use pandas and copy clipboard in Excel to import data
Ec =list(pd.read_clipboard(header = None)['0'])
'''

Ec = [0.0,
 -0.45,
 -1.5,
 -2.56,
 -2.74,
 -3.41,
 -3.99,
 -4.49,
 -4.38,
 -4.32,
 -5.23,
 -4.95,
 -6.07,
 -4.9,
 -5.11,
 -5.12,
 -6.26,
 -5.96,
 -5.91,
 -6.65,
 -7.77,
 -7.51,
 -9.31,
 -9.19,
 -9.05,
 -10.98,
 -9.97,
 -9.95,
 -11.88,
 -11.64,
 -11.48,
 -11.2,
 -13.33,
 -13.01,
 -14.41,
 -14.15,
 -15.79,
 -15.62,
 -15.16,
 -16.9,
 -18.14,
 -18.9,
 -18.82,
 -20.24,
 -21.9,
 -22.7,
 -23.49,
 -23.36,
 -23.08]

#%% CO adsorption regression
'''
CO adsorbed structure index
'''
COindex = [1,
 2,
 3,
 4,
 5,
 6,
 7,
 11,
 20,
 21,
 23,
 26,
 30,
 33,
 35,
 37,
 40,
 41,
 42,
 44,
 45,
 46,
 46,
 46]

COindex = list(np.array(COindex) - 1)

'''
CO adsorption energies
'''

ECO = [-2.4,
 -2.98,
 -3.06,
 -2.94,
 -2.58,
 -3.03,
 -2.23,
 -2.28,
 -2.2,
 -2.36,
 -2.03,
 -2.31,
 -2.29,
 -2.43,
 -2.38,
 -2.31,
 -2.25,
 -2.29,
 -2.33,
 -2.33,
 -2.08,
 -2.0,
 -1.91,
 -1.43]

'''
CO adsorption sites
'''

COsites = [[0],
           [0,1],
           [0,1],
           [0,1],
           [0,1,21],
           [0,1],
           [0,6,29],
           [1,6,29],
           [0,2,21],
           [2,21],
           [6,7,24],
           [2,14,27],
           [1,6,7],
           [6,7,24],
           [6,7,24],
           [6,7,24],
           [6,7,24],
           [6,7,24],
           [6,7,24],
           [6,7,24],
           [6,7,24],
           [6,7,24],
           [23,24],
           [7]]

