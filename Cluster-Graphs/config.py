# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:52:08 2018

@author: wangyf
"""
import lat_fun as lf
import numpy as np

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
Add z coordinate
'''
dz = 1
l1d = lf.add_z(l1, dz)
l2d = lf.add_z(l2, 2 * dz)
l3d = lf.add_z(l3, 3 * dz)
l4d = lf.add_z(l4, 4 * dz)


empty = 'grey'
filled = 'r'
occ = [empty, filled]

'''
only draw 1st nearest neighbors?
'''
NN1 = 1
mother = np.concatenate((l1d, l2d, l3d, l4d), axis=0)

Clusters = lf.clusters(occ, NN1)
Clusters.get_mother(mother)
Gm = Clusters.Gm

#%%
'''
Create Configurations
'''
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
          

Clusters.get_configs(config)
Gsv = Clusters.Gsv
          
'''




