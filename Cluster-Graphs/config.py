# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:52:08 2018

@author: wangyf
"""
import Lattice_fun as Laf
import numpy as np

'''
2D coorindates
'''

# first layer
l1 = np.array([(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5), (1, 3**0.5), 
               (1/2, -3**0.5/2), (3/2, -3**0.5/2), (0, -3**0.5), (1, -3**0.5), 
               (-1,0), (-1/2, -3**0.5/2), (-3/2, -3**0.5/2), (-1, 3**0.5), 
               (-1/2, 3**0.5/2), (-3/2, 3**0.5/2), (-1, -3**0.5)])
# second layer
l2 =  np.array([np.sum(l1[[13,14,15]],0),
                np.sum(l1[[2,4,14]],0),
                np.sum(l1[[2,3,5]],0),                
                np.sum(l1[[0,10,14]],0),
                np.sum(l1[[0,1,2]],0),
                np.sum(l1[[10,11,12]],0),
                np.sum(l1[[0,6,11]],0),
                np.sum(l1[[1,6,7]],0),
                np.sum(l1[[11,8,16]],0),
                np.sum(l1[[8,9,6]],0)])/3
      
# third layer
l3 = np.array([np.sum(l1[[0,2,14]],0),
               np.sum(l1[[0,10,11]],0),
               np.sum(l1[[0,1,6]],0),
               np.sum(l1[[6,8,11]],0)])/3
# fourth layer
l4 = l1[0]

l1d = np.concatenate((l1, np.array([np.ones(len(l1))]).T),axis =1)         


'''
Add z coordinate
'''




config = [[0],
          [0,1],
          [0,1,2],
          [0,1,2,3],
          [0,1,2,3,5],
          [0,1,2,3,5,6],
          [0,1,2,3,5,9],
          [0,1,2,3,5,11],
          [0,1,2,5,9,8,11],
          [0,1,2,3,5,9,7],
          [0,1,2,5,8,11,12],
          [0,1,2,8,9,11,12]]


empty = 'grey'
filled = 'r'
occ = [empty, filled]


Gm = Laf.gmothers(l1, occ)
Gm = Laf.gmothers(np.concatenate((l1,l2,l3), axis =0), occ)