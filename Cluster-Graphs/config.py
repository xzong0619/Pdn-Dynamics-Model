# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:52:08 2018

@author: wangyf
"""
import mlattice as ml


# first layer
mother = [(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5), (1, 3**0.5), 
          (1/2, -3**0.5/2), (3/2, -3**0.5/2), (0, -3**0.5), (1, -3**0.5), 
          (-1,0), (-1/2, -3**0.5/2), (-3/2, -3**0.5/2), (-1, 3**0.5), 
          (-1/2, 3**0.5/2), (-3/2, 3**0.5/2), (-1, -3**0.5)]
# second layer

# third layer

# fourth layer
          



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

Gm = ml.gmothers(mother)