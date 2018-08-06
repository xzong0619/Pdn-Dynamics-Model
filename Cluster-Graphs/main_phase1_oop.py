# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 11:39:26 2018

@author: wangyf
"""

'''
main phase 1 in OOP version
'''


import lat_fun as lf 
import numpy as np

print('Hello world')
mother = [(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5),
          (1/2, -3**0.5/2), (3/2, -3**0.5/2), (0, -3**0.5),
          (-1,0), (-1/2, -3**0.5/2), (-3/2, -3**0.5/2),
          (-1/2, 3**0.5/2), (-3/2, 3**0.5/2)]

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
'''
Energy of clusters in eV
'''

Ec = [0, -0.45, -1.5, -2.56, -3.41, -4.49, -4.38, -4.32, -6.07, -4.9, -5.11, -5.12]


     
empty = 'grey'
filled = 'r'
occ = [empty, filled]

Clusters = lf.clusters(occ)
Clusters.get_mother(mother)
Gm = Clusters.mother
  
'''
Creat 12 configurations
'''
Clusters.get_configs(config)
Gsv = Clusters.Gsv


