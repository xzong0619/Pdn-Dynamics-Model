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
l4 = np.array([l1[0]])

      
'''
Add z coordinate
'''
dz = 1
l1d = lf.add_z(l1, dz)  
l2d = lf.add_z(l2, 2*dz)
l3d = lf.add_z(l3, 3*dz)
l4d = lf.add_z(l4, 4*dz)



empty = 'grey'
filled = 'r'
occ = [empty, filled]

'''
only draw 1st nearest neighbors?
'''
NN1 = 0
mother  = np.concatenate((l1d,l2d,l3d,l4d), axis =0)

Clusters = lf.clusters(occ, NN1)
Clusters.get_mother(mother)
Gm = Clusters.mother

#%%
'''
Create Configurations
'''


