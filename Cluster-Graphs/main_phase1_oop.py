# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 11:39:26 2018

@author: wangyf
"""

'''
main phase 1 in OOP version
'''


import lattice_functions as lf 
import numpy as np


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

'''
only draw 1st nearest neighbors?
'''
NN1 = 1
 
Clusters = lf.clusters(occ, NN1)
Clusters.get_mother(mother)
Gm = Clusters.mother
  
'''
Creat 12 configurations
'''
Clusters.get_configs(config)
Gsv = Clusters.Gsv

#%%
'''
Creat 7 clusters
'''

cmother = [(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5)]
ccluster = [[0], [0,1], [0,3], [1,4],[0,1,2],[1,2,4], [0,1,3]]

Clusters.get_clusters(cmother, ccluster)
Gcv = Clusters.Gcv

#%% Stattistical analysis
'''
creat pi matrix
size of number of configuration * numbers of clusters
'''
ns = len(config)  # number of configurations 

Cal = lf.calculations(occ)
J, pi =  Cal.get_J(Ec, Gsv ,Gcv)      
MSE = np.sum(np.power((np.dot(pi,J) - Ec),2))/ns 



#%%
'''
Calculate Leave One Out CV score
'''
Ec_predict = np.zeros(ns)


for i in range(ns):
    '''
    i is the index to be removed from the list
    '''
    Ec_LOO =  lf.LeaveOneOut(Ec,i)
    Gsv_LOO = lf.LeaveOneOut(Gsv,i)
    
    J_LOO = Cal.get_J(Ec_LOO, Gsv_LOO, Gcv)[0]
    Ec_predict[i] = np.dot(pi[i],J_LOO)
    
CV = np.sum(np.power((Ec_predict - Ec),2))/ns    





