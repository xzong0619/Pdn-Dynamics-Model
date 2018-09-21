# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 11:27:10 2018

@author: wangyf
"""

#test the reverse problem of pi matrix 

import lattice_functions as lf
import numpy as np
from itertools import combinations


l1 = np.array([(0, 0), (1, 0), (1 / 2, 3**0.5 / 2)])
l2 = np.array([np.sum(l1[[0, 1, 2]], 0)/3])
dz = 1 
l1d = lf.add_z(l1, dz)
l2d = lf.add_z(l2, dz*2)

mother = np.concatenate((l1d, l2d), axis=0)

config = [[2,1,3]]
Ec = [0]
empty = 'grey'
filled = 'r'
occ = [empty, filled]

'''
only draw 1st nearest neighbors?
'''
NN1 = 1
'''
Draw mother/conifgurations/clusters?
'''
draw = [1, 1, 0]


Clusters = lf.clusters(occ, NN1, draw)
Clusters.get_mother(mother, dz)
Gm = Clusters.Gm
Clusters.get_configs(config)
Gsv = Clusters.Gsv


sub = lf.subgraphs(mother, dz)
Gcv1 = sub.get_s2(1)
Gcv2 = sub.get_s2(2)
Gcv3 = sub.get_s2(3)

Gcv = Gcv1+Gcv2+Gcv3
niso = []
for Gi in Gcv:
    niso.append(len(Gi))        
niso = np.array(niso)   
    
Cal = lf.calculations(occ)

delta =  Cal.get_delta_l(Gsv[0] ,Gcv[0])  
pi = Cal.get_pi_matrix(Gsv, Gcv)

'''
Now let's do reverse engineering! 
'''

num = np.multiply(pi,niso)[0]

# Construct a new structure
new_nodes = []
for i, Gi in enumerate(Gcv):
    if Gi in Gcv1:
        n_nodes = int(num[i])
        print(n_nodes)
        new_nodes.append(Gi[:n_nodes])


