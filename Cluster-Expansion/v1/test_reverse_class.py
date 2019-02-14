# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:40:06 2018

@author: wangyf
"""

import reverse_graph as rg
import lattice_functions as lf
import numpy as np


l1 = np.array([(0, 0), (1, 0), (1 / 2, 3**0.5 / 2)])
l2 = np.array([np.sum(l1[[0, 1, 2]], 0)/3])
dz = 1 
l1d = lf.add_z(l1, dz)
l2d = lf.add_z(l2, dz*2)

mother = np.concatenate((l1d, l2d), axis=0)

config = [[0,1]]
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

'''
Create Class object
'''

re_graph = rg.reverse(occ,Gcv)
re_graph.test(Gsv)
gs = re_graph.gsuper

'''
Draw the new configurations
'''
Clusters.get_configs(gs)
newGsv = Clusters.Gsv