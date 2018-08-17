# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:52:08 2018

@author: wangyf
"""
from structure_constants import mother, config, Ec
import lattice_functions as lf
import numpy as np
import pandas as pd
import pickle


empty = 'grey'
filled = 'r'
occ = [empty, filled]

'''
only draw 1st nearest neighbors?
'''
NN1 = 0
'''
Draw mother/conifgurations/clusters?
'''
draw = [1, 0, 0]


Clusters = lf.clusters(occ, NN1, draw)
Clusters.get_mother(mother)
Gm = Clusters.Gm


#%%
'''
Create Configurations
'''
Clusters.get_configs(config)
Gsv = Clusters.Gsv


#%%
'''
Create clusters
'''
sub = lf.subgraphs(mother)
c1 = sub.get_s(1)
c2 = sub.get_s(2)
c3 = sub.get_s(3)


Clusters.get_clusters(mother, c1)
Gcv1 = Clusters.Gcv

Clusters.get_clusters(mother, c2)
Gcv2 = Clusters.Gcv

Clusters.get_clusters(mother, c3)
Gcv3 = Clusters.Gcv

pickle.dump([Gm, Gsv, Gcv1, Gcv2, Gcv3], open('clusters.p','wb'))

