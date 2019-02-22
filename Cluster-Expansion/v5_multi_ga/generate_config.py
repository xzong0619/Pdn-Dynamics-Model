# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:52:08 2018

@author: wangyf
"""

from structure_constants import mother, dz
import lattice_functions as lf
import numpy as np
import pandas as pd
import pickle
import json

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
Clusters.get_mother(mother, dz)
Gm = Clusters.Gm

with open('ES_iso.json') as f:
    ES_data = json.load(f)
    
Ec = ES_data['E_iso']
config = ES_data['config_iso']
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
sub = lf.subgraphs(mother, dz)
Gcv1 = sub.get_s2(1)
Gcv2 = sub.get_s2(2)
Gcv3 = sub.get_s2(3)
Gcv4 = sub.get_s2(4)

#%%
pickle.dump([Gm, Gsv, Gcv1, Gcv2, Gcv3, Gcv4], open('clusters.p','wb'))

