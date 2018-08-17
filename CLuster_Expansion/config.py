# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:52:08 2018

@author: wangyf
"""
from config_constants import mother, config, Ec
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

Clusters = lf.clusters(occ, NN1)
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
cclusters = c1


Clusters.get_clusters(mother, cclusters)
Gcv = Clusters.Gcv


#%% Stattistical analysis
'''
creat pi matrix
size of number of configuration * numbers of clusters
'''

ns = len(config)  # number of configurations 

Cal = lf.calculations(occ)
#
J, pi =  Cal.get_J(Ec, Gsv ,Gcv)      

pickle.dump([J, pi], open('dump.p','wb'))


#x = pickle.load(open("dump.p", "rb"))


#%%

np.save('pi3', pi, allow_pickle = True)
np.load('pi3.npy')


