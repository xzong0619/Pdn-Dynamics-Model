# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 14:36:10 2018

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

[Gm, Gsv, Gcv1, Gcv2, Gcv3] = pickle.load(open("clusters.p", "rb"))


#%% Stattistical analysis
'''
creat pi matrix
size of number of configuration * numbers of clusters
'''


Gcv = Gcv1+Gcv2+Gcv3
Cal = lf.calculations(occ)
pi =  Cal.get_pi_matrix(Gsv ,Gcv)      

#%%

np.save('pi3', pi, allow_pickle = True)
#pi = np.load('pi3.npy')
