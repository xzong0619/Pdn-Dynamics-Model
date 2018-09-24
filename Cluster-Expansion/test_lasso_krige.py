# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:55:29 2018

@author: wangyf
"""
import sys
import os
import numpy as np
from matplotlib import pyplot as plt

HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Google Drive (wangyf@udel.edu)','Udel', 'Research', 'Code', 'pyKridging'))

import pyKriging
from pyKriging.krige import kriging

import test_lasso as tlasso
import reverse_graph as rg
import lattice_functions as lf




X = tlasso.pi_nonzero
y = tlasso.y


#%%
'''
Do kriging
'''
k = kriging(X, y)
k.train(optimizer='ga')
newpoints = k.infill(1, method = 'ei')


#%% 
'''
Generate reverse graph
'''
occ = [empty, filled]
Gcv = tlasso.Gcv_nonzero
re_graph = rg.reverse(occ,Gcv)