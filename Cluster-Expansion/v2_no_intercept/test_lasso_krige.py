# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:55:29 2018

@author: wangyf
"""
import sys
import os
import numpy as np
import pickle
from matplotlib import pyplot as plt

HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Google Drive (wangyf@udel.edu)','Udel', 'Research','Kriging', 'Code', 'pyKridging'))

from pyKriging.krige import kriging

[Gcv_nonzero, J_nonzero, 
 intercept, MSE_test, MSE_train,
 pi_nonzero,y] =  pickle.load(open("lasso.p", "rb"))

import reverse_graph as rg
import lattice_functions as lf
from config import Gsv


# one note is that when k>n, kriging will not work, pay attention to the results of LASSO
X = pi_nonzero # n*k matrix
y = y # n*1 vector 


#%%
'''
Do kriging
'''
k = kriging(X, y)
k.train(optimizer='ga')  
newpoints = k.infill(10, method = 'ei')

#%%
np.save('kriging_pts', newpoints)
#%% 
'''
'''
#Generate reverse graph
'''
empty = 'grey'
filled = 'r'
occ = [empty, filled]
Gcv = tlasso.Gcv_nonzero

rgf = rg.reverse(occ,Gcv)
rgf.get_supergraph(newpoints)

ncluster = rgf.ncluster

ncluster_int = np.round(ncluster).astype(int)
#%%
#com_list, com_list = rgf.get_all_com(Gcv, ncluster)
#%% take a test

testnc= rgf.test([Gsv[30]])




#%%
from scipy.special import comb

x = comb(39,19, exact = False)

'''