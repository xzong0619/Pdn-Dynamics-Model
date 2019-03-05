# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 14:37:54 2019

@author: wangyf
"""

'''
Test on sparse PCA
'''

import numpy as np
from sklearn.datasets import make_friedman1
from sklearn.decomposition import SparsePCA


#X, _ = make_friedman1(n_samples=100, n_features=10, random_state=0) 
#SPCA= SparsePCA(alpha = 1, n_components=10, random_state=0)
#SPCA.fit(X) 
#X_spca = SPCA.transform(X)
#
## Covriance matrix from PCA
#cov_pc = np.cov(X_spca.T) 
#eig_vals, eig_vecs = np.linalg.eig(cov_pc)
#
#eig_vecs_spc = SPCA.components_
#nzeros = np.mean(SPCA.components_ == 0) 
#
#
#
#print(nzeros)

# Number of data points
n = 3

# Number of parameters
p = 10

V1 = np.random.randn(n) * np.sqrt(290)
V2 = np.random.randn(n) * np.sqrt(300)
eps = np.random.randn(n) 
V3 = -0.3 * V1 + 0.925 * V2 + eps

X = np.zeros((n, p))
for i in [0,1,2,3]:
    X[:,i] = V1 + np.random.randn(n) 
for i in [4,5,6,7]:
    X[:,i] = V2 + np.random.randn(n) 
for i in [8,9]:
    X[:,i] = V3 + np.random.randn(n) 
    
    
    