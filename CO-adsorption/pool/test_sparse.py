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


X, _ = make_friedman1(n_samples=100, n_features=10, random_state=0) 
SPCA= SparsePCA(alpha = 1, n_components=10, random_state=0)
SPCA.fit(X) 
X_spca = SPCA.transform(X)

# Covriance matrix from PCA
cov_pc = np.cov(X_spca.T) 
eig_vals, eig_vecs = np.linalg.eig(cov_pc)

eig_vecs_spc = SPCA.components_
nzeros = np.mean(SPCA.components_ == 0) 



print(nzeros)