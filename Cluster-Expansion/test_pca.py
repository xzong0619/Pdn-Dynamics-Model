# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 12:00:34 2018

@author: wangyf
"""

from get_cn import df
from structure_constants import ECO
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler 
from sklearn import linear_model 
from sklearn.model_selection import cross_val_predict 
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import PolynomialFeatures 


import pandas as pd
import numpy as np


X = df.iloc[:,:].values
y = np.array(ECO)

#%% Plot the trend 
feature_dict = {0: '1st Coordination Number',
                1: '2nd Coordination Number',
                2: 'General Coordination Number',
                3: '1st Ce Coordination Number',
                4: '2nd Ce Coordination Number',
                5: 'General Coordination Number',
                6: 'Distance from the support'}

with plt.style.context('seaborn-whitegrid'):
    
    for cnt in range(7):
        plt.figure(figsize=(6, 4))
        plt.scatter(X[:,cnt],y)
        plt.xlabel(feature_dict[cnt])
        plt.ylabel('CO Adsorption Energy (eV)')
    plt.legend(loc='upper right', fancybox=True, fontsize=8)
    plt.tight_layout()
    plt.show()

#%% PCA 

X_std = StandardScaler().fit_transform(X)
cov_mat = np.cov(X_std.T)
eig_vals, eig_vecs = np.linalg.eig(cov_mat)

print('Eigenvectors \n%s' %eig_vecs)
print('\nEigenvalues \n%s' %eig_vals)

tot = sum(eig_vals)
var_exp = [(i / tot)*100 for i in sorted(eig_vals, reverse=True)]
cum_var_exp = np.cumsum(var_exp)
                        
                        
with plt.style.context('seaborn-whitegrid'):
    plt.figure(figsize=(6, 4))

    plt.bar(range(7), var_exp, alpha=0.5, align='center',
            label='individual explained variance')
    plt.step(range(7), cum_var_exp, where='mid',
             label='cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal components')
    plt.xticks(np.arange(7), 
                ['PC %i'%(w+1) for w in range(7)])
    plt.legend(loc='best')
    plt.tight_layout()

#%% PCA use sklearn
pca = PCA()    
Xreg = pca.fit_transform(X_std) 

eig_vals_p = pca.explained_variance_
eig_vecs_p = pca.components_  # eigenvector
var_exp_p = pca.explained_variance_ratio_

#%% Plot the normalized desciptor loading
'''
Need to be completed
'''

#%% Regression








