# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 12:00:34 2018

@author: wangyf
"""

from structure_constants import ECO
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler 
from sklearn import linear_model 
from sklearn.model_selection import cross_val_predict 
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import PolynomialFeatures 
from sklearn.pipeline import Pipeline

import pandas as pd
import numpy as np
import matplotlib as mat

mat.rc('font', size=12)

X = np.load('Xmatrix.npy')
y = np.array([-2.4,
                -2.98,
                -3.06,
                -2.58,
                -3.03,
                -2.28,
                -2.2,
                -2.36,
                -2.03,
                -2.31,
                -2.29,
                -2.43,
                -2.38,
                -2.31,
                -2.25,
                -2.29,
                -2.33,
                -2.33,
                -2.08,
                -2])

#%% Plot the trend 

nPC = X.shape[1]


feature_dict = {0: '1st Coordination Number',
                1: '2nd Coordination Number',
                2: 'General Coordination Number',
                3: '1st Ce Coordination Number',
                4: '2nd Ce Coordination Number',
                5: 'General Coordination Number',
                6: 'Distance from the support',
                7: 'Bader Charge',
                8: 'Number of Pd atoms'}


with plt.style.context('seaborn-whitegrid'):
    
    for cnt in range(nPC):
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

    plt.bar(range(nPC), var_exp, alpha=0.5, align='center',
            label='individual explained variance')
    plt.step(range(nPC), cum_var_exp, where='mid',
             label='cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal components')
    plt.xticks(np.arange(nPC), 
                ['PC %i'%(w+1) for w in range(nPC)])
    plt.legend(loc='best')
    plt.tight_layout()

#%% PCA use sklearn
pca = PCA()    
Xpc = pca.fit_transform(X_std) 

eig_vals_p = pca.explained_variance_
eig_vecs_p = pca.components_  # eigenvector
var_exp_p = pca.explained_variance_ratio_

#%% Plot the normalized desciptor loading
'''
Need to be completed
'''
# select the number of PCs

pc = 4

ind = 0
yvals = []
ylabels = []
bar_vals = []
space = 0.3
   


descriptors = ['CN1', 'CN2', 'GCN', 'CeCN1', 'CeCN2', 'CeGCN', 'Z', 'q', 'nPd']
cm = ['r', 'coral', 'pink',  'orange', 'y', 'gold', 'lightblue', 'lime', 'grey']
fig = plt.figure(figsize=(10,6))


ax = fig.add_subplot(111)
n = len(descriptors)
width = (1 - space) / (len(descriptors))
indeces = np.arange(0, pc) + 0.5  

# Create a set of bars at each position
for i, pci in enumerate(eig_vecs_p[:pc]):
    
    vals = pci/np.sum(np.absolute(pci))
    
    pos = width*np.arange(n) + i 
    ax.bar(pos, vals, width=width, label=str(i+1), color = cm) 
    

    

        
linex = np.arange(np.arange(0, pc).min() -0.5  , np.arange(0, pc).max()+ 2)

ax.set_xticks(indeces)
ax.set_xticklabels(list(np.arange(0,pc)+1))
ax.set_ylabel("Normalized Descriptoor Loading")
ax.set_xlabel("Principal Component #")    
'''
ax.legend(descriptors, bbox_to_anchor = (1.05, 1),loc= 'upper left', prop={'size':10},frameon=False)  
leg = ax.get_legend()
for c in range(n):
    leg.legendHandles[c].set_color(cm[c])
'''
plt.plot(linex, linex*0, c = 'k', lw = 0.8)
plt.show()


 
#%% Regression
# Create linear regression object
pc = 4
Xreg = Xpc[:,:pc]
def fit_linear_regression(X, y, degree):
    return Pipeline([("polynomial_features", PolynomialFeatures(degree=degree,
                                                                include_bias=False)),
                     ("linear_regression", linear_model.LinearRegression())]
                    ).fit(X, y) 
        
degree = 2
regr = linear_model.LinearRegression()

# Fit
regr.fit(Xreg, y)

# Calibration
y_c = regr.predict(Xreg)

# Cross-validation
y_cv = cross_val_predict(regr, Xreg, y, cv=10)

# Calculate scores for calibration and cross-validation
score_c = r2_score(y, y_c)
score_cv = r2_score(y, y_cv)

# Calculate mean square error for calibration and cross validation
mse_c = mean_squared_error(y, y_c)
mse_cv = mean_squared_error(y, y_cv)

estimator  = fit_linear_regression(Xreg, y, degree)
regr_poly = estimator.named_steps['linear_regression']
coefs = regr_poly.coef_
poly = estimator.named_steps['polynomial_features']
terms = poly.get_feature_names(['x1','x2','x3','x4','x5','x6','x7'])

y_poly = estimator.predict(Xreg)
score_poly = r2_score(y_poly, y_c)
mse_poly = mean_squared_error(y_poly, y_c)

fig, ax = plt.subplots()
ax.scatter(y, y_poly, facecolor = 'r', s  = 60, edgecolors=(0, 0, 0))
ax.plot([y.min(), y.max()], [y.min(), y.max()], 'k--', lw=2)
ax.set_xlabel('Measured')
ax.set_ylabel('Predicted')
plt.show()

xi = np.arange(len(coefs))
fig, ax = plt.subplots()
plt.bar(xi, coefs)
linex = np.arange(xi.min()-1, xi.max()+2)
plt.plot(linex, linex*0, c = 'k')
plt.xticks(xi, terms, rotation=45 )
plt.ylabel("Regression Coefficient Value (eV)")
plt.xlabel("Regression Coefficient")  
plt.show()
