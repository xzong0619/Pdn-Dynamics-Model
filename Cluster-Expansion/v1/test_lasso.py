# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 16:37:24 2018

@author: wangyf
"""

from sklearn.linear_model import LassoCV
from sklearn.linear_model import lasso_path
from sklearn.model_selection import RepeatedKFold
import sklearn.cross_validation as cv

import numpy as np
from structure_constants import Ec
import pickle 

import matplotlib.pyplot as plt 
import matplotlib

[Gm, Gsv, Gcv1, Gcv2, Gcv3] = pickle.load(open("clusters.p", "rb"))
Gcv = Gcv1+Gcv2+Gcv3

x = np.load('pi3.npy')
X = np.ones((x.shape[0], x.shape[1]+1))
X[:,1:] = x        
    
#%%

y = np.array(Ec)


X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.1, random_state=0)

rkf = RepeatedKFold(n_splits = 10, n_repeats = 10)

lasso_cv  = LassoCV(cv = rkf)
lasso_cv.fit(X_train, y_train)
alpha = lasso_cv.alpha_
alphas = lasso_cv.alphas_

alphas, coef_path, _ = lasso_path(X_train, y_train, alphas = lasso_cv.alphas_, fit_intercept=False)
#cv_scores = cv.cross_val_score(lasso_cv,X_train, y_train)


coefs = lasso_cv.coef_
y_predict = lasso_cv.predict(X_test)

MSE_test = sum((y_test - y_predict)**2)/len(y_test)
MSE_train = sum((y_train - lasso_cv.predict(X_train))**2)/len(y_train)

MSE_path = np.mean(lasso_cv.mse_path_, axis = 1)

#%%
J_index = np.nonzero(coefs)[0] -1 #for the matrix before we add column 1
J_nonzero = coefs[np.nonzero(coefs)[0]] 
pi_nonzero = X[:, np.nonzero(coefs)[0]]
n_nonzero = []

Gcv_nonzero = []
for i in J_index:
    Gcv_nonzero.append(Gcv[i])
for i in range(coef_path.shape[1]):
    n_nonzero.append(len(np.nonzero(coef_path[:,i])[0]))
    
#%%
'''
'''
#plot alphas vs MSE along the path
'''

plt.figure(figsize=(40,40))
fig, ax = plt.subplots()

plt.plot(-np.log10(alphas),np.log10(MSE_path),
         label='Average across the folds', linewidth=2)  
plt.axvline(-np.log10(alpha), linestyle='--' , color='r', linewidth=3,
            label='alpha: CV estimate') 

plt.legend(frameon=False)
plt.xlabel("-log10(alphas)")
plt.ylabel("log10(Mean Square Error (eV))")    
font = {'family' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)

plt.show()   
plt.savefig('a_vs_cv.png')


'''
#plot alphas vs nonzero coefficents
'''

plt.figure(figsize=(20,20))
fig, ax = plt.subplots()

plt.plot(-np.log10(alphas),n_nonzero,
         label='Average across the folds', linewidth=2)     
plt.axvline(-np.log10(alpha), linestyle='--' , color='r', linewidth=3,
            label='alpha: CV estimate') 
plt.legend(frameon=False)

plt.xlabel("-log10(alphas)")
plt.ylabel("Number of Nonzero Coefficients ")    

matplotlib.rc('font', **font)

plt.show() 
plt.savefig('a_vs_n.png')

'''
#plot parity plot
'''
y_predict_all = lasso_cv.predict(X)


plt.figure(figsize=(20,20))

fig, ax = plt.subplots()
ax.scatter(Ec,y_predict_all, s=60, facecolors='none', edgecolors='r')

plt.xlabel("Eads DFT (eV)")
plt.ylabel("Eads Model Prediction(eV)")


lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
ax.set_xlim(lims)
ax.set_ylim(lims)

font = {'family' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)
plt.show()
plt.savefig('parity.png')


''' 
    