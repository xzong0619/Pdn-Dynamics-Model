# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 16:37:24 2018

@author: wangyf
"""

from sklearn.linear_model import LassoCV
from sklearn.linear_model import lasso_path
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import mean_squared_error
import sklearn.cross_validation as cv

import numpy as np
import pickle 
import json

import matplotlib.pyplot as plt 
import matplotlib

fit_int_flag = False

with open('ES_iso.json') as f:
    ES_data = json.load(f)
    
Ec = ES_data['E_iso']
config = ES_data['config_iso']

[Gm, Gsv, Gcv1, Gcv2, Gcv3, Gcv4] = pickle.load(open("clusters.p", "rb"))
Gcv = Gcv1+Gcv2+Gcv3

x = np.load('pi4.npy')
X = x

if not fit_int_flag:
    X = np.ones((x.shape[0], x.shape[1]+1)) #the first column of pi matrix is set a 1, to be the intercept
    X[:,1:] = x     
   
# the number of Pd atoms in each structure
NPd_list = np.array([len(x) for x in config])

#%%

y = np.array(Ec)


X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.1, random_state=0)

NPd_test = []
NPd_train = []
for i in y_test: NPd_test.append(NPd_list[np.where(y==i)[0][0]])
for i in y_train: NPd_train.append(NPd_list[np.where(y==i)[0][0]])
                               
NPd_test = np.array(NPd_test)
NPd_train = np.array(NPd_train)                              
                                 
rkf = RepeatedKFold(n_splits = 10, n_repeats = 10, random_state=0)

lasso_cv  = LassoCV(cv = rkf, max_iter = 1e7, tol = 0.0001, fit_intercept=fit_int_flag, random_state=5)
lasso_cv.fit(X_train, y_train)
alpha = lasso_cv.alpha_
alphas = lasso_cv.alphas_

alphas, coef_path, _ = lasso_path(X_train, y_train, alphas = lasso_cv.alphas_, fit_intercept=fit_int_flag)
#cv_scores = cv.cross_val_score(lasso_cv,X_train, y_train)


coefs = lasso_cv.coef_
y_predict_test = lasso_cv.predict(X_test)
y_predict_train = lasso_cv.predict(X_train)

MSE_test = sum(((y_test - y_predict_test)/NPd_test)**2)/len(y_test)
MSE_train = sum(((y_train - y_predict_train)/NPd_train)**2)/len(y_test)

MSE_path = np.mean(lasso_cv.mse_path_, axis = 1)
intercept = lasso_cv.intercept_

#%%
Tol = 1e-7
J_index = np.where(abs(coefs)>=Tol)[0]
#J_index = np.nonzero(coefs)[0] #-1 #for the matrix before we add column 1
n_coef = len(J_index)
J_nonzero = coefs[J_index] 
pi_nonzero = X[:, J_index]
n_nonzero = [] # number of nonzeor coefficients in the trajectory
for i in range(coef_path.shape[1]):
    n_nonzero.append(len(np.nonzero(coef_path[:,i])[0]))

Gcv_nonzero = []

if not fit_int_flag: 
    intercept = J_nonzero[0]
    J_nonzero = J_nonzero[1:]
    for i in J_index[1:]:
        Gcv_nonzero.append(Gcv[i]) 

else:       
    for i in J_index:
        Gcv_nonzero.append(Gcv[i]) 
#%%
'''
Save Gcv_nonzero and n_nonzero for further use
''' 
pickle.dump([Gcv_nonzero, J_nonzero, 
             intercept, MSE_test, MSE_train, 
             pi_nonzero, y], open('lasso.p','wb'))
    
#%%
    
def real_predict(x, intercept, J_nonzero):
    # x is the column in pi matrix 
    y = np.dot(x, J_nonzero) + intercept
    
    return y
    
def plot_lasso():
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
    plt.xlabel(r'$-log10(\lambda)$')
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
    
    plt.xlabel(r'$-log10(\lambda)$')
    plt.ylabel("Number of Nonzero Coefficients ")    
    
    matplotlib.rc('font', **font)
    
    plt.show() 
    plt.savefig('a_vs_n.png')
    
    '''
    #plot parity plot
    '''
    #y_predict_all = lasso_cv.predict(X)
    y_predict_all = real_predict(pi_nonzero, intercept, J_nonzero)
    
    plt.figure(figsize=(20,20))
    
    fig, ax = plt.subplots()
    ax.scatter(Ec,y_predict_all, s=60, facecolors='none', edgecolors='r')
    
    plt.xlabel("DFT Cluster Energy (eV)")
    plt.ylabel("Predicted Cluster Energy (eV)")
    
    
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


    