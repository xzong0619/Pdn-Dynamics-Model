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

'''
fit_int_flag, explained: 
    
True - fit intercept automatically from LASSO regression
     - the intercept is obtained from LASSO.intercept_
False- fit intercept by adding additional first column of pi matrix 
     - the intercept is obtained from the first cluster interaction
     - the signficant cluster interactions need to be adjusted based on index
'''

fit_int_flag = False

with open('ES_iso.json') as f:
    ES_data = json.load(f)
    
Ec = ES_data['E_iso']
config = ES_data['config_iso']

# the number of Pd atoms in each structure
NPd_list = np.array([len(x) for x in config])


[Gm, Gsv, Gcv1, Gcv2, Gcv3, Gcv4] = pickle.load(open("clusters.p", "rb"))
Gcv = Gcv1+Gcv2+Gcv3

# Prepare for the pi matrix - X
x = np.load('pi4.npy')
X = x
# Add the first column of ones to pi matrix to fit intercept manually 
if not fit_int_flag:
    X = np.ones((x.shape[0], x.shape[1]+1)) #the first column of pi matrix is set a 1, to be the intercept
    X[:,1:] = x     

# Prepare for the true output values - y
y = np.array(Ec)
   
#%% LASSO Regression 

# save 10% of data point to the test set
X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.1, random_state=0)

NPd_test = []
NPd_train = []
for i in y_test: NPd_test.append(NPd_list[np.where(y==i)[0][0]])
for i in y_train: NPd_train.append(NPd_list[np.where(y==i)[0][0]])
                               
NPd_test = np.array(NPd_test)
NPd_train = np.array(NPd_train)                              

# Cross-validation scheme                                  
rkf = RepeatedKFold(n_splits = 5, n_repeats = 2, random_state=0)

lasso_cv  = LassoCV(cv = rkf, max_iter = 1e7, tol = 0.0001, fit_intercept=fit_int_flag, random_state=5)
lasso_cv.fit(X_train, y_train)

# the optimal alpha
alpha = lasso_cv.alpha_
# the alpha path
alphas = lasso_cv.alphas_

alphas, coef_path, _ = lasso_path(X_train, y_train, alphas = lasso_cv.alphas_, fit_intercept=fit_int_flag)

# Coefficients for each term
coefs = lasso_cv.coef_
# The original intercepts 
intercept = lasso_cv.intercept_

# Access the errors 
y_predict_test = lasso_cv.predict(X_test)
y_predict_train = lasso_cv.predict(X_train)

RMSE_test = np.sqrt(sum(((y_test - y_predict_test)/NPd_test)**2)/len(y_test))
RMSE_train = np.sqrt(sum(((y_train - y_predict_train)/NPd_train)**2)/len(y_test))

RMSE_path = np.sqrt(np.mean(lasso_cv.mse_path_, axis = 1))

#%% Select the significant cluster interactions 
# The tolerance for zero coefficients
Tol = 1e-7
# The indices for non-zero coefficients/significant cluster interactions 
J_index = np.where(abs(coefs)>=Tol)[0]
# The number of non-zero coefficients/significant cluster interactions  
n_nonzero = len(J_index)
# The values of non-zero coefficients/significant cluster interactions  
J_nonzero = coefs[J_index] 
pi_nonzero = X[:, J_index]

# The number of nonzero coefficients in the trajectory
n_nonzero_v = []
for i in range(coef_path.shape[1]):
    n_nonzero_v.append(len(np.nonzero(coef_path[:,i])[0]))

# Pick the significant clusters
Gcv_nonzero = []

# Adjust for the manual intercept fitting
if not fit_int_flag:
    
    intercept = J_nonzero[0]
    n_nonzero = n_nonzero - 1 
    J_nonzero = J_nonzero[1:]
    pi_nonzero = pi_nonzero[:,1:]
    for i in J_index[1:]:
        # take out the first one and adjust the indices by -1 
        Gcv_nonzero.append(Gcv[i-1]) 
else:       
    for i in J_index:
        Gcv_nonzero.append(Gcv[i]) 
#%%
'''
Save Gcv_nonzero and J_nonzero to pickle for further use
''' 
pickle.dump([Gcv_nonzero, J_nonzero, intercept, 
             RMSE_test, RMSE_train, pi_nonzero], open('lasso.p','wb'))
    
    
#%% Plot the regression results
    
def predict_y(x, intercept, J_nonzero):
    
    # x is the column in pi matrix or the pi matrix 
    y = np.dot(x, J_nonzero) + intercept
    # the results should be the same as y = lasso_cv.predict(X)
    return y
    
def plot_lasso():
    
    '''
    #plot alphas vs RMSE along the path
    '''
    font = {'family' : 'normal', 'size'   : 15}
    
    fig = plt.figure(figsize=(6, 4))
    
    plt.plot(-np.log10(alphas),np.log10(RMSE_path),
             label='Average across the folds', linewidth=2)  
    plt.axvline(-np.log10(alpha), linestyle='--' , color='r', linewidth=3,
                label='alpha: CV estimate') 
    
    plt.legend(frameon=False,loc='best')
    plt.xlabel(r'$-log10(\lambda)$')
    plt.ylabel("log10(RMSE (eV))")    
    plt.tight_layout()
    matplotlib.rc('font', **font)
    fig.savefig('lasso_a_vs_cv.png')
    plt.show()   
       
    
    '''
    #plot alphas vs the number of nonzero coefficents
    '''
    
    fig = plt.figure(figsize=(6, 4))
    
    plt.plot(-np.log10(alphas),n_nonzero_v,
             label='Average across the folds', linewidth=2)     
    plt.axvline(-np.log10(alpha), linestyle='--' , color='r', linewidth=3,
                label='alpha: CV estimate') 
    plt.legend(frameon=False, loc='best')
    plt.xlabel(r'$-log10(\lambda)$')
    plt.ylabel("Number of Nonzero Coefficients ")    
    plt.tight_layout()
    matplotlib.rc('font', **font)
    fig.savefig('lasso_a_vs_n.png')
    plt.show() 
    
    
    '''
    #plot parity plot
    '''
    y_predict_all = lasso_cv.predict(X)
    #y_predict_all = predict_y(pi_nonzero, intercept, J_nonzero)
    
    plt.figure(figsize=(6,4))
    
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

    matplotlib.rc('font', **font)
    plt.tight_layout()
    fig.savefig('lasso_parity.png')
    plt.show()
    

#%% Ridge regression
from sklearn.linear_model import RidgeCV

# Compute paths

n_alphas = 2
alphas = np.logspace(-5, 0, n_alphas)

coef_path_ridge = []
for a in alphas:
    ridge_cv = RidgeCV(alpha=a, cv = rkf, max_iter = 1e7, tol = 0.0001, fit_intercept=fit_int_flag, random_state=5)
    ridge_cv.fit(X_train, y_train)
    coef_path_ridge.append(ridge_cv.coef_)
    cvm = ridge_cv.cv_values_ 























