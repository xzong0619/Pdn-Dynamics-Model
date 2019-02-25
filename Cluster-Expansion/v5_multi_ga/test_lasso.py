# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 16:37:24 2018

@author: wangyf
"""

from sklearn.linear_model import LassoCV, lasso_path, Lasso
from sklearn.model_selection import RepeatedKFold, cross_val_score, train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_validate

import numpy as np
import pickle 
import json

import matplotlib.pyplot as plt 
import matplotlib

font = {'family' : 'normal', 'size'   : 12}
matplotlib.rc('font', **font)


#%% User define functions 
'''
Plot the regression results
'''
def cal_RMSE_path(alphas, model):
    
    MSE_path = []
    
    for i, ai in enumerate(alphas):
        print('{} % done'.format(100*(i+1)/len(alphas)))
        estimator = model(alpha = ai,  max_iter = 1e7, tol = 0.01, fit_intercept=fit_int_flag, random_state=5)
        scores  = cross_val_score(estimator, X_train, y_train, cv = rkf, scoring=('neg_mean_squared_error'))
        test_scores = np.abs(scores)
        MSE_path.append(test_scores)
    
    RMSE_path = np.sqrt(np.array(MSE_path))
    
    return RMSE_path

      
def predict_y(x, intercept, J_nonzero):
    
    # x is the column in pi matrix or the pi matrix 
    y = np.dot(x, J_nonzero) + intercept
    # the results should be the same as y = lasso_cv.predict(X)
    return y

def plot_ridge_path(alpha, alphas, RMSE_path, model_name):
    
    fig = plt.figure(figsize=(6, 4))
    
    #plt.plot(-np.log10(alphas), np.log10(RMSE_path), ':', linewidth= 0.8)
    plt.plot(-np.log10(alphas), np.mean(RMSE_path, axis = 1), 
             label='Average across the folds', linewidth=2)  
    plt.axvline(-np.log10(alpha), linestyle='--' , color='r', linewidth=3,
                label='alpha: CV estimate') 
    
    plt.legend(frameon=False,loc='best')
    plt.xlabel(r'$-log10(\lambda)$')
    plt.ylabel("RMSE/cluster(ev)")    
    plt.tight_layout()
    fig.savefig(model_name + '_a_vs_cv.png')
    plt.show()   
    
    
def plot_path(alpha, alphas, RMSE_path, coef_path, model, model_name):
    
    '''
    Calculate the number of nonzero coefficient along the path
    '''
    # The number of nonzero coefficients in the trajectory
    n_nonzero_v = []
    for i in range(coef_path.shape[1]):
        n_nonzero_v.append(len(np.nonzero(coef_path[:,i])[0]))
     
    
    '''
    #plot alphas vs RMSE along the path
    '''
    
    fig = plt.figure(figsize=(6, 4))
    
    plt.plot(-np.log10(alphas), np.sqrt(RMSE_path), ':', linewidth= 0.8)
    plt.plot(-np.log10(alphas), np.mean(np.sqrt(RMSE_path), axis = 1), 
             label='Average across the folds', linewidth=2)  
    plt.axvline(-np.log10(alpha), linestyle='--' , color='r', linewidth=3,
                label='alpha: CV estimate') 
    
    plt.legend(frameon=False,loc='best')
    plt.xlabel(r'$-log10(\lambda)$')
    plt.ylabel("RMSE/cluster(ev)")    
    plt.tight_layout()
   
    fig.savefig(model_name + '_a_vs_cv.png')
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

    fig.savefig(model_name + '_a_vs_n.png')
    plt.show() 
    
    
    '''
    #plot parity plot
    '''
    y_predict_all = model.predict(X)
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

    plt.tight_layout()
    fig.savefig(model_name + '_parity.png')
    plt.show()


#%% Import data
    
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
   


#%% Preparation before regression
# Train test split, save 10% of data point to the test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=0)

NPd_test = []
NPd_train = []
for i in y_test: NPd_test.append(NPd_list[np.where(y==i)[0][0]])
for i in y_train: NPd_train.append(NPd_list[np.where(y==i)[0][0]])
                               
NPd_test = np.array(NPd_test)
NPd_train = np.array(NPd_train)                              

# The alpha grid used for plotting path
alphas_grid = np.logspace(0, -3, 20)

# Cross-validation scheme                                  
rkf = RepeatedKFold(n_splits = 5, n_repeats = 2 , random_state = 0)


#%% LASSO regression
'''   
# LassoCV to obtain the best alpha, the proper training of Lasso
'''
lasso_cv  = LassoCV(cv = rkf,  max_iter = 1e7, tol = 0.0001, fit_intercept=fit_int_flag, random_state=5)
lasso_cv.fit(X_train, y_train)

# the optimal alpha from lassocv
lasso_alpha = lasso_cv.alpha_
# Coefficients for each term
lasso_coefs = lasso_cv.coef_
# The original intercepts 
lasso_intercept = lasso_cv.intercept_

# Access the errors 
y_predict_test = lasso_cv.predict(X_test)
y_predict_train = lasso_cv.predict(X_train)

# error per cluster
lasso_RMSE_test = np.sqrt(mean_squared_error(y_test, y_predict_test))
lasso_RMSE_train = np.sqrt(mean_squared_error(y_train, y_predict_train))

# error per atom
lasso_RMSE_test_atom = np.sqrt(mean_squared_error(y_test/NPd_test, y_predict_test/NPd_test))
lasso_RMSE_train_atom = np.sqrt(mean_squared_error(y_train/NPd_train, y_predict_train/NPd_train))

# mse path is the test mse path
# RMSE_path = np.sqrt(lasso_cv.mse_path_.mean(axis = 1))

'''
#Use alpha grid prepare for lassopath
'''
#lasso_path to get alphas and coef_path
lasso_alphas, lasso_coef_path, _ = lasso_path(X_train, y_train, alphas = alphas_grid, fit_intercept=fit_int_flag)
lasso_RMSE_path = cal_RMSE_path(alphas_grid, Lasso)    
plot_path(lasso_alpha, lasso_alphas, lasso_RMSE_path, lasso_coef_path, lasso_cv, 'lasso')


#%% Select the significant cluster interactions 

# The tolerance for zero coefficients
Tol = 1e-7
# The indices for non-zero coefficients/significant cluster interactions 
J_index = np.where(abs(lasso_coefs)>=Tol)[0]
# The number of non-zero coefficients/significant cluster interactions  
n_nonzero = len(J_index)
# The values of non-zero coefficients/significant cluster interactions  
J_nonzero = lasso_coefs[J_index] 
pi_nonzero = X[:, J_index]

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

'''
Save Gcv_nonzero and J_nonzero to pickle for further use
''' 
pickle.dump([Gcv_nonzero, J_nonzero, intercept, 
             lasso_RMSE_test_atom, lasso_RMSE_train_atom, pi_nonzero], open('lasso.p','wb'))


    
#%% Ridge regression
        
from sklearn.linear_model import Ridge
'''   
# cal_RMSE_path to obtain the best alpha, the proper training of ridge
'''
ridge_RMSE_path = cal_RMSE_path(alphas_grid, Ridge)  
ridge_alpha =  alphas_grid[np.argmin(np.mean(ridge_RMSE_path, axis =1))]
plot_ridge_path(ridge_alpha, alphas_grid, ridge_RMSE_path, 'ridge')

ridge = Ridge(alpha = ridge_alpha,  fit_intercept=fit_int_flag)
ridge.fit(X_train, y_train)

# Access the errors 
y_predict_test = ridge.predict(X_test)
y_predict_train = ridge.predict(X_train)

ridge_RMSE_test = np.sqrt(mean_squared_error(y_test, y_predict_test))
ridge_RMSE_train = np.sqrt(mean_squared_error(y_train, y_predict_train))

ridge_RMSE_test_atom = np.sqrt(mean_squared_error(y_test/NPd_test, y_predict_test/NPd_test))
ridge_RMSE_train_atom = np.sqrt(mean_squared_error(y_train/NPd_train, y_predict_train/NPd_train))


#%% Elastic net regression
'''
l1 ratio = 1 - lasso, l1 
l1 ratio = 0 - ridge, l2
'''
from sklearn.linear_model import ElasticNetCV
      
def l1_enet(ratio):
    
    '''
    input l1 ratio and return the model, non zero coefficients and cv scores
    training elastic net properly
    '''
    enet_cv  = ElasticNetCV(cv = rkf, l1_ratio=ratio,  max_iter = 1e7, tol = 0.0001, fit_intercept=fit_int_flag, random_state=5)
    enet_cv.fit(X_train, y_train)
    
    # the optimal alpha
    enet_alpha = enet_cv.alpha_
    enet_coefs = enet_cv.coef_
    n_nonzero = len(np.where(abs(enet_coefs)>=1e-7)[0])
    # Access the errors 
    y_predict_test = enet_cv.predict(X_test)
    y_predict_train = enet_cv.predict(X_train)
    
    # error per cluster
    enet_RMSE_test = np.sqrt(mean_squared_error(y_test, y_predict_test))
    enet_RMSE_train = np.sqrt(mean_squared_error(y_train, y_predict_train))

    # error per atom
    enet_RMSE_test_atom = np.sqrt(mean_squared_error(y_test/NPd_test, y_predict_test/NPd_test))
    enet_RMSE_train_atom = np.sqrt(mean_squared_error(y_train/NPd_train, y_predict_train/NPd_train))

    return enet_alpha, n_nonzero, enet_RMSE_test, enet_RMSE_train, enet_RMSE_test_atom, enet_RMSE_train_atom

# The vector of l1 ratio
l1s = [.1, .5, .7, .9, .95, .99, 1]
enet = []
enet_alphas = []
enet_n  = []
enet_RMSE_test = []
enet_RMSE_train = []
enet_RMSE_test_atom = []
enet_RMSE_train_atom = []


for i, l1i in enumerate(l1s):
    print('{} % done'.format(100*(i+1)/len(l1s)))
    ai, n, RMSE_test, RMSE_train, RMSE_test_atom, RMSE_train_atom = l1_enet(l1i)
    enet_alphas.append(ai)
    enet_n.append(n)
    
    enet_RMSE_test.append(RMSE_test)
    enet_RMSE_train.append(RMSE_train)
    
    enet_RMSE_test_atom.append(RMSE_test_atom)
    enet_RMSE_train_atom.append(RMSE_train_atom)


#%% Plot elastic net results
# expand the vector, put the result of ridge to the first
l1_ratio_v = np.array([0] + l1s)
enet_n_v  = np.array([X.shape[1]] + enet_n)
enet_RMSE_test_v = np.array([ridge_RMSE_test] + enet_RMSE_test)
plt.figure(figsize=(6,4))
fig, ax1 = plt.subplots()
ax1.plot(l1_ratio_v, enet_RMSE_test_v, 'b-')
ax1.set_xlabel('L1 Ratio')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('RMSE/cluster(ev)', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(l1_ratio_v, enet_n_v, 'r')
ax2.set_ylabel('Number of Nonzero Coefficients', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
fig.savefig('elastic_net.png')





