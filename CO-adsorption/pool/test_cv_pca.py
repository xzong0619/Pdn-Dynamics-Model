# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 17:43:46 2018

@author: wangyf
"""

'''
Perform Cross Validation PCA regression
Leave one out cross-validation 
Choose the number of PC with lowest test CV score
'''
from sklearn.model_selection import RepeatedKFold, cross_validate, LeaveOneOut
import sklearn.cross_validation as cv

nPC = 10
loo = LeaveOneOut() 
test_r2 = []
test_RMSEs = []
train_r2 = []
train_RMSEs = []

for nPCi in np.arange(0,nPC): 
    Xreg = Xpc[:,:nPCi+1]
    scores  = cross_validate(pc2_estimator, Xreg, y, cv=loo,
                                scoring=('r2', 'neg_mean_squared_error'),
                                return_train_score=True)
    
    train_scores = np.sqrt(np.abs(scores['train_neg_mean_squared_error'])) 
    train_score_mean = np.mean(train_scores)
    train_r2.append(np.mean(scores['train_r2']))
    
    test_scores = np.sqrt(np.abs(scores['test_neg_mean_squared_error'])) 
    test_score_mean = np.mean(test_scores)
    test_r2.append(np.mean(scores['test_r2']))
    
    train_RMSEs.append(train_score_mean)
    test_RMSEs.append(test_score_mean)
#%%
nplot = 8
plt.figure()    
plt.plot(np.arange(1,nPC+1)[:nplot], test_RMSEs[:nplot], c = 'r', label = 'Test')
plt.plot(np.arange(1,nPC+1)[:nplot], train_RMSEs[:nplot], c = 'b', label = 'Train')
plt.ylabel('RMSE(eV)')
plt.xlabel('nPC')
plt.xticks(np.arange(1,nplot+1))
plt.legend(loc= 'best', frameon=False)   
 
plt.figure() 
#plt.plot(np.arange(1,nPC+1)[:nplot], test_r2[:nplot], c = 'r', label = 'Test')
plt.plot(np.arange(1,nPC+1)[:nplot], train_r2[:nplot], c = 'b', label = 'Train')
plt.ylabel('r2')
plt.xlabel('nPC')
plt.xticks(np.arange(1,nplot+1))
plt.legend(loc= 'best', frameon=False)    