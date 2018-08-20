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

x = np.load('pi3.npy')
X = np.ones((x.shape[0], x.shape[1]+1))
X[:,1:] = x        
    
#%%

y = np.array(Ec)


X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.1, random_state=0)

rkf = RepeatedKFold(n_splits = 10, n_repeats = 5)

lasso_cv  = LassoCV(cv = rkf)
lasso_cv.fit(X_train, y_train)

#alphas, coef_path, _ = lasso_path(X_train, y_train, fit_intercept=False)
alpha = lasso_cv.alpha_
alphas = lasso_cv.alphas_
coefs = lasso_cv.coef_
y_predict = lasso_cv.predict(X_test)

MSE_predict = sum((y_test - y_predict))**2/len(y_test)
MSE_train = sum((y_train - lasso_cv.predict(X_train)))**2/len(y_train)

#%%
J_index = np.nonzero(coefs)[0] -1
J_nonzero = coefs[np.nonzero(coefs)]

