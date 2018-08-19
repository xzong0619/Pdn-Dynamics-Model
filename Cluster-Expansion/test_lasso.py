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

X = np.load('pi3.npy')
y = np.array(Ec)


X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.1, random_state=0)


lasso_cv  = LassoCV()

rkf = RepeatedKFold(n_splits = 5, n_repeats = 1)

lasso_cv.fit(X_train, y_train)
alphas, coef_path, _ = lasso_path(X_train, y_train, fit_intercept=False)
alpha = lasso_cv.alpha_

y_predict = lasso_cv.predict(X_test)