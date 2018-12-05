# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 17:43:46 2018

@author: wangyf
"""

'''
Perform Cross Validation PCA regression
'''
from sklearn.model_selection import RepeatedKFold, cross_validate
import sklearn.cross_validation as cv

degree = 2

def fit_linear_regression(X, y, degree):
    '''
    # Create linear regression object
    '''
    return Pipeline([("polynomial_features", PolynomialFeatures(degree=degree,
                                                                include_bias=False)),
                     ("linear_regression", linear_model.LinearRegression())]
                    ).fit(X, y)  



estimator =  Pipeline([("polynomial_features", PolynomialFeatures(degree=degree,
                                                                include_bias=False)),
                     ("linear_regression", linear_model.LinearRegression())])
   
# Split data into test and train set
# We only touch the training set from now on
X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.1, random_state=0)


# PCA analysis 
X_std = StandardScaler().fit_transform(X_train)
pca = PCA()    
Xpc = pca.fit_transform(X_std) 
nPC =  Xpc.shape[1]
    
# Use repeated k fold cross validation
rkf = RepeatedKFold(n_splits = 10, n_repeats = 10)

# Make a loop to change the hyperparameter
# The hypermeter in our case is number of PCs
test_mses = []
train_mses = []

for nPCi in range(1,nPC):
    Xreg = Xpc[:,:nPCi]
    scores  = cross_validate(estimator, Xreg, y_train, cv=rkf,
                                scoring=('r2', 'neg_mean_squared_error'),
                                return_train_score=True)
    test_mses.append(-np.mean(scores['test_neg_mean_squared_error']))
    train_mses.append(-np.mean(scores['train_neg_mean_squared_error']))

plt.figure()    
plt.plot(range(1,nPC), test_mses, c = 'r')
plt.plot(range(1,nPC), train_mses, c = 'b')    

nPC = 3
Xreg = X_train[:,:nPC]
pcr_estimator  = fit_linear_regression(Xreg, y_train, degree)
y_predict_test = pcr_estimator.predict(X_test[:,:nPC])
mse_pca, score_pca = parity_plot(y_test, y_predict_test)