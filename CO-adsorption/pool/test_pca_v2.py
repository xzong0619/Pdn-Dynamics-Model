# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 23:26:24 2018

@author: yifan
"""

from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler 
from sklearn import linear_model 
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import PolynomialFeatures 
from sklearn.pipeline import Pipeline
from scipy.stats import norm
from sklearn.model_selection import RepeatedKFold, cross_validate


import pandas as pd
import numpy as np
import pickle

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt 
import matplotlib


matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['ytick.major.width'] = 2


[dem, Eads, descriptors, filename_list, sitetype_list] = pickle.load(open("pca_data.p", "rb"))
X = dem
y = Eads

#%% PCA parameters

nDescriptors = X.shape[1]
# select the number of PCs to plot in the bar graph
pc = min(7, len(descriptors))
# select the number of PCs to plot in regression
pc_reg = min(7, pc) 


#%% Plot the trend for each discriptor
def plot_discriptors():
    '''
    Plot the trend for each discriptor
    '''
    with plt.style.context('seaborn-whitegrid'):
        
        for cnt in range(nDescriptors):
            plt.figure(figsize=(6, 4))
            plt.scatter(X[:,cnt],y)
            plt.xlabel(descriptors[cnt])
            plt.ylabel('CO Adsorption Energy (eV)')
        plt.legend(loc='upper right', fancybox=True, fontsize=8)
        plt.tight_layout()
        plt.show()


def plot_discriptors_st():    
    '''
    Plot the trend for each discriptor with site type color coded
    '''
    with plt.style.context('seaborn-whitegrid'):
        for cnt in range(nDescriptors):
            plt.figure(figsize=(6, 4))
            for site, col in zip(('top', 'bridge', 'hollow'),
                        ('red', 'green', 'blue')):
                indices = np.where(np.array(sitetype_list) == site)[0]
                plt.scatter(X[:,cnt][indices],
                            y[indices],
                            label=site,
                            facecolor = col, 
                            alpha = 0.5,
                            s  = 60)
                
                plt.xlabel(descriptors[cnt])
                plt.ylabel('CO Adsorption Energy (eV)')
                
            plt.legend(bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)
            plt.tight_layout()
            plt.show() 

def plot_PdC_distances(): 
    '''
    Plot the histogram of PdC bond length
    '''
    PdCs = ['Pd1C', 'Pd2C', 'Pd3C'] 
    PdCi = []
    for PdC in PdCs: PdCi.append(descriptors.index(PdC))
    
    with plt.style.context('seaborn-whitegrid'):
        for cnt in PdCi:
            plt.figure(figsize=(8, 2))
            for site, col in zip(('top', 'bridge', 'hollow'),
                            ('red', 'green', 'blue')):
                indices = np.where(np.array(sitetype_list) == site)[0]
                plt.hist(X[:,cnt][indices],
                         bins= 160,
                         range= (1.5, 4),
                         label=site,
                         color = col,
                         alpha=0.5,)
            plt.xlabel(descriptors[cnt])
            plt.ylabel('Count')
            plt.xlim((1.5,4))
            plt.ylim((0,10))
            plt.legend(bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)
            plt.tight_layout()
            plt.show()


def plot_PdC1_distances(): 
    '''
    Plot the histogram of PdC bond length in a short x range
    '''
    PdCs = ['Pd1C', 'Pd2C', 'Pd3C'] 
    PdCi = []
    for PdC in PdCs: PdCi.append(descriptors.index(PdC))
    
    with plt.style.context('seaborn-whitegrid'):
        for cnt in PdCi:
            plt.figure(figsize=(6, 3))
            for site, col in zip(('top', 'bridge', 'hollow'),
                            ('red', 'green', 'blue')):
                indices = np.where(np.array(sitetype_list) == site)[0]
                plt.hist(X[:,cnt][indices],
                         bins= 60,
                         range= (1.6, 2.4),
                         label=site,
                         color = col,
                         alpha=0.5,)
            plt.xlabel(descriptors[cnt])
            plt.ylabel('Count')
            plt.xlim((1.6,2.4))
            plt.ylim((0,10))
            plt.legend(bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)
            plt.tight_layout()
            plt.show()

plot_discriptors_st()
plot_PdC_distances()
plot_PdC1_distances()

        
#%% PCA 
# Normalize the data
X_std = StandardScaler().fit_transform(X)
# Covriance matrix of original data
cov_mat = np.cov(X_std.T) 

# PCA use sklearn
pca = PCA()    
Xpc = pca.fit_transform(X_std) 
# Covriance matrix from PCA
cov_pc = np.cov(Xpc.T) 


# Plot Covariance structure 
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10,6))
ax1.set_title('X')
im1 = ax1.imshow(cov_mat)
fig.colorbar(im1, ax = ax1, shrink = 0.5)
ax2.set_title('X_PCA')
im2 = ax2.imshow(cov_pc)
fig.colorbar(im2, ax = ax2, shrink = 0.5)

#eig_vals, eig_vecs = np.linalg.eig(cov_mat)
eig_vals = pca.explained_variance_ #eigenvalues 
eig_vecs = pca.components_  # eigenvector
#tot = sum(eig_vals)
#var_exp = [(i / tot)*100 for i in sorted(eig_vals, reverse=True)]
var_exp = pca.explained_variance_ratio_ #explained variance ratio
cum_var_exp = np.cumsum(var_exp) #cumulative variance ratio
#print('Eigenvectors \n%s' %eig_vecs)
#print('\nEigenvalues \n%s' %eig_vals)

                      
with plt.style.context('seaborn-whitegrid'):
    plt.figure(figsize=(6, 4))

    plt.bar(range(nDescriptors), var_exp, alpha=0.5, align='center',
            label='individual explained variance')
    plt.step(range(nDescriptors), cum_var_exp, where='mid',
             label='cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal components')
    plt.xticks(np.arange(nDescriptors), 
                ['PC %i'%(w+1) for w in range(nDescriptors)])
    plt.legend(loc='best')
    plt.tight_layout()

#%% Plot the normalized desciptor loading
'''
Need to be completed
'''
ind = 0
yvals = []
ylabels = []
bar_vals = []
space = 0.3

cm = ['r', 'coral', 'pink',  'orange', 'y', 'gold', 'lightblue', 'lime', 'grey', 'green', 'brown'][:len(descriptors)]
fig = plt.figure(figsize=(10,6))


ax = fig.add_subplot(111)
n = len(descriptors)
width = (1 - space) / (len(descriptors))
indeces = np.arange(0, pc) + 0.5  

# Create a set of bars at each position
for i, pci in enumerate(eig_vecs[:pc]):
    
    vals = pci/np.sum(np.absolute(pci))
    
    pos = width*np.arange(n) + i 
    ax.bar(pos, vals, width=width, label=str(i+1), color = cm) 
        
linex = np.arange(np.arange(0, pc).min() -0.5  , np.arange(0, pc).max()+ 2)

ax.set_xticks(indeces)
ax.set_xticklabels(list(np.arange(0,pc)+1))
ax.set_ylabel("Normalized Descriptoor Loading")
ax.set_xlabel("Principal Component #")    
  
# Add legend using color patches
patches = []
for c in range(n):
    patches.append(mpatches.Patch(color=cm[c]))
plt.legend(patches, descriptors,
           bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)

plt.plot(linex, linex*0, c = 'k', lw = 0.8)
plt.show()

#%% Regression
def parity_plot_st(yobj,ypred, method):
    '''
    Plot the parity plot of y vs ypred
    return R2 score and MSE for the model
    colorcode different site types
    '''
    RMSE = np.sqrt(np.mean((yobj - ypred)**2))
    r2 = r2_score(yobj, ypred)
    fig, ax = plt.subplots()
    for site, col in zip(('top', 'bridge', 'hollow'),
                    ('red', 'green', 'blue')):
        indices = np.where(np.array(sitetype_list) == site)[0]
        ax.scatter(yobj[indices],
                    ypred[indices],
                    label=site,
                    facecolor = col, 
                    alpha = 0.5,
                    s  = 60)
    ax.plot([yobj.min(), yobj.max()], [yobj.min(), yobj.max()], 'k--', lw=2)
    ax.set_xlabel('DFT-Calculated ')
    ax.set_ylabel('Model Prediction')
    plt.title(r'{}, RMSE-{:.2}, $r^2$ -{:.2}'.format(method, RMSE, r2))
    plt.legend(bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)
    plt.show()
    
    return RMSE, r2

def error_distribution(yobj, ypred, method):
    
    '''
    Plot the error distribution
    return the standard deviation of the error distribution
    '''
    fig, ax = plt.subplots(figsize=(6,4))
    ax.hist(yobj - ypred,density=1, alpha=0.5, color='steelblue')
    mu = 0
    sigma = np.std(yobj - ypred)
    x_resid = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
    ax.plot(x_resid,norm.pdf(x_resid, mu, sigma), color='r')
    plt.title(r'{}, $\sigma$-{:.2}'.format(method, sigma))
    
    return sigma

def cross_validation(X, y, estimator): 
    '''
    Cross-validation
    '''
    rkf = RepeatedKFold(n_splits = 10, n_repeats = 10)
    scores  = cross_validate(estimator, X, y, cv=rkf,
                                scoring=('neg_mean_squared_error'),
                                return_train_score=True)
    # RMSE for repeated 10 fold test data 
    
    train_scores = np.sqrt(np.abs(scores['train_score'])) 
    train_score_mean = np.mean(train_scores)
    train_score_std = np.std(train_scores)
    
    test_scores = np.sqrt(np.abs(scores['test_score'])) 
    test_score_mean = np.mean(test_scores)
    test_score_std = np.std(test_scores)
    
    return [train_score_mean, train_score_std, test_score_mean, test_score_std]


def linear_regression(degree):
    '''
    # Create linear regression object
    '''
    return Pipeline([("polynomial_features", PolynomialFeatures(degree=degree,
                                                                include_bias=False)),
                     ("linear_regression", linear_model.LinearRegression())]

                    )

def regression_pipeline(X, y, estimator, method):

    y_predict = estimator.predict(X)
    RMSE, r2 = parity_plot_st(y, y_predict, method)
    scores = cross_validation(X, y, estimator)
    
    return RMSE, r2, scores
    
#def detect_outliers(y, y_predict, threshold = 0.5):
#    '''
#    detect the outlier
#    '''
#    diff = abs(y-y_predict)
#    outlier_index = np.where(diff>threshold)[0]
#    outliers_metal = []
#    outliers_surface = []
#    print('\nOutliers with AE in eV:')
#    for i in outlier_index: 
#        outliers_metal.append(data['Metal'][i])
#        outliers_surface.append(data['surface'][i])
#        print('{} {} - {:.2}'.format(outliers_metal[-1],outliers_surface[-1], diff[i]))
#        
#    return outliers_metal, outliers_surface

def ploy_coef(estimator, nreg):
    '''
    Plot the magnitude of each parameters
    '''
    regr = estimator.named_steps['linear_regression']
    coefs = regr.coef_
    intercept = regr.intercept_
    poly = estimator.named_steps['polynomial_features']
    feature_names = []
    for i in range(nreg): feature_names.append('x'+ str(i+1))
    terms = poly.get_feature_names(feature_names)
    
    
    xi = np.arange(len(coefs))
    fig, ax = plt.subplots()
    plt.bar(xi, coefs)
    linex = np.arange(xi.min()-1, xi.max()+2)
    plt.plot(linex, linex*0, c = 'k')
    plt.xticks(xi, terms, rotation=45, fontsize = 8 )
    plt.ylabel("Regression Coefficient Value (eV)")
    plt.xlabel("Regression Coefficient")  
    plt.show()
    
    return intercept, coefs

Xreg = Xpc[:,:pc_reg]
degree = 2

pc2_estimator  = linear_regression(degree)
pc2_estimator.fit(Xreg, y)
y_pc2 = pc2_estimator.predict(Xreg)

RMSE_pc2, r2_pc2, scores_pc2 = regression_pipeline(Xreg, y, pc2_estimator, 'PC2')
sigma_pc2 = error_distribution(y, y_pc2, 'PC2')
detect_outliers(y, y_pc2)
intercept_pc2, coefs_pc2 = ploy_coef(pc2_estimator, pc_reg)


#%%Try different model of regression
'''
PCR first order
'''

pc1_estimator  = linear_regression(1)
pc1_estimator.fit(Xreg, y)
y_pc1 = pc1_estimator.predict(Xreg)
RMSE_pc1, r2_pc1, scores_pc1 = regression_pipeline(Xreg, y, pc1_estimator, 'PC1')
intercept_pc1, coefs_pc1 = ploy_coef(pc1_estimator, pc_reg)



'''
1st order linear regression
'''

first  = linear_regression(1)
first.fit(X, y)
y_first = first.predict(X)
RMSE_first, r2_first, scores_first = regression_pipeline(X, y, first, 'First Order')
intercept_first, coefs_first = ploy_coef(first, X.shape[1])

'''
2nd order linear regression
'''
second  = linear_regression(2)
second.fit(X, y)
y_second = second.predict(X)
RMSE_second, r2_second, scores_second = regression_pipeline(X, y, second, 'Second Order')
intercept_second, coefs_second = ploy_coef(second,  X.shape[1])


#%%
'''
Compare different regression models 
'''

regression_method = ['PCR 1st Order', 'PCR 2nd Order', '1st Order', '2nd Order']
means_train = np.array([scores_pc1[0], scores_pc2[0], scores_first[0], scores_second[0]])
std_train = np.array([scores_pc1[1], scores_pc2[1], scores_first[1], scores_second[1]])
means_test = np.array([scores_pc1[2], scores_pc2[2], scores_first[2], scores_second[2]])
std_test = np.array([scores_pc1[3], scores_pc2[3], scores_first[3], scores_second[3]])
base_line = 0
x_pos = np.arange(len(regression_method))
opacity = 0.8
bar_width = 0.35
plt.figure(figsize=(8,6))
rects1 = plt.bar(x_pos, means_train - base_line, bar_width, yerr=std_train,
                alpha=opacity, color='lightblue',
                label='Train')
rects2 = plt.bar(x_pos+bar_width, means_test - base_line, bar_width, yerr=std_test,  
                alpha=opacity, color='salmon',
                label='Test')
#plt.ylim([-1,18])
plt.xticks(x_pos+bar_width/2, regression_method, rotation=0)
plt.xlabel('Regression Method')
plt.ylabel('RMSE (eV)')
plt.legend(loc= 'upper right', frameon=False)


#%%PLS regression
from sklearn.cross_decomposition import PLSRegression

N = pc_reg
PLS = PLSRegression(n_components = N, tol=1e-8) #<- N_components tells the model how many sub-components to select
PLS.fit(X,y) #<- we have to pass y into the fit function now
yhat_PLS = PLS.predict(X)[:,0] #<- the prediction here is a column vector
# make a parity plot
mse_PLS, score_PLS = parity_plot_st(y, yhat_PLS, 'PLS')
sigma_PLS = error_distribution(y, yhat_PLS, 'PLS')

