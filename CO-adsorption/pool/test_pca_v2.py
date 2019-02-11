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
from sklearn.model_selection import RepeatedKFold, cross_validate, LeaveOneOut
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from sklearn.mixture import GaussianMixture

import pandas as pd
import numpy as np
import pickle

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt 
import matplotlib
from matplotlib.patches import Ellipse


matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['ytick.major.width'] = 2



#%% Processs descriptor_data.cvs
fdata = pd.read_csv('descriptor_data.csv')

#possible descriptors
descriptors =  ['NPd', 'CN1', 'CN2','GCN', 'Z', 'Charge', 'Nsites', 'Pd1C', 'Pd2C', 'Pd3C', 'CeNN1', 'ONN1'] #10 in total
descriptors_g =  ['NPd', 'CN1', 'Z', 'Nsites',  'CeNN1', 'ONN1'] #6 geometric descriptors
#descriptors = ['NPd', 'CN1', 'Z', 'Charge',  'Pd1C', 'Pd2C', 'Pd3C'] 
#descriptors =  ['CN1', 'Z', 'Charge',  'Pd1C', 'Pd2C', 'Nsites']
#descriptors =  ['CN1', 'Z', 'Charge',  'Pd1C', 'Pd2C', 'Pd3C']
#descriptors =  ['GCN', 'Z', 'Charge',  'Pd1C', 'Pd2C', 'Pd3C']
dem1 =  np.array(fdata.loc[:,descriptors], dtype = float)
Eads = np.array(fdata.loc[:,'Eads'], dtype = float)
filename_list = list(fdata.loc[:,'Filename'])
sitetype_list = list(fdata.loc[:,'SiteType'])

#%%Count the number of sites 
ntop = (np.array(sitetype_list) == 'top').astype(int).sum()
nbridge = (np.array(sitetype_list) == 'bridge').astype(int).sum()
nhollow = (np.array(sitetype_list) == 'hollow').astype(int).sum()

# show one example -Pd10 hollow
fdata.loc[1]

X = dem1
y = Eads

#%% PCA parameters

nDescriptors = X.shape[1]
# select the number of PCs to plot in the bar graph
pc_draw = min(5, len(descriptors))
# select the number of PCs to plot in regression
pc_reg = min(7, len(descriptors)) 


#%% Plot the trend for each discriptor
def plot_discriptors():
    '''
    Plot the trend for each discriptor
    '''
        
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
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(15,6))
ax1.set_title('X')
ax1.set_xticks(range(0,nDescriptors))
ax1.set_xticklabels(descriptors, rotation = 40)
ax1.set_yticks(range(0,nDescriptors))
ax1.set_yticklabels(descriptors)
im1 = ax1.imshow(cov_mat)

fig.colorbar(im1, ax = ax1, shrink = 0.5)
ax2.set_title('X_PCA')
ax2.set_xticks(range(0,nDescriptors))
ax2.set_xticklabels(descriptors, rotation = 40)
ax2.set_yticks(range(0,nDescriptors))
ax2.set_yticklabels(descriptors)
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

'''
Scree plot for PCs
'''
                      
plt.figure(figsize=(6, 4))

plt.bar(range(nDescriptors), var_exp, alpha=0.5, align='center',
        label='individual explained variance')
plt.step(range(nDescriptors), cum_var_exp, where='mid',
         label='cumulative explained variance')
plt.ylabel('Explained variance ratio')
plt.xlabel('Principal components')
plt.xticks(np.arange(nDescriptors), 
            ['PC%i'%(w+1) for w in range(nDescriptors)])
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
space = 0.4

cm = ['r', 'coral', 'pink',  'orange',  'gold', 'y','lightgreen', 'lightblue',  'c', 'mediumpurple', 'brown', 'grey', 'orchid'][:len(descriptors)]
fig = plt.figure(figsize=(12,6))


ax = fig.add_subplot(111)
n = len(descriptors)
width = (1 - space) / (len(descriptors))
indeces = np.arange(0, pc_draw) + 0.5  

# Create a set of bars at each position
for i, pci in enumerate(eig_vecs[:pc_draw]):
    
    vals = pci/np.sum(np.absolute(pci))
    
    pos = width*np.arange(n) + i 
    ax.bar(pos, vals, width=width, label=str(i+1), color = cm, alpha = 1) 
        
linex = np.arange(np.arange(0, pc_draw).min() -0.1  , np.arange(0, pc_draw).max()+ 1)

ax.set_xticks(indeces)
ax.set_xticklabels(list(np.arange(0,pc_draw)+1))
ax.set_ylabel("Normalized Descriptoor Loading")
ax.set_xlabel("Principal Component #")    
  
# Add legend using color patches
patches = []
for c in range(n):
    patches.append(mpatches.Patch(color=cm[c]))
plt.legend(patches, descriptors, 
           bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)

plt.plot(linex, linex*0, c = 'k', lw = 1.5)
plt.show()

#%% Regression
'''
cluster plot functions
'''
def plot_kmeans(kmeans, X,  ax=None):
    labels = kmeans.fit_predict(X)

    # plot the input data
    ax = ax or plt.gca()
    ax.axis('equal')
    #ax.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis', zorder=2)

    # plot the representation of the KMeans model
    centers = kmeans.cluster_centers_
    radii = [cdist(X[labels == i], [center]).max()
             for i, center in enumerate(centers)]
    for c, r in zip(centers, radii):
        ax.add_patch(plt.Circle(c, r, fc='#CCCCCC', lw=3, alpha=0.5, zorder=1))
        ax.plot(c[0], c[1], marker='*', color='k', mec='w', markersize=20)

        
def draw_ellipse(position, covariance, ax=None, **kwargs):
    """Draw an ellipse with a given position and covariance"""
    ax = ax or plt.gca()

    #Convert covariance to principal axes
    if covariance.shape == (2, 2):
        U, s, Vt = np.linalg.svd(covariance)
        angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
        width, height = 2 * np.sqrt(s)
    else:
        angle = 0
        width, height = 2 * np.sqrt(covariance)
    
    # Draw the Ellipse
    for nsig in range(1, 4):
        ax.add_patch(Ellipse(position, nsig * width, nsig * height,
                             angle, **kwargs))

        
def plot_gmm(gmm, X, label=True, ax=None):
    ax = ax or plt.gca()
    labels = gmm.fit(X).predict(X)
#    if label:
#        ax.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis', zorder=2)
#    else:
#        ax.scatter(X[:, 0], X[:, 1], s=40, zorder=2)
#    ax.axis('equal')
#    
    w_factor = 0.15 / gmm.weights_.max()
    for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
        draw_ellipse(pos, covar, alpha=w * w_factor)



def parity_plot_st(yobj, ypred, method):
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
    
    X_cluster = np.transpose(np.vstack((yobj, ypred)))
    gmm = GaussianMixture(n_components= 1, random_state = 0)
    plot_gmm(gmm, X_cluster)
    
#    kmeans = KMeans(n_clusters = 1, random_state=0)
#    plot_kmeans(kmeans, X_cluster)


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
    #rkf = RepeatedKFold(n_splits = 3, n_repeats = 10) #change this to leave one out
    loo = LeaveOneOut()
    scores  = cross_validate(estimator, X, y, cv = loo,
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

def cross_validation_r2(X, y, estimator): 
    '''
    Cross-validation
    '''
    #rkf = RepeatedKFold(n_splits = 3, n_repeats = 10) #change this to leave one out
    loo = LeaveOneOut()
    scores  = cross_validate(estimator, X, y, cv = loo,
                                scoring=('r2'),
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
    '''
    Standard regression procedure including predict, plot parity plot and 
    perform cross validation
    '''

    y_predict = estimator.predict(X)
    if len(y_predict.shape) > 1: 
        y_predict = y_predict[:,0]
    RMSE, r2 = parity_plot_st(y, y_predict, method)
    scores = cross_validation(X, y, estimator)
    r2s = cross_validation_r2(X, y, estimator)
    
    return RMSE, r2, scores, r2s
    
def detect_outliers(y, y_predict, threshold = 0.2):
    '''
    detect the outlier
    '''
    diff = abs(y-y_predict)
    outlier_index = np.where(diff>threshold)[0]
    outliers_files = []
    outliers_sites = []
    print('\nOutliers with AE in eV:')
    for i in outlier_index: 
        outliers_files.append(filename_list[i])
        outliers_sites.append(sitetype_list[i])
        print('{} {} {} - {:.2}'.format(i,filename_list[i],sitetype_list[i], diff[i]))
        
    return outliers_files, outliers_sites

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

RMSE_pc2, r2_pc2, scores_pc2, r2s_pc2 = regression_pipeline(Xreg, y, pc2_estimator, 'PC2')
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
RMSE_pc1, r2_pc1, scores_pc1, r2s_pc1 = regression_pipeline(Xreg, y, pc1_estimator, 'PC1')
intercept_pc1, coefs_pc1 = ploy_coef(pc1_estimator, pc_reg)



'''
1st order linear regression
'''

first  = linear_regression(1)
first.fit(X, y)
y_first = first.predict(X)
RMSE_first, r2_first, scores_first, r2s_first = regression_pipeline(X, y, first, 'First Order')
intercept_first, coefs_first = ploy_coef(first, X.shape[1])

'''
2nd order linear regression
'''
second  = linear_regression(2)
second.fit(X, y)
y_second = second.predict(X)
RMSE_second, r2_second, scores_second, r2s_second = regression_pipeline(X, y, second, 'Second Order')
intercept_second, coefs_second = ploy_coef(second,  X.shape[1])

'''
PLS regression 
'''
from sklearn.cross_decomposition import PLSRegression

PLS = PLSRegression(n_components = pc_reg, tol=1e-8) #<- N_components tells the model how many sub-components to select
PLS.fit(X,y) 
y_PLS = PLS.predict(X)[:,0] #<- the prediction here is a column vector
RMSE_PLS, r2_PLS, scores_PLS, r2s_PLS = regression_pipeline(X, y, PLS, 'PLS')

#%%
'''
Use geomertric descriptors only
'''
Xg =  np.array(fdata.loc[:,descriptors_g], dtype = float)
X_std_g = StandardScaler().fit_transform(Xg)
Xpc_g = pca.fit_transform(X_std_g)

'''
1st order pcg
'''
pcg_estimator_first  = linear_regression(2)
pcg_estimator_first.fit(Xpc_g, y)
y_pcg_first = pcg_estimator_first.predict(Xpc_g)
RMSE_pcg_first, r2_pcg_first, scores_pcg_first, r2s_pcg_first = regression_pipeline(Xpc_g, y, pcg_estimator_first, 'PCg')
detect_outliers(y, y_pcg_first)
intercept_pcg_first, coefs_pcg_first = ploy_coef(pcg_estimator_first, Xpc_g.shape[1])
#pickle.dump(pcg_estimator_second, open('g_estimator_second','wb'))

'''
2nd order pcg
'''
pcg_estimator_second  = linear_regression(2)
pcg_estimator_second.fit(Xpc_g, y)
y_pcg_second = pcg_estimator_second.predict(Xpc_g)
RMSE_pcg_second, r2_pcg_second, scores_pcg_second, r2s_pcg_second = regression_pipeline(Xpc_g, y, pcg_estimator_second, 'PCg')
detect_outliers(y, y_pcg_second)
intercept_pcg_second, coefs_pcg_second = ploy_coef(pcg_estimator_second, Xpc_g.shape[1])
#pickle.dump(pcg_estimator_second, open('g_estimator_second','wb'))

#%%
'''
Compare different regression models 
'''

'''
Plot CV RMSE
'''
regression_method = ['PLS', 'PCR 1st', 'PCR 2nd', 'PC G 1st',  'PC G 2nd', 'Poly 1st', 'Poly 2nd']
scores_mx = np.array([scores_PLS, scores_pc1, scores_pc2, scores_pcg_first, scores_pcg_second, scores_first, scores_second])
means_train = np.array(scores_mx[:,0])
std_train = np.array(scores_mx[:,1])
means_test = np.array(scores_mx[:,2])
std_test = np.array(scores_mx[:,3])

r2s = np.array([r2_PLS, r2_pc1, r2_pc2, r2_pcg_first, r2_pcg_second, r2_first, r2_second])

base_line = 0
x_pos = np.arange(len(regression_method))
opacity = 0.9
bar_width = 0.25
fig, ax1 = plt.subplots(figsize=(8,6))
rects1 = plt.bar(x_pos, means_train - base_line, bar_width, #yerr=std_train,
                alpha = opacity, color='lightblue',
                label='Train')
rects2 = plt.bar(x_pos+bar_width, means_test - base_line, bar_width, #yerr=std_test,  
                alpha = opacity, color='salmon',
                label='Test')
rects3 = plt.bar(x_pos+bar_width*2, r2s - base_line, bar_width, #yerr=std_test,  
                alpha = opacity, color='lightgreen',
                label='r2')
#plt.ylim([-1,18])

plt.xticks(x_pos+bar_width/2, regression_method, rotation=40)
plt.xlabel('Regression Method')
plt.ylim([0,1])
plt.legend(loc= 'best', frameon=False)

ax1.set_ylabel('CV RMSE (eV)')
ax2 = ax1.twinx()
ax2.set_ylabel('r2')

#'''
#Plot r2
#'''
#r2s_mx = np.array([r2s_PLS, r2s_pc1, r2s_pc2, r2s_pcg_first, r2s_pcg_second, r2s_first, r2s_second])
#means_train = np.array(r2s_mx[:,0])
#std_train = np.array(r2s_mx[:,1])
#means_test = np.array(r2s_mx[:,2])
#std_test = np.array(r2s_mx[:,3])
#base_line = 0
#x_pos = np.arange(len(regression_method))
#opacity = 0.8
#bar_width = 0.35
#plt.figure(figsize=(8,6))
#rects1 = plt.bar(x_pos, means_train - base_line, bar_width, #yerr=std_train,
#                alpha=opacity, color='lightblue',
#                label='Train')
#rects2 = plt.bar(x_pos+bar_width, means_test - base_line, bar_width, #yerr=std_test,  
#                alpha=opacity, color='salmon',
#                label='Test')
##plt.ylim([-1,18])
#plt.xticks(x_pos+bar_width/2, regression_method, rotation=40)
#plt.xlabel('Regression Method')
#plt.ylabel('CV r^2')
#plt.ylim([0,1])
#plt.legend(loc= 'best', frameon=False)

