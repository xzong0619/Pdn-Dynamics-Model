# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 12:00:34 2018

@author: wangyf
"""
 
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler 
from sklearn import linear_model 
from sklearn.model_selection import cross_val_predict 
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import PolynomialFeatures 
from sklearn.pipeline import Pipeline
from scipy.stats import norm
import pickle
import numpy as np
import matplotlib as mat

mat.rc('font', size=12)
[dem, Eads, descriptors, filename_list, sitetype_list] = pickle.load(open("pca_data.p", "rb"))
X = dem
y = Eads

#%% PCA parameters

nPC = X.shape[1]
# select the number of PCs to plot in the bar graph
pc = len(descriptors)
# select the number of PCs to plot in regression
pc_reg = min(7, pc) 


#%% Plot the trend for each discriptor
def plot_discriptors():
    '''
    Plot the trend for each discriptor
    '''
    with plt.style.context('seaborn-whitegrid'):
        
        for cnt in range(nPC):
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
        for cnt in range(nPC):
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

plot_discriptors_st()

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

    plt.bar(range(nPC), var_exp, alpha=0.5, align='center',
            label='individual explained variance')
    plt.step(range(nPC), cum_var_exp, where='mid',
             label='cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal components')
    plt.xticks(np.arange(nPC), 
                ['PC %i'%(w+1) for w in range(nPC)])
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
def parity_plot(yobj,ypred, method = 'PCA'):
    '''
    Plot the parity plot of y vs ypred
    return R2 score and MSE for the model
    '''
    MSE = mean_squared_error(yobj, ypred)
    score = r2_score(yobj, ypred)
    fig, ax = plt.subplots()
    ax.scatter(yobj, ypred, facecolor = 'r', s  = 60, edgecolors=(0, 0, 0))
    ax.plot([yobj.min(), yobj.max()], [yobj.min(), yobj.max()], 'k--', lw=2)
    ax.set_xlabel('Objective')
    ax.set_ylabel('Predicted')
    plt.title(r'Method-{}, MSE-{:.2}, $r^2$ -{:.2}'.format(method, MSE, score))
    plt.show()
    
    return MSE, score

def parity_plot_st(yobj,ypred, method = 'PCA'):
    '''
    Plot the parity plot of y vs ypred
    return R2 score and MSE for the model
    colorcode different site types
    '''
    MSE = mean_squared_error(yobj, ypred)
    score = r2_score(yobj, ypred)
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
    ax.set_xlabel('Objective')
    ax.set_ylabel('Predicted')
    plt.title(r'Method-{}, MSE-{:.2}, $r^2$ -{:.2}'.format(method, MSE, score))
    plt.legend(bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)
    plt.show()
    
    return MSE, score

def error_distribution(yobj, ypred, method  = 'PCA'):
    
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
    plt.title(r'Method-{}, $\sigma$-{:.2}'.format(method, sigma))
    
    return sigma

def fit_linear_regression(X, y, degree):
    '''
    # Create linear regression object
    '''
    return Pipeline([("polynomial_features", PolynomialFeatures(degree=degree,
                                                                include_bias=False)),
                     ("linear_regression", linear_model.LinearRegression())]
                    ).fit(X, y)  
   

Xreg = Xpc[:,:pc_reg]
degree = 2

#regr = linear_model.LinearRegression()
## Fit
#regr.fit(Xreg, y)
#
## Calibration
#y_c = regr.predict(Xreg)
#
## Cross-validation
#y_cv = cross_val_predict(regr, Xreg, y, cv=20)
#
## Calculate scores for calibration and cross-validation
#score_c = r2_score(y, y_c)
#score_cv = r2_score(y, y_cv)
#
## Calculate mean square error for calibration and cross validation
#mse_c = mean_squared_error(y, y_c)
#mse_cv = mean_squared_error(y, y_cv)

estimator  = fit_linear_regression(Xreg, y, degree)
regr_pca = estimator.named_steps['linear_regression']
coefs = regr_pca.coef_
poly = estimator.named_steps['polynomial_features']
feature_names = []
for i in range(pc_reg): feature_names.append('x'+ str(i+1))
terms = poly.get_feature_names(feature_names)

y_pca = estimator.predict(Xreg)
mse_pca, score_pca = parity_plot_st(y, y_pca)
sigma_pca = error_distribution(y, y_pca)

'''
Plot the magnitude of each parameters
'''
xi = np.arange(len(coefs))
fig, ax = plt.subplots()
plt.bar(xi, coefs)
linex = np.arange(xi.min()-1, xi.max()+2)
plt.plot(linex, linex*0, c = 'k')
plt.xticks(xi, terms, rotation=45, fontsize = 8 )
plt.ylabel("Regression Coefficient Value (eV)")
plt.xlabel("Regression Coefficient")  
plt.show()

#%%detect the outlier
threshold = 0.2
diff = abs(y-y_pca)
outlier_index = np.where(diff>threshold)[0]
outlier_file = []
outlier_stype = []
for i in outlier_index: 
    outlier_file.append(filename_list[i])
    outlier_stype.append(sitetype_list[i])

#%%PLS regression
from sklearn.cross_decomposition import PLSRegression

N = pc_reg
PLS = PLSRegression(n_components = N, tol=1e-8) #<- N_components tells the model how many sub-components to select
PLS.fit(X,y) #<- we have to pass y into the fit function now
yhat_PLS = PLS.predict(X)[:,0] #<- the prediction here is a column vector
# make a parity plot
mse_PLS, score_PLS = parity_plot_st(y, yhat_PLS, 'PLS')
sigma_PLS = error_distribution(y, yhat_PLS, 'PLS')


