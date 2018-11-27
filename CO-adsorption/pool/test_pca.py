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
import pickle
import numpy as np
import matplotlib as mat

mat.rc('font', size=12)
[dem, Eads, labels] = pickle.load(open("pca_data.p", "rb"))
X = dem
y = Eads

#%% Plot the trend for each discriptor

nPC = X.shape[1]

with plt.style.context('seaborn-whitegrid'):
    
    for cnt in range(nPC):
        plt.figure(figsize=(6, 4))
        plt.scatter(X[:,cnt],y)
        plt.xlabel(labels[cnt])
        plt.ylabel('CO Adsorption Energy (eV)')
    plt.legend(loc='upper right', fancybox=True, fontsize=8)
    plt.tight_layout()
    plt.show()

#%% PCA 
# Normalize the data
X_std = StandardScaler().fit_transform(X)
# Covriance matrix 
cov_mat = np.cov(X_std.T) 

# Plot Covariance structure 
plt.figure()
ax1 = plt.subplot(121)
ax1.set_title('X')
ax1.imshow(cov_mat)

# PCA use sklearn
pca = PCA()    
Xpc = pca.fit_transform(X_std) 
# Covriance matrix from PCA
cov_pc = np.cov(Xpc.T) 
# Plot Covariance structure 
ax2 = plt.subplot(122)
ax2.set_title('X_PCA')
ax2.imshow(cov_pc)

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
# select the number of PCs to plot in the bar graph
pc = 5
ind = 0
yvals = []
ylabels = []
bar_vals = []
space = 0.3

descriptors = labels
cm = ['r', 'coral', 'pink',  'orange', 'y', 'gold', 'lightblue', 'lime', 'grey'][:len(descriptors)]
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
           bbox_to_anchor = (1.02, 1),loc= 'upper left', prop={'size':10},frameon=False)

plt.plot(linex, linex*0, c = 'k', lw = 0.8)
plt.show()

#%% Regression
# Create linear regression object
pc_reg = 7 #len(descriptors)
Xreg = Xpc[:,:pc_reg]
def fit_linear_regression(X, y, degree):
    return Pipeline([("polynomial_features", PolynomialFeatures(degree=degree,
                                                                include_bias=False)),
                     ("linear_regression", linear_model.LinearRegression())]
                    ).fit(X, y) 
        
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
regr_poly = estimator.named_steps['linear_regression']
coefs = regr_poly.coef_
poly = estimator.named_steps['polynomial_features']
terms = poly.get_feature_names(['x1','x2','x3','x4','x5','x6','x7'])

y_poly = estimator.predict(Xreg)
score_poly = r2_score(y, y_poly)
mse_poly = mean_squared_error(y, y_poly)

fig, ax = plt.subplots()
ax.scatter(y, y_poly, facecolor = 'r', s  = 60, edgecolors=(0, 0, 0))
ax.plot([y.min(), y.max()], [y.min(), y.max()], 'k--', lw=2)
ax.set_xlabel('Measured')
ax.set_ylabel('Predicted')
plt.show()

xi = np.arange(len(coefs))
fig, ax = plt.subplots()
plt.bar(xi, coefs)
linex = np.arange(xi.min()-1, xi.max()+2)
plt.plot(linex, linex*0, c = 'k')
plt.xticks(xi, terms, rotation=45 )
plt.ylabel("Regression Coefficient Value (eV)")
plt.xlabel("Regression Coefficient")  
plt.show()

#%%PLS regression
from sklearn.cross_decomposition import PLSRegression

N = 7
PLS = PLSRegression(n_components = N, tol=1e-8) #<- N_components tells the model how many sub-components to select
PLS.fit(X,y) #<- we have to pass y into the fit function now
yhat_PLS = PLS.predict(X)[:,0] #<- the prediction here is a column vector

fig, ax = plt.subplots()

# make a parity plot
ax.scatter(y,yhat_PLS)
ax.plot(y,y,ls='--',color='k')

# print R^2
print(PLS.score(X,y))

# print MAE
print(np.mean(np.abs(y-yhat_PLS)))
