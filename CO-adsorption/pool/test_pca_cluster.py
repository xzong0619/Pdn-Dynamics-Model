# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 07:32:06 2019

@author: wangyf
"""

'''
Test PCA data with clustering technqiue 
'''


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
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt 
import matplotlib


matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['ytick.major.width'] = 2



#%% Processs descriptor_data.cvs
fdata = pd.read_csv('descriptor_data.csv')

#possible descriptors
descriptors =  ['NPd', 'CN1', 'CN2','GCN', 'Z', 'Charge', 'Nsites', 'Pd1C', 'Pd2C', 'Pd3C'] #10 in total
descriptors_g =  ['CN1', 'Z', 'Nsites',   'Pd1C', 'Pd2C', 'Pd3C'] #6 geometric descriptors
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


#%% Cluster section
'''
cluster plot functions
'''
def plot_kmeans(kmeans, X,  ax=None):
    labels = kmeans.fit_predict(X)

    # plot the input data
    ax = ax or plt.gca()
    ax.axis('equal')
    ax.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis', zorder=2)

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
    for nsig in range(1, 3):
        ax.add_patch(Ellipse(position, nsig * width, nsig * height,
                             angle, **kwargs))

        
def plot_gmm(gmm, X, label=True, ax=None):
    ax = ax or plt.gca()
    labels = gmm.fit(X).predict(X)
    if label:
        ax.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis', zorder=2)
    else:
        ax.scatter(X[:, 0], X[:, 1], s=40, zorder=2)
    ax.axis('equal')
    
    w_factor = 0.2 / gmm.weights_.max()
    for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
        draw_ellipse(pos, covar, alpha=w * w_factor)

'''
Plot important PC vs each other and observe clustering
'''
pc1 = Xpc[:,0]
pc2 = Xpc[:,1]
pc3 = Xpc[:,2]


# pc1 vs pc2
plt.figure(figsize=(6, 4))
for site, col in zip(('top', 'bridge', 'hollow'),
            ('red', 'green', 'blue')):
    indices = np.where(np.array(sitetype_list) == site)[0]
    plt.scatter(pc1[indices],
                pc2[indices],
                label=site,
                facecolor = col, 
                alpha = 0.5,
                s  = 60)
    
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend(bbox_to_anchor = (1.01, 1),loc= 'upper left', frameon=False)
plt.tight_layout()  
plt.show() 


plt.figure(figsize=(6,4))
plt.scatter(pc1, pc2, c = y, s = 80)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.colorbar(shrink = 1.8)


random_state =2
X_cluster = Xpc[:,0:2]
kmeans = KMeans(n_clusters = 3, random_state=random_state)
kmeans.fit(X_cluster)
y_predict = kmeans.predict(X_cluster)
plt.figure(figsize=(5,4))
plot_kmeans(kmeans, X_cluster)
plt.xlabel('PC1')
plt.ylabel('PC2')

gmm = GaussianMixture(n_components=3, random_state=random_state)
plt.figure(figsize=(5,4))
plot_gmm(gmm, X_cluster)
plt.xlabel('PC1')
plt.ylabel('PC2')


# Add more points
import random

random.seed(random_state)
def get_random_num(n, r):
    rset = []
    for i in range(n):
        rset.append(random.random() * (r[1] - r[0]) + r[0])    
    rset = np.array(rset)
    
    return rset
    

range1 = [-4, 0]
range2  = [-1, 2]

rset1 =get_random_num(10, range1)
rset2 =get_random_num(10, range2)
rset = np.transpose(np.vstack((rset1, rset2)))
X_cluster_new = np.vstack((X_cluster, rset))

plt.figure(figsize=(5,4))
plot_kmeans(kmeans, X_cluster_new)
plt.xlabel('PC1')
plt.ylabel('PC2')

plt.figure(figsize=(5,4))
plot_gmm(gmm, X_cluster_new)
plt.xlabel('PC1')
plt.ylabel('PC2')















# pc1 vs pc3
plt.figure(figsize=(6, 4))
for site, col in zip(('top', 'bridge', 'hollow'),
            ('red', 'green', 'blue')):
    indices = np.where(np.array(sitetype_list) == site)[0]
    plt.scatter(pc1[indices],
                pc3[indices],
                label=site,
                facecolor = col, 
                alpha = 0.5,
                s  = 60)
    
plt.xlabel('PC1')
plt.ylabel('PC3')
plt.legend(bbox_to_anchor = (1.01, 1),loc= 'upper left', frameon=False)
plt.tight_layout()  
plt.show() 


plt.figure(figsize=(6,4))
plt.scatter(pc1, pc3, c = y, s = 80)
plt.xlabel('PC1')
plt.ylabel('PC3')
plt.colorbar(shrink = 1.8)


random_state =2
X_cluster = Xpc[:,[0,2]]
kmeans = KMeans(n_clusters = 3, random_state=random_state)
kmeans.fit(X_cluster)
y_predict = kmeans.predict(X_cluster)
plt.figure(figsize=(5,4))
plot_kmeans(kmeans, X_cluster)
plt.xlabel('PC1')
plt.ylabel('PC3')

gmm = GaussianMixture(n_components= 3, random_state=random_state)
plt.figure(figsize=(5,4))
plot_gmm(gmm, X_cluster)
plt.xlabel('PC1')
plt.ylabel('PC3')
# Cluster for parity plot
# Cluster for PC1-PC2
# Cluster for parameters



