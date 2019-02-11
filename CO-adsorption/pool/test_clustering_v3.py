# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 21:11:52 2019

@author: yifan
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist

from sklearn.datasets.samples_generator import make_blobs

def plot_kmeans(kmeans, X, n_clusters=4, rseed=0, ax=None):
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
                                
                                
X, y_true = make_blobs(n_samples=400, centers=4,
                       cluster_std=0.60, random_state=0)
X = X[:, ::-1] # flip axes for better plotting

#%% K means
kmeans = KMeans(4, random_state=0)
labels = kmeans.fit(X).predict(X)
plt.figure()
plot_kmeans(kmeans, X)

rng = np.random.RandomState(13)
X_stretched = np.dot(X, rng.randn(2, 2))
plt.figure()
plot_kmeans(kmeans, X_stretched)

#%% Gaussian Mixture models
from sklearn.mixture import GaussianMixture
gmm = GaussianMixture(n_components=4).fit(X)
labels = gmm.predict(X)
plt.figure()
ax = plt.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis')

from matplotlib.patches import Ellipse

def draw_ellipse(position, covariance, ax=None, option = 'a', **kwargs):
    """Draw an ellipse with a given position and covariance"""
    ax = ax or plt.gca()
    if option == 'a':
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
    if option == 'b':    
        for n, color in enumerate(colors):
            if gmm.covariance_type == 'full':
                covariances = gmm.covariances_[n][:2, :2]
            elif gmm.covariance_type == 'tied':
                covariances = gmm.covariances_[:2, :2]
            elif gmm.covariance_type == 'diag':
                covariances = np.diag(gmm.covariances_[n][:2])
            elif gmm.covariance_type == 'spherical':
                covariances = np.eye(gmm.means_.shape[1]) * gmm.covariances_[n]
            v, w = np.linalg.eigh(covariances)
            u = w[0] / np.linalg.norm(w[0])
            angle = np.arctan2(u[1], u[0])
            angle = 180 * angle / np.pi  # convert to degrees
            v = 2. *  np.sqrt(2.) * np.sqrt(v)
            ax.add_patch(Ellipse(gmm.means_[n, :2], v[0], v[1],
                                      180 + angle, color=color))
        
def plot_gmm(gmm, X, option = 'a',label=True, ax=None):
    ax = ax or plt.gca()
    labels = gmm.fit(X).predict(X)
    if label:
        ax.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis', zorder=2)
    else:
        ax.scatter(X[:, 0], X[:, 1], s=40, zorder=2)
    ax.axis('equal')
    
    w_factor = 0.2 / gmm.weights_.max()
    for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
        draw_ellipse(pos, covar, alpha=w * w_factor, option = option)
        
colors = ['navy', 'turquoise', 'darkorange']




gmm = GaussianMixture(n_components=4, random_state=42)
plt.figure()
plot_gmm(gmm, X, option = 'a')


plt.figure()
plot_gmm(gmm, X, option = 'b')


