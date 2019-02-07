#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 15:26:16 2019

@author: wangyifan
"""

import pylab as plt
import numpy as np
from sklearn.datasets import make_blobs, make_moons, make_biclusters, make_circles
from sklearn.metrics import silhouette_samples, silhouette_score, calinski_harabaz_score

random_state = 1
k = 4
N = 300

X_blobs, y_blobs = make_blobs(N,centers=k,random_state=random_state)
transformation = [[0.6, -0.6], [-0.4, 0.8]]
X_aniso = np.dot(X_blobs, transformation)

X_circles, y_circles = make_circles(N,noise=0.1, factor=0.4, random_state=random_state)

X_moons, y_moons = make_moons(N, noise=0.1, random_state=random_state)

fig, axes = plt.subplots(1,4, figsize=(15,4))

axes[0].scatter(X_blobs[:,0], X_blobs[:,1], c=y_blobs)
axes[1].scatter(X_aniso[:,0], X_aniso[:,1], c=y_blobs)
axes[2].scatter(X_circles[:,0], X_circles[:,1], c=y_circles)
axes[3].scatter(X_moons[:,0], X_moons[:,1], c=y_moons)


X = X_blobs
y = y_blobs

#cluster_centers = ([-0.5,0], [0.5,0])
cluster_centers = ([-0.5,0], [-4,3])
fig, ax = plt.subplots()
ax.scatter(X[:,0], X[:,1], color='k', alpha=0.2)
colors = {0:'r', 1:'g', 2:'b',3:'m'}
for i,ci in enumerate(cluster_centers):
    ax.plot(ci[0], ci[1], marker='*', markersize='12', color=colors[i])
    
    
X_1d = X[:,0]

fig, ax = plt.subplots()
ax.hist(X_1d, bins=10)

from sklearn.cluster import KMeans

model = KMeans(n_clusters=4, random_state=random_state)
model.fit(X)
y_predict = model.predict(X)
centers = model.cluster_centers_

fig, axes = plt.subplots(1,2,figsize=(8,4))
axes[0].scatter(X[:,0], X[:,1], c=y_predict, cmap='RdBu')
axes[1].scatter(X[:,0], X[:,1], c=y)
for center in centers:
    x_i = center[0]
    y_i = center[1]
    axes[0].plot(x_i, y_i, marker='*', color='k', mec='w', markersize=20)
    

k = 25
model = KMeans(n_clusters=k, random_state=random_state)
X_1d_skl = X_1d.reshape(-1,1)
model.fit(X_1d_skl)
y_predict = model.predict(X_1d_skl)
centers = model.cluster_centers_
unique, counts = np.unique(y_predict, return_counts=True) #<- useful way to count occurences!

fig, ax = plt.subplots()
ax.hist(X_1d, bins=k, alpha=0.5)
ax.scatter(centers,counts)