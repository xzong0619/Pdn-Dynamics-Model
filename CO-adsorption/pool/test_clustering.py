# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 12:09:22 2019

@author: wangyf
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


def dist(pt1, pt2):
    "Euclidean distance between two points"
    #note that this can also be performed with np.linalg.norm(pt1-pt2)
    return np.sqrt(sum([(xi-yi)**2 for xi, yi in zip(pt1, pt2)]))

def expected_assignment(pt, cluster_centers):
    dists = [dist(pt,ci) for ci in cluster_centers] #<- find distance to each center
    min_index = dists.index(min(dists)) #<- find the index (cluster) with the minimum dist
    return min_index

def new_centers(cluster_points, centers):
    centers = list(centers)
    for i,ci in enumerate(cluster_points):
        if ci != []:
            centers[i] = np.mean(ci, axis=0)
    return centers


X = X_blobs
y = y_blobs

#cluster_centers = ([-0.5,0], [0.5,0])
cluster_centers = ([-0.5,0], [-4,3])
fig, ax = plt.subplots()
ax.scatter(X[:,0], X[:,1], color='k', alpha=0.2)
colors = {0:'r', 1:'g', 2:'b',3:'m'}
for i,ci in enumerate(cluster_centers):
    ax.plot(ci[0], ci[1], marker='*', markersize='12', color=colors[i])
    
    
    
fig, ax = plt.subplots()
# Plot old centers
for i,ci in enumerate(cluster_centers):
    ax.plot(ci[0], ci[1], marker='+', markersize='12', color=colors[i])
    
# Which cluster do we "expect" each point to belong to?
clusters = [[],[],[],[]]
for pt in X:
    cluster_idx = expected_assignment(pt, cluster_centers)
    clusters[cluster_idx].append(pt)
    
# What centers best represent these new assignments?
cluster_centers = new_centers(clusters, cluster_centers)

# Plot new assignments
for i, ci in enumerate(clusters):
    for pt in ci:
        ax.plot(pt[0], pt[1], marker='o', color=colors[i], alpha=0.2)
        
# Plot new centers
for i,ci in enumerate(cluster_centers):
    print(ci)
    ax.plot(ci[0], ci[1], marker='*', markersize='12', color=colors[i])
    
#%%
from sklearn.cluster import KMeans

model = KMeans(n_clusters=2, random_state=random_state)
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
    