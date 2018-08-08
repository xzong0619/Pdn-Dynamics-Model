# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 11:35:39 2018

@author: wangyf
"""
from config import mother
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py
import plotly.graph_objs as go

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
color = []

for i in range(len(mother)):
    if mother[i,2] == 1: color.append('r')
    if mother[i,2] == 2: color.append('b')
    if mother[i,2] == 3: color.append('g')
    if mother[i,2] == 4: color.append('y')
x = mother[:,0]
y = mother[:,1]
z = mother[:,2]

for i in range(len(mother)): #plot each point + it's index as text above
    ax.text(x[i],y[i],z[i],  '%s' % (str(i)), size=10, color='k') 
    ax.scatter(x[i],y[i],z[i],c=color[i], marker='o', s=10) 
    
 
rd = 90

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.view_init(azim=rd)

plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c=color, marker='o', s= 500)
ax.view_init(azim=rd)



    