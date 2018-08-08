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
    
    
ax.scatter(mother[:,0], mother[:,1], mother[:,2],c= color, marker='o', s=500)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.view_init(azim=60)

for i, txt in enumerate(n):
    ax.annotate(txt, (z[i], y[i]))
    
plt.show()

#%%
trace1 = go.Scatter3d(
    x=mother[:,0],
    y=mother[:,1],
    z=mother[:,2],
    mode='markers',
    marker=dict(
        size=20,
        line=dict(
            color=color,
            width=0.5
        ),
        opacity=0.8
    )
)
        
layout = go.Layout(
    margin=dict(
        l=0,
        r=0,
        b=0,
        t=0
    )
)        
fig = go.Figure(data=[trace1], layout=layout)
py.iplot(fig, filename='simple-3d-scatter')