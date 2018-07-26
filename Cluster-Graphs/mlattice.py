# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 11:08:30 2018

@author: wangyf
"""
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import math

def two_points_D(A,B):
    '''
    Calculates the distance between two points A and B
    '''
    d = math.sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)
    return d


def drawing(G):
    '''
    takes in graph g, draw it
    '''
    color = nx.get_node_attributes(G,'color')
    pos=nx.get_node_attributes(G,'pos')
    plt.figure() 
    nx.draw(G, pos, with_labels=True, node_color = list(color.values()))
    
def gmothers(mother):
    
    nm = len(mother)
    Gm = nx.Graph()
    
    for i in range(nm):
        Gm.add_node(i, pos = mother[i], color = empty)
    
    edge_d = []
    edge = []
    
    # Add all egdes and calculate the edge distance
    for i in range(nm):
        for j in np.arange(i+1,nm):
            edge.append((i,j))
            edge_d.append(two_points_D(mother[i],mother[j]))
            
    ne = len(edge)
    for i in range(ne):
        #if edge_d[i] <= 1:
        Gm.add_edges_from([edge[i]], length = edge_d[i])
    drawing(Gm)
    plt.title('%d lattice points' %nm)
    
    return Gm
    
def gconfigurations(mother, son):
    '''
    takes in mother coordinate list and son's index number
    returns the shaded son graph
    '''
    Gs = nx.Graph()
    nm = len(mother)
    ns = len(son)
    
    
    edge_d = []
    edge = []
    for i in range(nm):
        for j in np.arange(i+1,nm):
            edge.append((i,j))
            edge_d.append(two_points_D(mother[i],mother[j]))
            
    ne = len(edge)
    for i in range(ne):
        Gs.add_edges_from([edge[i]], length = edge_d[i])
    for i in range(nm):
        Gs.add_node(i, pos = mother[i], color = empty)
    for i in range(ne):
        if edge_d[i] <= 1:
            Gs.add_edges_from([edge[i]], length = edge_d[i])    
    for si in range(ns):
        Gs.node[son[si]]['color'] = filled
    
    drawing(Gs)
    plt.title('Pd %d' %ns)

    return Gs

def gclusters(cmother, son):
    ns = len(son)
    mother = []
    for i in range(ns):
        mother.append(cmother[i])
    Gc = gmothers(mother)
    return Gc
    
'''
main
'''
    

mother = [(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5),
          (1/2, -3**0.5/2), (3/2, -3**0.5/2), (0, -3**0.5),
          (-1,0), (-1/2, -3**0.5/2), (-3/2, -3**0.5/2),
          (-1/2, 3**0.5/2), (-3/2, 3**0.5/2)]
config = [[0],
          [0,1],
          [0,1,2],
          [0,1,2,3],
          [0,1,2,3,5],
          [0,1,2,3,5,6],
          [0,1,2,3,5,9],
          [0,1,2,3,5,11],
          [0,1,2,5,9,8,11],
          [0,1,2,3,5,9,7],
          [0,1,2,5,8,11,12],
          [0,1,2,8,9,11,12]]
         
empty = 'grey'
filled = 'r'

Gm = gmothers(mother)
    
        
'''
Creat 12 configurations
'''
ns = len(config)
Gsv = []
for si in range(ns):
    son = config[si]
    Gs = gconfigurations(mother,son)
    Gsv.append(Gs)



'''
Creat 6 clusters
'''

cmother = [(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5)]
cconfig = [[0], [0,1], [0,3], [1,4],[0,1,2],[1,2,4], [0,1,3]]
cedge = 


cm = gmothers(cmother)
cns = len(cconfig)
Gcv = []
for si in range(cns):
    cson = cconfig[si]
    Gc = gclusters(cmother,cson)
    Gcv.append(Gc)       
        
        