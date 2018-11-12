# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 11:55:36 2018

@author: wangyf
"""

import numpy as np



import lattice_functions as lf
from structure_constants import mother, dz


NN1_edges1 = [(0, 1), (0, 2), (0, 6), (0, 10), (0, 11), (0, 14), (0, 20), (0, 21), 
             (0, 23), (0, 27), (0, 28), (0, 29), (1, 2), (1, 3), (1, 6), (1, 7), 
             (1, 21), (1, 24), (1, 29), (2, 3), (2, 4), (2, 5), (2, 14), (2, 18), 
             (2, 19), (2, 21), (2, 27), (3, 5), (3, 19), (4, 5), (4, 13), (4, 14), 
             (4, 18), (5, 19), (6, 7), (6, 8), (6, 9), (6, 11), (6, 23), (6, 24), 
             (6, 26), (6, 29), (7, 9), (7, 24), (8, 9), (8, 11), (8, 16), (8, 25), 
             (8, 26), (9, 26), (10, 11), (10, 12), (10, 14), (10, 15), (10, 20), 
             (10, 22), (10, 28), (11, 12), (11, 16), (11, 22), (11, 23), (11, 25), 
             (11, 28), (12, 16), (12, 22), (13, 14), (13, 15), (13, 17), (14, 15), 
             (14, 17), (14, 18), (14, 20), (14, 27), (15, 17), (16, 25), (17, 18), 
             (17, 20), (18, 19), (18, 20), (18, 21), (18, 27), (18, 30), (19, 21), 
             (20, 21), (20, 22), (20, 23), (20, 27), (20, 28), (20, 30), (20, 31), 
             (20, 34), (21, 23), (21, 24), (21, 27), (21, 29), (21, 30), (21, 32), 
             (21, 34), (22, 23), (22, 25), (22, 28), (22, 31), (23, 24), (23, 25), 
             (23, 26), (23, 28), (23, 29), (23, 31), (23, 32), (23, 33), (23, 34), 
             (24, 26), (24, 29), (24, 32), (25, 26), (25, 33), (26, 33), (27, 28), 
             (27, 29), (27, 30), (27, 34), (28, 29), (28, 31), (28, 34), (29, 32), 
             (29, 34), (30, 31), (30, 32), (30, 34), (30, 35), (31, 32), (31, 33), 
             (31, 34), (31, 35), (32, 33), (32, 34), (32, 35), (34, 35)]

def inverse_ab(ab):
    a = ab[0]
    b = ab[1]
    return (b,a)

NN1_edges2 = []
for edge in NN1_edges1:
    NN1_edges2.append(inverse_ab(edge))
    
NN1_edges = NN1_edges1 + NN1_edges2
      

#%%  
def get_possible_edges(nodes):
    
    edges = []
    nn = len(nodes) 
    for i in range(nn):
            for j in np.arange(i+1,nn):
                edges.append((nodes[i],nodes[j]))
    return edges
    
def get_edge_score(edges):
    
    count = 0
    for edge in edges:
        if edge in NN1_edges: count = count+1
    ratio = count/len(edges)
    
    return ratio



def connect_score(occ_nodes):
    
    occ_edges = get_possible_edges(occ_nodes)
    score = get_edge_score(occ_edges)
    
    return 1-score

def connect_score_2(occ_nodes):
    
    count = 0
    for i, nodei in enumerate(occ_nodes):
        nNN1 = 0 
        edges = []
        for j,nodej in enumerate(occ_nodes):
            if not nodei == nodej: edges.append((nodei, nodej))
        for edge in edges:
            if edge in NN1_edges: nNN1 = nNN1+1  
        #print(nNN1)
            
        if nNN1 >= 3: count = count +1 

    #print(count)
    score = 1 - count/len(occ_nodes)
    
    return score
'''
occ_nodes = [0,1,2,21]
occ_edges = get_possible_edges(occ_nodes)
score = get_edge_score(occ_edges)

score2 = connect_score_2(occ_nodes)
'''