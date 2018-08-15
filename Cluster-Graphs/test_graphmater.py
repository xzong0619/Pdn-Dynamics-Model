# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 13:26:00 2018

@author: wangyf
"""

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism as iso
from pprint import pprint
import numpy as np
import lattice_functions as lf
def get_occupancy(G, i):
    
    '''
    Get the occupancy from the graph G for node i 
    Occupied is 1 and unoccupied is 0
    '''
    if G.nodes[i]['color'] == empty: o = 0
    if G.nodes[i]['color'] == filled: o = 1 
    
    return o

G1 =  Gsv[7]
plt.figure()
lf.drawing(G1)
G2 = Gcv[4]
plt.figure()
lf.drawing(G2)

GM = iso.GraphMatcher(G1, G2, edge_match=iso.numerical_edge_match(['length'],[1.0]))
x= [y for y in GM.subgraph_isomorphisms_iter()]

niso =len(x)

subg = list()
for i in range(niso):    
    subg.append(tuple(x[i].keys()))

subi = []
subs = []
for i in range(niso):
    subi.append([])
    for j in range(len(subg[i])):
        subi[i].append(get_occupancy(G1,subg[i][j]))   
    subs.append(np.product(subi[i]))
delta = np.sum(subs)/niso











