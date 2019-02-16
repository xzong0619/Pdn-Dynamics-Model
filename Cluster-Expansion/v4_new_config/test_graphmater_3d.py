# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 16:22:52 2018

@author: wangyf
"""

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism as iso
from pprint import pprint
import numpy as np
import lattice_functions as lf
import pickle
from structure_constants import mother, dz, config, Ec, node_layer_dict
import json

[Gm, Gsv, Gcv1, Gcv2, Gcv3] = pickle.load(open("clusters.p", "rb"))

layer = 1
node_index = []
for i in range(layer):
    node_index = node_index + node_layer_dict[i]

mother = mother[np.array(node_index)]

#%%
def initialize_Clusters_object():
    
    empty = 'grey'
    filled = 'r'
    occ = [empty, filled]
    
    '''
    only draw 1st nearest neighbors?
    '''
    NN1 = 0
    '''
    Draw mother/conifgurations/clusters?
    '''
    draw = [0, 0, 0]
    
    
    Clusters = lf.clusters(occ, NN1, draw)
    Clusters.get_mother(mother, dz)
    
    return Clusters

Clusters = initialize_Clusters_object()



empty = 'grey'
filled = 'r'
occ = [empty, filled]

G1 =  Clusters.Gm
#plt.figure()
#lf.drawing(G1)

Clusters.get_clusters(mother, [[0,1,2]]) #one in layer 3 and one in layer 4
Gcv = Clusters.Gcv
G2 = Gcv[0]
#plt.figure()
#lf.drawing(G2)

#%%
if len(G2) > 1:
    
    GMz= iso.GraphMatcher(G1, G2, edge_match= iso.categorical_edge_match(['z'],[1.0])  )
    GMl = iso.GraphMatcher(G1, G2, edge_match= iso.numerical_edge_match(['length'],[1.0])  )
    x = [y for y in GMz.subgraph_isomorphisms_iter() if y in GMl.subgraph_isomorphisms_iter()]
    
else:
    GMn = iso.GraphMatcher(G1, G2, node_match= iso.categorical_edge_match(['z'],[1]) )
    x = [y for y in GMn.subgraph_isomorphisms_iter()]


niso = len(x)
iso_indices = [list(xi.keys()) for xi in x]
iso_indices = [list(yi) for yi  in list(set(xi) for xi in iso_indices)]

iso_indices = [list(xi) for xi in np.unique(iso_indices, axis = 0)]


#%%
#Take all configurations with points fall on the left panel
#check if x of all the point > 0
iso_indices_pos = []
for iso_i in iso_indices:
    if np.any(mother[np.array(iso_i)][:,0] >= 0):
        iso_indices_pos.append([int(xi) for xi in iso_i])

with open('iso.json', 'w') as outfile:
    json.dump(iso_indices_pos, outfile)
#Cal = lf.calculations(occ)
#delta = Cal.get_delta(G1,G2)
