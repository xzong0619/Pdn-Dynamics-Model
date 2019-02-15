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

from config import mother, Gm, Clusters, Gsv

empty = 'grey'
filled = 'r'
occ = [empty, filled]

G1 =  Gsv[0]
plt.figure()
lf.drawing(G1)

Clusters.get_clusters(mother, [[0]]) #one in layer 3 and one in layer 4
Gcv = Clusters.Gcv
G2 = Gcv[0]
plt.figure()
lf.drawing(G2)

if len(G2) > 1:
    
    GMz= iso.GraphMatcher(G1, G2, edge_match= iso.categorical_edge_match(['z'],[1.0])  )
    GMl = iso.GraphMatcher(G1, G2, edge_match= iso.numerical_edge_match(['length'],[1.0])  )
    x = [y for y in GMz.subgraph_isomorphisms_iter() if y in GMl.subgraph_isomorphisms_iter()]
    
else:
    GMn = iso.GraphMatcher(G1, G2, node_match= iso.categorical_edge_match(['z'],[1]) )
    x = [y for y in GMn.subgraph_isomorphisms_iter()]

niso = len(x)

#Cal = lf.calculations(occ)
#delta = Cal.get_delta(G1,G2)
