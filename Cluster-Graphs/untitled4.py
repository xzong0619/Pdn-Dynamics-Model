# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 13:26:00 2018

@author: wangyf
"""

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism as iso

G1 =  Gsv[7]
print(G1.nodes(data=True))
plt.figure()
drawing(G1)
G2 = Gcv[4]
print(G2.nodes(data=True))
plt.figure()
drawing(G2)

GM = iso.GraphMatcher(G1, G2, edge_match=iso.categorical_edge_match(['length'],[]))
x= [y for y in GM.subgraph_isomorphisms_iter()]

