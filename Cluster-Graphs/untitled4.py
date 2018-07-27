# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 13:26:00 2018

@author: wangyf
"""

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism as iso
from pprint import pprint

G1 =  Gsv[7]
plt.figure()
drawing(G1)
G2 = Gcv[4]
plt.figure()
drawing(G2)

GM = iso.GraphMatcher(G1, G2, edge_match=iso.numerical_edge_match(['length'],[2.0]))
x= [y for y in GM.subgraph_isomorphisms_iter()]

print(len(x))
pprint(x)