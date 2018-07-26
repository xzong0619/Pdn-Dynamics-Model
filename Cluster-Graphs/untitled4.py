# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 13:26:00 2018

@author: wangyf
"""

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism as iso

G1 =  nc_v[7]
print(G1.nodes(data=True))
plt.figure()
nx.draw(G1, with_labels=True)
G2 = gv[2]
print(G2.nodes(data=True))
plt.figure()
nx.draw(G2, with_labels=True)

GM = iso.GraphMatcher(Gm, G1)
x= [y for y in GM.subgraph_isomorphisms_iter()]

