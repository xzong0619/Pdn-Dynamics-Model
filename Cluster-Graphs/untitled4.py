# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 13:26:00 2018

@author: wangyf
"""

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism as iso

G1 =  nc_v[5]
print(G1.nodes(data=True))
plt.figure()
nx.draw(G1, with_labels=True)
G2 = gv[2]
print(G2.nodes(data=True))
plt.figure()
nx.draw(G2, with_labels=True)
GM = iso.GraphMatcher(G1,G2)
x= []
for mapping in GM.subgraph_isomorphisms_iter():
        print (mapping)
        x.append(mapping)
#x = [y for y in ]
#z = GM.mapping()
#x = [y for y in GM.isomorphisms_iter()]
#print(z)