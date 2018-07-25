# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 13:26:00 2018

@author: wangyf
"""

import networkx as nx
from networkx.algorithms import isomorphism as isomo

G1 =  nc_v[2]
plt.figure()
nx.draw(G1, with_labels=True)
G2 = gv[1]
plt.figure()
nx.draw(G2, with_labels=True)
GM = isomo.GraphMatcher(G1,G2)
x = GM.isomorphisms_iter()