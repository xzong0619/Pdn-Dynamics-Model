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

from config import mother, Gm, Clusters


G1 =  Gm
plt.figure()
lf.drawing(G1)

Clusters.get_clusters(mother, [[34,35]]) #one in layer 3 and one in layer 4
Gcv = Clusters.Gcv
G2 = Gcv[0]
plt.figure()
lf.drawing(G2)

GM = iso.GraphMatcher(G1, G2, edge_match=   [iso.numerical_edge_match(['length'],[1]), iso.categorical_edge_match(['z'],[1])    ] )
x= [y for y in GM.subgraph_isomorphisms_iter()]

niso =len(x)

