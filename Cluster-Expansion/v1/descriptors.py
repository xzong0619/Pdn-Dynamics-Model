# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 17:33:26 2018

@author: wangyf
"""

import sys
import os
import networkx as nx

# working on all platforms, using cross-platform home path

HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Cluster-Graphs'))

import Cluster_structure as cs
from nc import nc


Pd1 = nx.Graph()
Pd1.add_node(1) 
nc_edge = []
nc_v = [Pd1]
gv = [Pd1]
edgev = []

for i in range(len(nc)):
    nc_edge.append(cs.neighbor_string_list([nc[i]]))    

for i in range(len(nc)):
    nc_v.append(nx.Graph(nc_edge[i])) 
    
    
dev = ['1-2',
       '1-2 2-3',
       '1-2 1-3 2-3',
       '1-2 1-3 2-3 4-2 4-3',
       '1-2 1-3 1-4',
       '1-2 2-3 1-3 1-4']


for i in range(len(dev)):
    edgev.append(cs.neighbor_string_list([dev[i]]))


for i in range(len(edgev)):
    gv.append(nx.Graph(edgev[i]))
    
cs.plot_graph(gv)