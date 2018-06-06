# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 12:09:40 2018

@author: wangyf
"""

import Cluster_structure as cs
import networkx as nx
import matplotlib.pyplot as plt
import itertools as iter
import numpy as np


def sort_list(edge_list, poss_edge):
    
    # check if (n, n-1) in the list
    def bubble_sort(a, b, edge_list):
        
        m = len(edge_list)
        
        for i in range(m):
            
            for j in range(1,m-i):
                
                if (a,b) in edge_list[j] or (b,a) in edge_list[j]:
                    
                    if not (a,b) in edge_list[j-1] or not (b,a) in edge_list[j-1]:
                        
                        edge_list[j-1], edge_list[j] = edge_list[j], edge_list[j-1]
        return edge_list
                    
    # check if (n-1, n-2) in the list                
    def bubble_sort_v2(a,  b,  edge_list):
        m = len(edge_list)
        
        for i in range(m):
            
            for j in range(1,m-i):
                
                if (a,b) in edge_list[j] or (b,a) in edge_list[j]:
                    if (a, b-1) in edge_list[j] or (b-1, a) in edge_list[j]:
                        if (a, b) in edge_list[j-1] or (b, a) in edge_list[j-1]:
                            if not (a, b-1) in edge_list[j-1] or (b-a) in edge_list[j-1]:
                                edge_list[j-1], edge_list[j] = edge_list[j], edge_list[j-1]
        return edge_list
                        
    n = len(nx.Graph(edge_list[0]))
    el = bubble_sort(n, n-1, edge_list)
    el = bubble_sort_v2(n, n-1, el)
    
    
    i_list = []
    ni = len(edge_list)
    
    for i in range(ni):
        i_list.append(edge_list.index(el[i]))
    poss_edge = list(np.array(poss_edge)[i_list])
        
    return el, poss_edge
    


n = 3
old_edge = [(1, 2), (1, 3) ,(2, 3)]
old_wt_i = []

cs.plot_graph([old_edge])


n = 4

old_edge = [(1, 2), (1, 3) ,(2, 3)]
old_wt_i = []
Pd4 = cs.Pdn(n, old_edge, old_wt_i)
Pd4.new_graph()
ni_u = Pd4.ni_u
new_H = Pd4.H
new_wt_i = Pd4.wt_i
print('Pd%d done' %n)
cs.plot_graph(Pd4.H)

n = 5
Pd5 = cs.Pdn(n, new_H[0], new_wt_i[0])
Pd5.new_graph()
ni_u = Pd5.ni_u
new_H = Pd5.H
new_wt_i = Pd5.wt_i
print('Pd%d done' %n)
cs.plot_graph(Pd5.H)
