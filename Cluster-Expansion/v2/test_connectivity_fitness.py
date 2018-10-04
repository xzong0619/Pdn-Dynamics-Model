# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 11:55:36 2018

@author: wangyf
"""

import GA_functions as GA
import networkx as nx
def initialize_Clusters_object():
    
    empty = 'grey'
    filled = 'r'
    occ = [empty, filled]
    
    '''
    only draw 1st nearest neighbors?
    '''
    NN1 = 1
    '''
    Draw mother/conifgurations/clusters?
    '''
    draw = [0, 1, 0]
    
    
    Clusters = lf.clusters(occ, NN1, draw)
    Clusters.get_mother(mother, dz)
    
    return Clusters
def get_color(Gsv):
    
    return [Gsv.nodes[i]['color'] for i in Gsv.nodes]
    
def get_occupied_nodes(colorsv):
    nodes = len(colorsv)
    boov =  [color == 'r' for color in colorsv]
    occ_nodes = np.arange(nodes)[np.array(boov)]

    return occ_nodes
    
def get_possible_edges(nodes):
    
    edges = []
    nn = len(nodes) 
    for i in range(nn):
            for j in np.arange(i+1,nn):
                edges.append((nodes[i],nodes[j]))
    return edges
    
def get_edge_score(edges, Gsv):
    
    count = 0
    for edge in edges:
        if edge in Gsv.edges: count = count+1
    score = count/len(edges)
    
    return score

Clusters = initialize_Clusters_object()
occ = Clusters.occupancy
'''
ind = [1,1,1,1]
ind_Gsv = GA.individual_config(ind, Clusters)[0]
'''
ind_list = [24,21,1,2,3]
Clusters.get_configs([ind_list])
ind_Gsv = Clusters.Gsv[0]
node_color = get_color(ind_Gsv)
occ_nodes = get_occupied_nodes(node_color)
occ_edges = get_possible_edges(occ_nodes)
score = get_edge_score(occ_edges, ind_Gsv)