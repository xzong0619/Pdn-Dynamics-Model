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
    NN1 = 0
    '''
    Draw mother/conifgurations/clusters?
    '''
    draw = [0, 0, 0]
    
    
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
    

def get_neighbor_colors(node, Gsv):
    
    neighbors = list(nx.all_neighbors(Gsv, i))
    edges = 
    
Clusters = initialize_Clusters_object()
occ = Clusters.occupancy

ind = [0,1,1,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0]
ind_Gsv = GA.individual_config(ind, Clusters)[0]

node_color = get_color(ind_Gsv)
occ_nodes = get_occupied_nodes(node_color)