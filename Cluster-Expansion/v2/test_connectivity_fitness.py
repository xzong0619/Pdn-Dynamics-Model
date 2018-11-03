# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 11:55:36 2018

@author: wangyf
"""
<<<<<<< HEAD


import networkx as nx
import numpy as np
import GA_functions as GA
=======
import numpy as np
import GA_functions as GA


>>>>>>> aa3295d97abc46f1e6bcdf3dd7e92d6fbc929afd
import lattice_functions as lf
from structure_constants import mother, dz

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
<<<<<<< HEAD
    draw = [0,0, 0]
=======
    draw = [0, 0, 0]
>>>>>>> aa3295d97abc46f1e6bcdf3dd7e92d6fbc929afd
    
    
    Clusters = lf.clusters(occ, NN1, draw)
    Clusters.get_mother(mother, dz)
    
    return Clusters

def get_color(Gsv):
    
    return [Gsv.nodes[i]['color'] for i in Gsv.nodes]
    
<<<<<<< HEAD
def get_occupied_nodes(colorsv, occ):
    nodes = len(colorsv)
    boov =  [color == occ[1] for color in colorsv]
=======
def get_occupied_nodes(colorsv, Clusters):
    nodes = len(colorsv)
    boov =  [color == Clusters.occupancy[1] for color in colorsv]
>>>>>>> aa3295d97abc46f1e6bcdf3dd7e92d6fbc929afd
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
<<<<<<< HEAD
    ratio = count/len(edges)
    
    return 1-ratio

def connect_score(individual):
 
    Clusters = initialize_Clusters_object()
    ind_Gsv = GA.individual_config(individual, Clusters)[0]
    node_color = get_color(ind_Gsv)
    occ_nodes = get_occupied_nodes(node_color, Clusters.occupancy)
=======
    percentage = count/len(edges)
    
    return 1-percentage


def connect_score(individual):
    
    Clusters = initialize_Clusters_object()
    ind_Gsv = GA.individual_config(individual, Clusters)[0]
    node_color = get_color(ind_Gsv)
    occ_nodes = get_occupied_nodes(node_color, Clusters)
>>>>>>> aa3295d97abc46f1e6bcdf3dd7e92d6fbc929afd
    occ_edges = get_possible_edges(occ_nodes)
    score = get_edge_score(occ_edges, ind_Gsv)
    
    return score



<<<<<<< HEAD
=======







>>>>>>> aa3295d97abc46f1e6bcdf3dd7e92d6fbc929afd
'''
ind = [1,1,1,1]
ind_Gsv = GA.individual_config(ind, Clusters)[0]
'''
<<<<<<< HEAD

'''
Clusters = initialize_Clusters_object()
=======
'''
>>>>>>> aa3295d97abc46f1e6bcdf3dd7e92d6fbc929afd
ind_list = [24,21,1,2,3]
Clusters.get_configs([ind_list])
ind_Gsv = Clusters.Gsv[0]
node_color = get_color(ind_Gsv)
<<<<<<< HEAD
occ_nodes = get_occupied_nodes(node_color, Clusters.occupancy)
=======
occ_nodes = get_occupied_nodes(node_color, Clusters)
>>>>>>> aa3295d97abc46f1e6bcdf3dd7e92d6fbc929afd
occ_edges = get_possible_edges(occ_nodes)
score = get_edge_score(occ_edges, ind_Gsv)
'''