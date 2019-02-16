# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 16:22:52 2018

@author: wangyf
"""

import networkx as nx
from networkx.algorithms import isomorphism as iso
import matplotlib.pyplot as plt
from pprint import pprint

import json
import numpy as np
import lattice_functions as lf
import pickle
from structure_constants import mother, dz, config, Ec, node_layer_dict

[Gm, Gsv, Gcv1, Gcv2, Gcv3] = pickle.load(open("clusters.p", "rb"))


def get_iso_config(config_list, i_config, drawing_flag = False):    
    '''
    Take each configuration (the configuration list and its index)
    check isomorphoric subgraphs
    and return a dictionary of nodes list
    '''
    config_i = config_list[i_config]
    n_nodes = len(config_i)
    n_layers = lf.get_layers(mother, dz, config_i)
    
    node_index = []
    for i in range(n_layers):
        node_index = node_index + node_layer_dict[i]
    
    sub_mother = mother[np.array(node_index)]
    
    
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
        Clusters.get_mother(sub_mother, dz)
        
        return Clusters
    
    Clusters = initialize_Clusters_object()

    empty = 'grey'
    filled = 'r'
    occ = [empty, filled]
    
    # Generate the mothe graph
    G1 =  Clusters.Gm
    # Generate the configuration graph
    Clusters.get_clusters(sub_mother, [config_i]) #one in layer 3 and one in layer 4
    Gcv = Clusters.Gcv
    G2 = Gcv[0]
    
    if drawing_flag == True:
        plt.figure()
        plt.title('G1')
        lf.drawing(G1)
        
        plt.figure()
        plt.title('G2')
        lf.drawing(G2)
    
    # Detect isomorphirc subgraphs 
    # the matching graphs are stored as key-value dictionary pairs
    if len(G2) > 1:
        
        GMz= iso.GraphMatcher(G1, G2, edge_match= iso.categorical_edge_match(['z'],[1.0])  )
        GMl = iso.GraphMatcher(G1, G2, edge_match= iso.numerical_edge_match(['length'],[1.0])  )
        iso_matches = [y for y in GMz.subgraph_isomorphisms_iter() if y in GMl.subgraph_isomorphisms_iter()]
        
    else:
        GMn = iso.GraphMatcher(G1, G2, node_match= iso.categorical_edge_match(['z'],[1]) )
        iso_matches = [y for y in GMn.subgraph_isomorphisms_iter()]
    
    # Use set and np.unqiue to eliminate the repeated graphs
    # return the graphs in indices list
    niso = len(iso_matches)
    iso_indices = [list(xi.keys()) for xi in iso_matches]
    iso_indices = [list(yi) for yi  in list(set(xi) for xi in iso_indices)]
    iso_indices = [list(xi) for xi in np.unique(iso_indices, axis = 0)]
    
    #Take all configurations with points fall on the right panel
    # as cluster expansion can handle sysmetric graphs
    #check if x of all the point > 0
    iso_indices_pos = []
    for iso_i in iso_indices:
        if np.any(mother[np.array(iso_i)][:,0] >= 0):
            iso_indices_pos.append([int(xi) for xi in iso_i])
    
    # save to a json file
    output_dict = {'configuration': config_i,
                   'index': i_config,
                   'n_nodes': n_nodes,
                   'n_layers': n_layers,
                   'n_iso': niso,
                   'iso_graph_list': iso_indices_pos}
    
    with open('iso_config_' + str(i_config) +'.json', 'w') as outfile:
        json.dump(output_dict, outfile)

#%% Main part of the function
get_iso_config(config, 1)  