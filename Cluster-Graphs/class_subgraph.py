# -*- coding: utf-8 -*-
"""
Created on Sun Aug 12 16:29:25 2018

@author: yifan
"""
from itertools import combinations
import numpy as np
import lat_fun as lf



class subgraphs():
    '''
    generate subgraph list with the nodes numbers under the mother graph
    '''
    
    def __init__(self, mother):
        
       self.index= np.arange(len(mother)) # generate the index of nodes
       self.mother = mother
       
    @staticmethod
    def layer_tuple(mother, ci):
        
        '''
        takes in a combo of index and returns tuple of layers they are in 
        '''
    
        n = len(ci)
        index = []
        for i in range(n):
            index.append(ci[i])
            
        layers = []
        
        for i in range(n):
            layers.append(mother[index[i]][2])
            
        layers= tuple(layers)
            
        return layers
    
    @staticmethod   
    def distance_tuple(mother, ci):
        '''
        takes in a combo of index and returns sorted distance between nodes
        '''
        n = len(ci)
        index = []
        
        for i in range(n):
            index.append(ci[i])
        
        combo = list(combinations(index,2))
        ncombo = len(combo) #0 for 1 node, 2 for 2 nodes, 3 for 3 nodes
          
        distances = []
        
        for i in range(ncombo):
            pt1 = mother[combo[i][0]]
            pt2 = mother[combo[i][1]]
            distances.append(lf.two_points_D(pt1, pt2))
            
        distances = tuple(sorted(distances))
        
        return distances
    
    def get_s(self, n_atoms):
        
        '''
        Input number of nodes in a subgraph
        Generate combinations among the nodes
        '''
        self.n_atoms = n_atoms
        
        combo = list(combinations(self.index, self.n_atoms))
        ncombo  = len(combo)
        
        
        '''
        generate the information list
        store the sorted distance of nodes in tuple 1
        + the layer each node is in in tuple 2
        '''
        
        info = [] 
        
        for i in range(ncombo):
            ci  = combo[i]
            
            distances = self.distance_tuple(self.mother, ci)
            layers = self.layer_tuple(self.mother, ci)
            
            info.append((distances, layers))
        
        info_set = list(set(info))
        
        index_list =[]
        
        for i in info_set:
            index_list.append(info.index(i))
            
        index_list.sort() # sort the list and take out those indices
        
        s_list = np.array(combo)[index_list]
            
        return s_list



