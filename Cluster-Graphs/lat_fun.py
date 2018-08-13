# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 11:08:30 2018

@author: wangyf
"""
'''
Objective Oriented Version of Lattice Building functon
'''

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import math
from networkx.algorithms import isomorphism as iso
from itertools import combinations

#%%
'''
Basic functions
'''

def two_points_D(A,B):
    
    '''
    Calculates the distance between two points A and B
    '''
    n = len(A)
    s = 0
    for i in range(n):
        s = s+ (A[i]-B[i])**2
    d = math.sqrt(s)
    return d


def drawing(G):
    '''
    takes in graph g, draw it
    '''
    color = nx.get_node_attributes(G,'color')
    pos=nx.get_node_attributes(G,'pos')
    plt.figure() 
    nx.draw(G, pos, with_labels=True, node_color = list(color.values()))
    
    
def LeaveOneOut(A, a):
    '''
    takes in a list A and returns a new list B by leaving ath element out
    '''     
    B = [x for i,x in enumerate(A) if i!=a]
    
    return B 

def add_z(v, z):
    '''
    takes in a np array and add z coordinates to it
    '''
    vd = np.concatenate((v, np.array([z*np.ones(len(v))]).T), axis =1) 

    return vd
    
#%%
'''
clusters object
'''

class clusters():
    
    def __init__(self, occupancy, NN1):
        
        '''
        takes in the occupancy color vector 
        occupancy[0] is the empty color
        occupancy[1] is the filled color
        '''
        
        self.occupancy = occupancy
        self.empty = self.occupancy[0]
        self.filled = self.occupancy[1]
        self.NN1 = NN1
        
    def gmothers(self, mother):
    
        '''
        takes in mother cooridate list 
        returns connected lattice graph
        '''
        self.mother = mother
        self.nm = len(mother)
        Gm = nx.Graph()
        
        for i in range(self.nm):
            Gm.add_node(i, pos = mother[i][:2], color = self.empty)
        
        self.edge_d = []
        self.edge = []
        
        # Add all egdes and calculate the edge distance
        for i in range(self.nm):
            for j in np.arange(i+1,self.nm):
                self.edge.append((i,j))
                self.edge_d.append(two_points_D(mother[i],mother[j]))
                
        self.ne = len(self.edge)
        for i in range(self.ne): 
            if self.NN1: # only draw 1st Nearest Neighbors 
                if self.edge_d[i] <= 1: 
                    Gm.add_edges_from([self.edge[i]], length = self.edge_d[i])
            else:
                Gm.add_edges_from([self.edge[i]], length = self.edge_d[i])
                
        drawing(Gm)
        plt.title('%d lattice points' %self.nm)
        return Gm
    
    
    def gconfigurations(self, son):
        
        '''
        takes in mother coordinate list and son's index number and occupancy vector
        returns the shaded son graph
        '''     
        ns = len(son)
        Gs = nx.Graph()

        for i in range(self.nm):
            Gs.add_node(i, pos = self.mother[i][:2], color = self.empty)

        for i in range(self.ne): 
            if self.NN1: # only draw 1st Nearest Neighbors 
                if self.edge_d[i] <= 1:
                    Gs.add_edges_from([self.edge[i]], length = self.edge_d[i])
            else:
                Gs.add_edges_from([self.edge[i]], length = self.edge_d[i])

        for si in range(ns):
            Gs.node[son[si]]['color'] = self.filled
        
        drawing(Gs)
        plt.title('Pd %d' %ns)
        
        return Gs
    

    def gclusters(self, cmother, cson):
    
        '''
        takes in clusters 
        return cluster graph objective
        '''
        Gc = nx.Graph()
        cns = len(cson)
        
        for i in range(cns):
            c = cson[i]
            Gc.add_node(i, pos = cmother[c][:2], color = self.filled)
            
        cedge_d = []
        cedge = []
        for i in range(cns):        
            for j in np.arange(i+1,cns):
                c  = cson[i]
                d = cson[j]
                cedge.append((i,j))
                cedge_d.append(two_points_D(cmother[c],cmother[d]))    
        
        cne = len(cedge)
        for i in range(cne):
           Gc.add_edges_from([cedge[i]], length = cedge_d[i])
            
        drawing(Gc)
        plt.title('Pd %d' %cns)
        return Gc
    

    def get_mother(self, mother):
        '''
        takes in mother coordinates list and 
        add mother attribute to the class
        '''
        
        self.Gm  = self.gmothers(mother)
        
        
    def get_configs(self, config):
        
        '''
        takes in configuration index list
        get a list of configurations as graph objects
        '''
        
        self.Gsv = []
        self.nconfig = len(config)
        
        for si in range(self.nconfig):
            son_i = config[si]
            Gs = self.gconfigurations(son_i)
            self.Gsv.append(Gs)
     
        
    def get_clusters(self, cmother, ccluster):
        
        '''
        takes in cluster coordinates list and cluster index list
        returns a list of clusters as graph objects
        '''
        
        self.nc = len(ccluster) # number of clusers
        self.Gcv = [] # list of clusters
        for si in range(self.nc):
            cson = ccluster[si] 
            Gc = self.gclusters(cmother,cson)
            self.Gcv.append(Gc)       
        

#%%        
class calculations():
    
    '''
    Perform statistical calculation for cluster expansion
    '''
    
    def __init__(self,occupancy):
        
        '''
        takes in the occupancy color vector 
        occupancy[0] is the empty color
        occupancy[1] is the filled color
        '''
        
        self.occupancy = occupancy
        self.empty = self.occupancy[0]
        self.filled = self.occupancy[1]
        
        
    def get_occupancy(self, G, i):
        
        '''
        Get the occupancy from the graph G for node i 
        Occupied is 1 and unoccupied is 0
        '''
        
        if G.nodes[i]['color'] == self.empty: o = 0
        if G.nodes[i]['color'] == self.filled: o = 1 
        
        return o    

    def get_delta(self, Gl,Gs):
        
        '''
        takes in larger graph Gl and smaller graph Gs
        find sub isomorphic graphs of Gs from Gl
        calculate the delta value in pi matrix 
        '''
        '''
        find subgraphs using edge distance match
        '''
        GM = iso.GraphMatcher(Gl, Gs, edge_match=iso.numerical_edge_match(['length'],[1.0]))
        '''
        list down total number of subgraphs niso
        '''
        x= [y for y in GM.subgraph_isomorphisms_iter()]
        
        niso =len(x)
        '''
        save subgraphs to a list
        '''
        subg = list()
        for i in range(niso):    
            subg.append(tuple(x[i].keys()))
        
        '''
        save product into a list 
        and caclulate the sum divid by total number of subgraphs
        '''
        subi = []
        subs = []
        for i in range(niso):
            subi.append([])
            for j in range(len(subg[i])):
                subi[i].append(self.get_occupancy(Gl,subg[i][j]))   
            subs.append(np.product(subi[i]))
        delta = np.sum(subs)/niso
        
        return delta
    
    
    
    def get_J(self, Ev, G1v, G2v):
        '''
        The function that gets 
            energy of configurations, Ev
            configuration graphs, G1v
            cluster graphs, G2v
        and returns the interaction energy in J vector and correlation matrix pi
        '''
        
        Ev = np.array(Ev)
        n1 = len(G1v)
        n2 = len(G2v)
        pi = np.zeros((n1,n2))
        
        for i in range(n1):
            for j in range(n2):
                pi[i][j] = self.get_delta(G1v[i],G2v[j])
                
        J = np.linalg.lstsq(pi, Ev)[0]
        
        return J, pi

#%%
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
