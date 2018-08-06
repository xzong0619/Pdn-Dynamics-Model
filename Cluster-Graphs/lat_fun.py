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
from numpy.linalg import inv

#%%
'''
Basic functions
'''

def two_points_D(A,B):
    '''
    Calculates the distance between two points A and B
    '''
    d = math.sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)
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
    
#%%
'''
clusters object
'''

class clusters():
    
    def __init__(self, occupancy):
        
        '''
        takes in the occupancy color vector 
        occupancy[0] is the empty color
        occupancy[1] is the filled color
        '''
        
        self.occupancy = occupancy
        self.empty = self.occupancy[0]
        self.filled = self.occupancy[1]
        
    def gmothers(self, mother):
    
        '''
        takes in mother cooridate list 
        returns connected lattice graph
        '''
        self.mother = mother
        self.nm = len(mother)
        Gm = nx.Graph()
        
        for i in range(self.nm):
            Gm.add_node(i, pos = mother[i], color = self.empty)
        
        self.edge_d = []
        self.edge = []
        
        # Add all egdes and calculate the edge distance
        for i in range(self.nm):
            for j in np.arange(i+1,self.nm):
                self.edge.append((i,j))
                self.edge_d.append(two_points_D(mother[i],mother[j]))
                
        self.ne = len(self.edge)
        for i in range(self.ne): 
            if self.edge_d[i] <= 1:
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
        print('Creating G')
        Gs = nx.Graph()
        print('Adding nodes')
    
        for i in range(self.nm):
            Gs.add_node(i, pos = self.mother[i], color = self.empty)

        print('Adding edges')
        for i in range(self.ne):
            if self.edge_d[i] <= 1:
                 Gs.add_edges_from([self.edge[i]], length = self.edge_d[i])   
        print('after')
        print(Gs.nodes)
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
            Gc.add_node(i, pos = cmother[c], color = self.filled)
            
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


