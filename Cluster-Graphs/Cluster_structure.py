#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 10:36:45 2018

@author: wangyifan
"""

import networkx as nx
import matplotlib.pyplot as plt
import itertools as iter
import numpy as np
import re

class Pdn: 
    
    def __init__(self,n,input_edge,input_wt_i):
        self.n = n
        self.old_edge = input_edge
        self.old_wt_i = input_wt_i
    
    '''
    # N+1 graph connection generation function 
    # Input n -number of nodes in the new graphs
    #       edgelist and weighted index arrary from n-1 graph
    # return ni_u - number of unique graphs
    #       new edge list
    #       weight 2 edge in the list
    '''
        
    def new_graph(self):
        self.new_edge_add = [] # new possible edges that can be added to the graph'''
        self.new_edge_full = [] # new possible full edge lists '''
        self.H  = [] # edge list without isophoric graphs, return value 1 '''
        self.old_nhx = Pdn.check_hexagon(self.old_edge)[0]
        self.nhx = []
        self.ns = [] #neighboring structure list
        
        poss_edge = [] #the possible edge that can form triangles by adding one node'''
        
        '''
        Generate the possible edge that can form triangles 
        by taking out one with two triangles on the both side
        '''
        
        for i in range(len(self.old_edge)):
            if not np.any(i == np.array(self.old_wt_i)):
                poss_edge.append(self.old_edge[i])
                
        #print(poss_edge)
                
        '''number of possible graphs including isophoric graphs'''
        ni = len(poss_edge) 
        for i in range(ni):   
            self.new_edge_add.append([(self.n, poss_edge[i][0]), (self.n, poss_edge[i][1])])
        for i in range(ni):
            self.new_edge_full.append(self.old_edge + self.new_edge_add[i])
            
        #print('before\n')
        #print(self.new_edge_full)
        new_edge_full_old = []
        for i in range(ni):
            new_edge_full_old.append(self.new_edge_full[i])
        self.new_edge_full = Pdn.sort_list(self.new_edge_full)
        #print('before\n')
        #print(new_edge_full_old)
        #print('after\n')
        #print(self.new_edge_full)
        poss_edge = Pdn.sort_poss_edge(new_edge_full_old, self.new_edge_full, poss_edge)
        
        
        #print(poss_edge)
        
        
        i_list = []
        for i in range(len(self.new_edge_full)):
            i_list.append(i)
        i_list = Pdn.clear_isomorphic(self.new_edge_full)

        '''
        Update new graph edge list and weight 2 edge list
        weight = 2, 2 triangles already formed on both sides of the edge 
        '''          
        self.ni_u = len(i_list) #  number of unique graphs'''
        self.wt_i = [[] for j in range(len(i_list))] # initialize wt index array'''
        for i in range(self.ni_u):
                self.H.append(self.new_edge_full[i_list[i]])
        for i in range(self.ni_u):
            self.wt_i[i].append(self.H[i].index(poss_edge[i_list[i]]))
            for j in range(len(self.old_wt_i)):
                self.wt_i[i].append(self.old_wt_i[j])
        ''' # Check if there is hexagon formed in the graph'''
        for i in range(self.ni_u):
            reconnect = Pdn.connect_hexagon(self.H[i], self.wt_i[i], self.old_nhx)
            self.H[i] = reconnect[0]
            self.wt_i[i] = reconnect[1]
            self.nhx.append(reconnect[2])
        for i in range(self.ni_u):
            self.ns.append(neighbor_list_string(self.H[i]))
    '''    
    # check isophorism function
    # Input index pair (a,b) and full edge list edge_l
    # if two graphs are isoomorphic (identical), return 1
    '''
    @staticmethod
    def check_isomorphic(index,edge_list):
        '''
        # a is the index in edge list of graph A
        # b is the index in edge list of graph B
        '''
        a = index[0]
        b = index[1]
        A = nx.Graph(edge_list[a])
        B = nx.Graph(edge_list[b])
        if nx.is_isomorphic(A,B):
            return 1
    '''
    Remove isomorphic graphs from the list
    Return only the indices unique graphs
    '''    
    @staticmethod
    def clear_isomorphic(edge_list):
        '''
        n - number of struture needed to check
        edge_list - the full list of all strutures
        '''
        ni = len(edge_list) # number of isomorphy
        i_list = [] #  index of possible graphs including isophoric graphs'''
        for i in range(ni):
            i_list.append(i)
        index_pair = list(iter.combinations(i_list,2)) # combination of possible index pair
        is_iso = [] # Array of whether each pair are isophoric, 1 = yes, None = no'''
        for i in range(len(index_pair)):
            is_iso.append(Pdn.check_isomorphic(index_pair[i],edge_list))
        '''Remove isophoric graphs from the overall index'''
        for i in range(len(index_pair)):
            if is_iso[i] == 1:
                if index_pair[i][1] in i_list:
                    i_list.remove(index_pair[i][1])
        return i_list
        
    @staticmethod    
    def check_hexagon(edge_list):
        G = nx.Graph(edge_list)
        n_nodes = len(G)
        nn = []  # the number of nearest neigbhor for each node'''
        for i in range(n_nodes):
            nn.append(len(list(G.neighbors(i+1))))
        nhx = nn.count(6)
        return nhx,nn
    '''    
    # check if there is hexagon in the graph
    # Input, the edge list of a graph/weight 2 index/# of hexagon embedded in the plot
    '''
    @staticmethod
    def connect_hexagon(edge_list, wt_i, old_nhx):
        '''
        check if there is more hex detected  
        (center nodes with 6 connections is greater than 1)
        
        select the side nodes which connects to the center 
        (only two will be selected)
        ''' 
        [nhx,nn] = Pdn.check_hexagon(edge_list)
    
        if nhx == old_nhx+1:
            #nb_nn = []
            s_node = []
            G = nx.Graph(edge_list)
            ct_node = [i+1 for i, x in enumerate(nn) if x == 6] #the center of hexa node
            ct_node = ct_node[-1] # take the uncounted center node'''
            ct_neighbors = list(G.neighbors(ct_node)) # find the neighbors of the center
            # only nodes with two connects need to be connected as side nodes
            for i in range(len(ct_neighbors)): 
                wt_poss_edge = (ct_neighbors[i],ct_node)
                if wt_poss_edge in edge_list:
                    if not edge_list.index(wt_poss_edge) in wt_i:
                        s_node.append(ct_neighbors[i])
                        
                wt_poss_edge = (ct_node,ct_neighbors[i])
                if wt_poss_edge in edge_list:
                    if not edge_list.index(wt_poss_edge) in wt_i:
                        s_node.append(ct_neighbors[i])
            edge_list.append((s_node[0],s_node[1]))
            
            wt_edge = (s_node[0],ct_node)
            if wt_edge in edge_list:
                wt_i.append(edge_list.index(wt_edge))
                
            wt_edge = (ct_node,s_node[0])
            if wt_edge in edge_list:
                wt_i.append(edge_list.index(wt_edge))
            
            wt_edge = (s_node[1],ct_node)
            if wt_edge in edge_list:
                wt_i.append(edge_list.index(wt_edge))
                
            wt_edge = (ct_node,s_node[1])
            if wt_edge in edge_list:
                wt_i.append(edge_list.index(wt_edge))    
                
                
            return (edge_list,wt_i,nhx)
            '''
            for i in range(len(ct_neighbors)): 
            
                nb_nn.append(len(list(G.neighbors(ct_neighbors[i]))))
            for i in range(len(ct_neighbors)):
                if nb_nn[i] == 2 :
                    s_node.append(ct_neighbors[i]) #the side nodes need to be connected
            if not len(s_node) == 2:
                return (edge_list,wt_i,nhx)
            else:     
                edge_list.append((s_node[0],s_node[1])) #  connect the two side nodes
                #  add weight to center and the side nodes edge
                wt_i.append(edge_list.index((s_node[0],ct_node)))
                wt_i.append(edge_list.index((s_node[1],ct_node)))
                return (edge_list,wt_i,nhx)
            '''
        else: 
            return (edge_list,wt_i,nhx)
    
    @staticmethod
    def sort_list(edge_list):
    
        # check if (n, n-1) in the list
        def bubble_sort(a, b, edge_list):
            
            m = len(edge_list)
            
            for i in range(m):
                
                for j in range(1,m-i):
                    
                    if (a,b) in edge_list[j] or (b,a) in edge_list[j]:
                        
                        if not (a,b) in edge_list[j-1] or not (b,a) in edge_list[j-1]:
                            
                            edge_list[j-1], edge_list[j] = edge_list[j], edge_list[j-1]
            return edge_list
                        
        # check if (n, n-2) in the list                
        def bubble_sort_v2(a,  b,  edge_list):
            m = len(edge_list)
            
            for i in range(m):
                
                for j in range(1,m-i):
                    
                    if (a,b) in edge_list[j] or (b,a) in edge_list[j]:
                        if (a, b-1) in edge_list[j] or (b-1, a) in edge_list[j]:
                            if (a, b) in edge_list[j-1] or (b, a) in edge_list[j-1]:
                                if not (a, b-1) in edge_list[j-1] or not (b-1, a) in edge_list[j-1]:
                                    edge_list[j-1], edge_list[j] = edge_list[j], edge_list[j-1]
            return edge_list
        
        # check if (n-1, n-2) in the list                
        def bubble_sort_v3(a,  b,  edge_list):
            m = len(edge_list)
            
            for i in range(m):
                
                for j in range(1,m-i):
                    
                    if (a,b) in edge_list[j] or (b,a) in edge_list[j]:
                        if (a, b-1) in edge_list[j] or (b-1, a) in edge_list[j]:
                            if (a-1, b-1) in edge_list[j] or (b-1, a-1) in edge_list[j]:
                                if (a, b) in edge_list[j-1] or (b, a) in edge_list[j-1]:
                                    if (a, b-1) in edge_list[j-1] or (b-1, a) in edge_list[j-1]:
                                        if not (a-1, b-1) in edge_list[j-1] or not (b-a) in edge_list[j-1]:
                                            edge_list[j-1], edge_list[j] = edge_list[j], edge_list[j-1]
            return edge_list
                        
        n = len(nx.Graph(edge_list[0]))
        el = bubble_sort(n, n-1, edge_list)
        el = bubble_sort_v2(n, n-1, el)
        el = bubble_sort_v3(n, n-1, el)
        
        return el
    
    @staticmethod    
    def sort_poss_edge(edge_list_old, edge_list_new, poss_edge_old):
    
        i_list = []
        poss_edge_new = []
        ni = len(edge_list_old)
        
        for i in range(ni):
            i_list.append(edge_list_old.index(edge_list_new[i]))
            
        #print(i_list)
        for i in range(ni):
            poss_edge_new.append(poss_edge_old[i_list[i]])
        
            
        return poss_edge_new

    


def neighbor_list_string(edge_list):
    
    '''
    # function Generate neighboring list for zacors
    '''
    
   
    ll = edge_list
    l = []
    for i in range(len(ll)):
        l.append('%d-%d ' %(ll[i][0],ll[i][1]))
    s = ''.join(l) 
    return s

def neighbor_string_list(s):
    '''
    # function takes neighboring string from zacros and 
    # converts it to neighboring edge list
    '''
    s = s[0]
    ll = []
    p = re.compile(r'\W+')
    B  = p.split(s)
    nb = int(len(B)/2)
    for i in range(nb):
        bi = 2* i 
        ll.append((int(B[bi]), int(B[bi+1])))
    return ll
    
def neighbor_weight(edge_list):
    '''
    # function takes neighboring edge list and 
    # returns the weighted edge index
    '''
    
    n_edge = len(edge_list)
    G = nx.Graph(edge_list)
    wt_i = []
    
    for i in range(n_edge):
        pt_a = edge_list[i][0]
        pt_b = edge_list[i][1]
        
        list_a = list(G.neighbors(pt_a))
        list_b = list(G.neighbors(pt_b))
        wt = len(list(set(list_a).intersection(list_b)))
        if wt == 2: 
            wt_i.append(i)
    
    return wt_i

    

class Parallel_growth:
    def __init__(self,n, ni_u, input_edge_list, input_wt_i_list):
        self.n = n
        self.ni_u = ni_u
        self.H_list = input_edge_list
        self.wt_i_list = input_wt_i_list    
    
    def unique(self):
        
        self.H_u = []
        self.wt_i_u = []
        self.H_all = []
        self.wt_i_all = []
        self.ns = []
        for i in range(self.ni_u):
            Pdi = Pdn(self.n, self.H_list[i], self.wt_i_list[i])
            Pdi.new_graph()
            nj_u = Pdi.ni_u
            for j in range(nj_u):
                self.H_all.append(Pdi.H[j])
                self.wt_i_all.append(Pdi.wt_i[j])
        i_list = Pdn.clear_isomorphic(self.H_all)
        
        # Select the unqiue graphs
        self.mi_u = len(i_list)
        for i in range(self.mi_u):
            self.H_u.append(self.H_all[i_list[i]])
            self.wt_i_u.append(self.wt_i_all[i_list[i]])
            self.ns.append(neighbor_list_string(self.H_all[i_list[i]])) 
            
'''
# Plot unique cluster graphs
'''

def plot_graph(edge_list):      
    ni = len(edge_list)
    n = len(nx.Graph(edge_list[-1]))
    
    for j in range(ni//2):        
            plt.figure()
            for i in range(2):
                plt.subplot(1,2, i+1)
                G = nx.Graph(edge_list[2*j+i])
                nx.draw(G,with_labels=True, alpha=.75)
                plt.title('S %d' %(2*j+i+1))
            plt.suptitle('Pd %d' %n)
            plt.savefig('Pd%d_%d_structures_%d.png' %(n,ni,j+1))
    if not ni%2 == 0:
        plt.figure()
        G = nx.Graph(edge_list[-1])
        nx.draw(G,with_labels=True, alpha=.75)
        plt.title('Pd %d' %n)
        plt.suptitle('S %d' %(ni))
        plt.savefig('Pd%d_%d_structures_%d.png' %(n,ni,ni//2+1))

def clean_up(edge_list):
        
        '''
        clean up identical structure
        n - number of struture needed to check
        edge_list - the full list of all strutures
        '''
        ni = len(edge_list) # number of isomorphy
        i_list = [] #  index of possible graphs including isophoric graphs'''
        for i in range(ni):
            i_list.append(i)
        index_pair = list(iter.combinations(i_list,2)) # combination of possible index pair
        is_iso = [] # Array of whether each pair are isophoric, 1 = yes, None = no'''
        for i in range(len(index_pair)):
            
            a = index_pair[i][0]
            b = index_pair[i][1]
            is_iso.append(edge_list[a] == edge_list[b])
            
        '''Remove identical graphs from the overall index'''
        for i in range(len(index_pair)):
            if is_iso[i] == 1:
                if index_pair[i][1] in i_list:
                    i_list.remove(index_pair[i][1])
        return i_list          
        





