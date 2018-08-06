# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 11:08:30 2018

@author: wangyf
"""
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import math
from networkx.algorithms import isomorphism as iso
from numpy.linalg import inv


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
    
def gmothers(mother, occupancy):
    
    '''
    takes in mother cooridate list and occupancy vector 
    occupancy[0] is the empty color
    occupancy[1] is the filled color
    returns connected lattice graph
    '''
    
    nm = len(mother)
    Gm = nx.Graph()
    
    for i in range(nm):
        # take only the first two coordinates x and y
        Gm.add_node(i, pos = mother[i], color = occupancy[0])
    
    edge_d = []
    edge = []
    
    # Add all egdes and calculate the edge distance
    for i in range(nm):
        for j in np.arange(i+1,nm):
            edge.append((i,j))
            edge_d.append(two_points_D(mother[i],mother[j]))
            
    ne = len(edge)
    for i in range(ne):
        if edge_d[i] <= 1:
            Gm.add_edges_from([edge[i]], length = edge_d[i])
    drawing(Gm)
    plt.title('%d lattice points' %nm)
    
    return Gm
    
def gconfigurations(mother, son, occupancy):
    '''
    takes in mother coordinate list and son's index number and occupancy vector
    returns the shaded son graph
    '''
    Gs = nx.Graph()
    nm = len(mother)
    ns = len(son)
    
    
    edge_d = []
    edge = []
    for i in range(nm):
        for j in np.arange(i+1,nm):
            edge.append((i,j))
            edge_d.append(two_points_D(mother[i],mother[j]))
            
    ne = len(edge)

    for i in range(nm):
        Gs.add_node(i, pos = mother[i], color = occupancy[0])
    for i in range(ne):
        if edge_d[i] <= 1:
            Gs.add_edges_from([edge[i]], length = edge_d[i])    
    for si in range(ns):
        Gs.node[son[si]]['color'] = occupancy[1]
    
    drawing(Gs)
    plt.title('Pd %d' %ns)

    return Gs

def gclusters(cmother, son, occupancy):
    ns = len(son)
    Gc = nx.Graph()
    for i in range(ns):
        c = son[i]
        Gc.add_node(i, pos = cmother[c], color = occupancy[1])
        
    edge_d = []
    edge = []
    for i in range(ns):        
        for j in np.arange(i+1,ns):
            c  = son[i]
            d = son[j]
            edge.append((i,j))
            edge_d.append(two_points_D(cmother[c],cmother[d]))    
    ne = len(edge)
    for i in range(ne):
        Gc.add_edges_from([edge[i]], length = edge_d[i])
        
    drawing(Gc)
    plt.title('Pd %d' %ns)
       
    return Gc
 
def get_occupancy(G, i, occupancy):
    
    '''
    Get the occupancy from the graph G for node i 
    Occupied is 1 and unoccupied is 0
    '''
    empty = occupancy[0]
    filled = occupancy[1]
    
    if G.nodes[i]['color'] == empty: o = 0
    if G.nodes[i]['color'] == filled: o = 1 
    
    return o    

def get_delta(Gl,Gs, occupancy):
    
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
            subi[i].append(get_occupancy(Gl,subg[i][j], occupancy))   
        subs.append(np.product(subi[i]))
    delta = np.sum(subs)/niso
    
    return delta



def get_J(Ev, G1v, G2v, occupancy):
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
            pi[i][j] = get_delta(G1v[i],G2v[j], occupancy)
            
    J = np.linalg.lstsq(pi, Ev)[0]
    
    return J, pi

def LeaveOneOut(A, a):
    '''
    takes in a list A and returns a new list B by leaving ath element out
    '''     
    B = [x for i,x in enumerate(A) if i!=a]
    
    return B 
    
