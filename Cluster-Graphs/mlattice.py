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
from pprint import pprint
from numpy.linalg import inv


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
    
def gmothers(mother):
    
    nm = len(mother)
    Gm = nx.Graph()
    
    for i in range(nm):
        Gm.add_node(i, pos = mother[i], color = empty)
    
    edge_d = []
    edge = []
    
    # Add all egdes and calculate the edge distance
    for i in range(nm):
        for j in np.arange(i+1,nm):
            edge.append((i,j))
            edge_d.append(two_points_D(mother[i],mother[j]))
            
    ne = len(edge)
    for i in range(ne):
        #if edge_d[i] <= 1:
        Gm.add_edges_from([edge[i]], length = edge_d[i])
    drawing(Gm)
    plt.title('%d lattice points' %nm)
    
    return Gm
    
def gconfigurations(mother, son):
    '''
    takes in mother coordinate list and son's index number
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
    for i in range(ne):
        Gs.add_edges_from([edge[i]], length = edge_d[i])
    for i in range(nm):
        Gs.add_node(i, pos = mother[i], color = empty)
    for i in range(ne):
        if edge_d[i] <= 1:
            Gs.add_edges_from([edge[i]], length = edge_d[i])    
    for si in range(ns):
        Gs.node[son[si]]['color'] = filled
    
    drawing(Gs)
    plt.title('Pd %d' %ns)

    return Gs

def gclusters(cmother, son):
    ns = len(son)
    Gc = nx.Graph()
    for i in range(ns):
        c = son[i]
        Gc.add_node(i, pos = cmother[c], color = filled)
        
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
 
def get_occupancy(G, i):
    
    '''
    Get the occupancy from the graph G for node i 
    Occupied is 1 and unoccupied is 0
    '''
    if G.nodes[i]['color'] == empty: o = 0
    if G.nodes[i]['color'] == filled: o = 1 
    
    return o    

def get_delta(Gl,Gs):
    
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
            subi[i].append(get_occupancy(Gl,subg[i][j]))   
        subs.append(np.product(subi[i]))
    delta = np.sum(subs)/niso
    
    return delta



def get_J(Ev, G1v, G2v):
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
            pi[i][j] = get_delta(G1v[i],G2v[j])
            
    J = np.linalg.lstsq(pi, Ev)[0]
    
    return J, pi

def LeaveOneOut(A, a):
    '''
    takes in a list A and returns a new list B by leaving ath element out
    '''     
    B = [x for i,x in enumerate(A) if i!=a]
    
    return B 
    
#%%
'''
main
'''
    

mother = [(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5),
          (1/2, -3**0.5/2), (3/2, -3**0.5/2), (0, -3**0.5),
          (-1,0), (-1/2, -3**0.5/2), (-3/2, -3**0.5/2),
          (-1/2, 3**0.5/2), (-3/2, 3**0.5/2)]
config = [[0],
          [0,1],
          [0,1,2],
          [0,1,2,3],
          [0,1,2,3,5],
          [0,1,2,3,5,6],
          [0,1,2,3,5,9],
          [0,1,2,3,5,11],
          [0,1,2,5,9,8,11],
          [0,1,2,3,5,9,7],
          [0,1,2,5,8,11,12],
          [0,1,2,8,9,11,12]]
'''
Energy of clusters in eV
'''

Ec = [0, -0.45, -1.5, -2.56, -3.41, -4.49, -4.38, -4.32, -6.07, -4.9, -5.11, -5.12]


     
empty = 'grey'
filled = 'r'

Gm = gmothers(mother)
    
        
'''
Creat 12 configurations
'''
ns = len(config)  # number of configurations 
Gsv = [] # list of configurations 
for si in range(ns):
    son = config[si]
    Gs = gconfigurations(mother,son)
    Gsv.append(Gs)


#%%
'''
Creat 6 clusters
'''

cmother = [(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5)]
cconfig = [[0], [0,1], [0,3], [1,4],[0,1,2],[1,2,4], [0,1,3]]


nc = len(cconfig) # number of clusers
Gcv = [] # list of clusters
for si in range(nc):
    cson = cconfig[si] 
    Gc = gclusters(cmother,cson)
    Gcv.append(Gc)       
        
#%% Stattistical analysis
'''
creat pi matrix
size of number of configuration * numbers of clusters
'''

J, pi =  get_J(Ec,Gsv,Gcv)      
MSE = np.sum(np.power((np.dot(pi,J) - Ec),2))/ns 



#%%
'''
Calculate Leave One Out CV score
'''
Ec_predict = np.zeros(ns)


for i in range(ns):
    '''
    i is the index to be removed from the list
    '''
    Ec_LOO =  LeaveOneOut(Ec,i)
    Gsv_LOO = LeaveOneOut(Gsv,i)
    
    J_LOO = get_J(Ec_LOO, Gsv_LOO, Gcv)[0]
    Ec_predict[i] = np.dot(pi[i],J_LOO)
    
CV = np.sum(np.power((Ec_predict - Ec),2))/ns    

    