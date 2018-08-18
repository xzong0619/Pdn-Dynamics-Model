# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 11:08:25 2018

@author: wangyf
"""

from structure_constants import mother, config, Ec
import numpy as np
import lattice_functions as lf
from itertools import combinations
from config import Gsv
import pickle


empty = 'grey'
filled = 'r'
occ = [empty, filled]

def unique_combo(combo, indices_list):
    
    Gcv = []
    nclusters = len(indices_list)
    
    for i in range(nclusters):
        Gcv.append([])
        niso = len(indices_list[i])
        for j in range(niso):
            Gcv[i].append(combo[indices_list[i][j]])
    
    return Gcv


def layer_tuple(mother, ci):
    
    n = len(ci)
    index = []
    for i in range(n):
        index.append(ci[i])
        
    layer = []
    
    for i in range(n):
        layer.append(int(mother[index[i]][2]))
        
    layer = tuple(layer)
        
    return layer
    
def distance_tuple(mother, ci):
    
    n = len(ci)
    index = []
    
    for i in range(n):
        index.append(ci[i])
    
    combo = list(combinations(index,2))
    ncombo = len(combo) #0 for 1 node, 2 for 2 nodes, 3 for 3 nodes
      
    distance= []
    
    for i in range(ncombo):
        pt1 = mother[combo[i][0]]
        pt2 = mother[combo[i][1]]
        distance.append(lf.two_points_D(pt1, pt2))
        
    distance = tuple(sorted(distance))
        
    return distance
    
def np_to_list(s_np):
    
    s_list = []
    for i in range(s_np.shape[0]):
        s_list.append(list(s_np[i]))
    return s_list
            
def get_occupancy(G, i):
    
    '''
    Get the occupancy from the graph G for node i 
    Occupied is 1 and unoccupied is 0
    '''
    
    if G.nodes[i]['color'] == empty: o = 0
    if G.nodes[i]['color'] == filled: o = 1 
    
    return o    
    

def get_delta(Gl, Gs):
    
    '''
    takes in larger graph Gl and smaller graph Gs
    find sub isomorphic graphs of Gs from Gl
    calculate the delta value in pi matrix 
    '''
    '''
    if there are more than 2 nodes in a cluster
    '''
    niso =len(Gs)
    
    '''
    save product into a list 
    and caclulate the sum divid by total number of subgraphs
    '''
    subi = []
    subs = []
    for i in range(niso):
        subi.append([])
        for j in range(len(Gs[i])):
            subi[i].append(get_occupancy(Gl,Gs[i][j]))   
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
    progress = 0
    
    for i in range(n1):
        for j in range(n2):
            pi[i][j] = get_delta(G1v[i],G2v[j])
            
            progress = progress + 1
            per = progress/n1/n2 *100
            print('%.2f %% done!' %per)
            
            
            
    J = np.linalg.lstsq(pi, Ev)[0]
    
    J = J
    pi = pi
    
    return J, pi
    


index = np.arange(len(mother))

c1 = list(combinations(index,1))

c2 = list(combinations(index,2))

c3 = list(combinations(index,3))

'''
add c1 
'''
sub = lf.subgraphs(mother)
d1 = sub.get_s(1)




'''
use c2 as an example
'''
d2 = []

for i in range(len(c2)):

    ci = c2[i]
    
    distances = distance_tuple(mother, ci)
    layers = layer_tuple(mother, ci)
    
    d2.append((distances, layers))
    
d2d = list(set(d2))

index_list =[]
indices_list = []
for i in d2d:
    index_list.append(d2.index(i))
index_list.sort()

for i in index_list:
    indices_list.append([a for a, x in enumerate(d2) if x == d2[i]])

d2_d = np.array(c2)[index_list]
d2_d = np_to_list(d2_d)

'''
creat pi matrix
size of number of configuration * numbers of clusters
'''
gcv2 = unique_combo(c2, indices_list)
 

J2, pi2 = get_J(Ec, Gsv ,gcv2 ) 
pi1 = np.load('pi2.npy')
'''
Now lets try c3
'''

d3 = []
for i in range(len(c3)):
    
    ci = c3[i]
    
    distances = distance_tuple(mother, ci)
    layers = layer_tuple(mother, ci)
 
    d3.append((distances, layers))
                        
d3d = list(set(d3))    
    
index_list3 =[]
indices_list3 = []
for i in d3d:
    index_list3.append(d3.index(i))
index_list3.sort()

for i in index_list3:
    indices_list3.append([a for a, x in enumerate(d3) if x == d3[i]])

d3_d = np.array(c3)[index_list3]    
d3_d = np_to_list(d3_d)    








    