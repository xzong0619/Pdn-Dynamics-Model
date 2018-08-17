# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 11:08:25 2018

@author: wangyf
"""

from config import mother
import numpy as np
import lattice_functions as lf
from itertools import combinations

def layer_tuple(mother, ci):
    
    n = len(ci)
    index = []
    for i in range(n):
        index.append(ci[i])
        
    layer_tuple = []
    
    for i in range(n):
        layer_tuple.append(mother[index[i]][2])
        
    layer_tuple = tuple(layer_tuple)
        
    return layer_tuple
    
def  distance_tuple(mother, ci):
    
    n = len(ci)
    index = []
    
    for i in range(n):
        index.append(ci[i])
    
    combo = list(combinations(index,2))
    ncombo = len(combo) #0 for 1 node, 2 for 2 nodes, 3 for 3 nodes
      
    distance_tuple = []
    
    for i in range(ncombo):
        pt1 = mother[combo[i][0]]
        pt2 = mother[combo[i][1]]
        distance_tuple.append(lf.two_points_D(pt1, pt2))
        
    distance_tuple = tuple(sorted(distance_tuple))
        
    return distance_tuple
    
    
index = np.arange(len(mother))

c1 = list(combinations(index,1))

c2 = list(combinations(index,2))

c3 = list(combinations(index,3))

'''
use c2 as an example
'''
d2 = []
for i in range(len(c2)):
    a = c2[i][0]
    b = c2[i][1]
    d2.append((lf.two_points_D(mother[a], mother[b]), (mother[a][2], mother[b][2])))
d2d = list(set(d2))

index_list =[]
for i in d2d:
    index_list.append(d2.index(i))
index_list.sort()

d2_d = np.array(c2)[index_list]


'''
Now lets try c3
'''
d3 = []
for i in range(len(c3)):
    a = c3[i][0]
    b = c3[i][1]
    c = c3[i][2]
    t1 = tuple(sorted((lf.two_points_D(mother[a], mother[b]),
                       lf.two_points_D(mother[b], mother[c]),
                       lf.two_points_D(mother[a], mother[c]))))
    d3.append((t1, (mother[a][2], mother[b][2], mother[c][2])))
                        
d3d = list(set(d3))    
    
index_list3 =[]
for i in d3d:
    index_list3.append(d3.index(i))
#index_list3.sort()

d3_d = np.array(c3)[index_list3]    
    
#%%
'''
Now lets try oop
'''

from class_subgraph import *

ob = subgraphs(mother)

x = ob.get_s(1)









    