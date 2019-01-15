# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 11:27:10 2018

@author: wangyf
"""

#test the reverse problem of pi matrix 

import lattice_functions as lf
import numpy as np
from itertools import combinations, product


l1 = np.array([(0, 0), (1, 0), (1 / 2, 3**0.5 / 2)])
l2 = np.array([np.sum(l1[[0, 1, 2]], 0)/3])
dz = 1 
l1d = lf.add_z(l1, dz)
l2d = lf.add_z(l2, dz*2)

mother = np.concatenate((l1d, l2d), axis=0)

config = [[0,1]]
Ec = [0]
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
draw = [1, 1, 0]


Clusters = lf.clusters(occ, NN1, draw)
Clusters.get_mother(mother, dz)
Gm = Clusters.Gm
Clusters.get_configs(config)
Gsv = Clusters.Gsv


sub = lf.subgraphs(mother, dz)
Gcv1 = sub.get_s2(1)
Gcv2 = sub.get_s2(2)
Gcv3 = sub.get_s2(3)

Gcv = Gcv1+Gcv2+Gcv3
niso = []
for Gi in Gcv:
    niso.append(len(Gi))        
niso = np.array(niso)   
    
Cal = lf.calculations(occ)

#delta =  Cal.get_delta_l(Gsv[0] ,Gcv[0])  
pi_vector = Cal.get_pi_matrix(Gsv, Gcv)

'''
Now let's do reverse engineering! 
'''
# num is number of clusters in the supergraph
num = np.multiply(pi_vector,niso)[0]
num = num.astype(int)



def get_com (L,n):
    
    com = tuple(combinations(L,n))
    ncom = len(com)
    
    return ncom,com

def get_all_com(Gcv, num):   

    ncom_list = []
    com_list = []
    
    for gi,ng in enumerate(num):
        
        if not ng==0:
            ncom,com =  get_com(Gcv[gi], ng)
            ncom_list.append(ncom)
            com_list.append(com)
            
    return ncom_list, com_list    

def possible_config(ncom_list, com_list):    
        
    n_possible = np.product(ncom_list)
    nG = len(ncom_list)
    
    combo_list = com_list[0]
    for gi in range(1,nG):
        combo_list = product(combo_list, com_list[gi])
        
     
    combo_list = list(combo_list)
        
    return n_possible,combo_list

def combine_list(ll):

    lnew  = ()
    lnodes = []
    for i in range(len(ll)):
        
        if isinstance(ll[i],tuple):
            lnew= lnew + ll[i]
        else:
            lnodes.append(ll[i])
    
    return lnew,lnodes     

def get_nodes(ll):
    x = ll
    y = []
    while len(x) > 0:
        x,ynew = combine_list(x)
        y = y+ynew
    y = tuple(set(y))
    
    return y

def get_unique_graphs(combo_list):

    y = []
    for gi in combo_list:
        y.append(get_nodes(gi))
    
    y = list(set(y))  
    
    return y
          



ncom_list, com_list = get_all_com(Gcv, num)
n_possible,combo_list = possible_config(ncom_list, com_list)
gs = get_unique_graphs(combo_list)

 