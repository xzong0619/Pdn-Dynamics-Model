# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:52:08 2018

@author: wangyf
"""
from structure_constants import mother, dz, config, Ec
import lattice_functions as lf
import numpy as np
import pandas as pd
import pickle
from itertools import combinations, product



node_layer_v = np.around(mother[:,2]/dz, decimals = 0).astype(int)
node_layer = dict()
n_layer = np.amax(node_layer_v)

for layer_i in range(n_layer):
    node_layer[layer_i +1] = list(np.where(node_layer_v == layer_i +1)[0] + 1)   
    
config_i = config[6]
config_layers = np.around(mother[np.array(config_i)][:,2]/dz, decimals = 0).astype(int)
layer_count = []

for layer_i in range(n_layer):  
    layer_count.append(len(np.where(config_layers == layer_i + 1)[0]))

layer_count = list(np.array(layer_count)[np.nonzero(layer_count)[0]])

combo_dict = dict()

for i,layer_count_i in enumerate(layer_count):
    layer_i  = i + 1
    combo_dict[layer_i] =  list(combinations(node_layer[layer_i], layer_count_i))

combo = combo_dict[1]
for i in range(1,len(layer_count)):
    combo = product(combo, combo_dict[i+1])
combo = list(combo)   
for i,combo_i in enumerate(combo):
    temp = ()
    for j, combo_ij in enumerate(combo_i):
        temp = temp + combo_i[j]
    combo[i] = temp
#    
#empty = 'grey'
#filled = 'r'
#occ = [empty, filled]
#
#'''
#only draw 1st nearest neighbors?
#'''
#NN1 = 0
#'''
#Draw mother/conifgurations/clusters?
#'''
#draw = [1, 0, 0]
#
#
#Clusters = lf.clusters(occ, NN1, draw)
#Clusters.get_mother(mother, dz)
#Gm = Clusters.Gm
#
#
##%%
#'''
#Create Configurations
#'''
#Clusters.get_configs(config)
#Gsv = Clusters.Gsv
#
#
##%%
#'''
#Create clusters
#'''
#sub = lf.subgraphs(mother, dz)
#Gcv1 = sub.get_s2(1)
#Gcv2 = sub.get_s2(2)
#Gcv3 = sub.get_s2(3)
#Gcv4 = sub.get_s2(4)
#
##pickle.dump([Gm, Gsv, Gcv1, Gcv2, Gcv3], open('clusters.p','wb'))
