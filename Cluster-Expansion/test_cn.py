# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:25:35 2018

@author: wangyf
"""
'''
Use No.26 Pd10a to test
'''
from structure_constants import mother, config, Ec
import lattice_functions as lf
import pickle
import numpy as np

def num_1NN(G, i):
    
    '''
    G is a networkx graph
    i is the index of node in G
    returns number of nearest neighbors of node i
    '''
    n_1NN = 0   
    list_1NN = []
    if G.nodes[i]['color'] == filled: 
        for j in list(G.neighbors(i)):
            if G.nodes[j]['color'] == filled:
                n_1NN = n_1NN + 1
                list_1NN.append(j)
    else:
        print('No atoms detected at this position') 
    return n_1NN, list_1NN

def num_2NN(G,i):
    
    '''
    G is a networkx graph
    i is the index of node in G where CO adsorbs
    returns a list of numbers of 2nd nearest neighbors of node i
    Avoid double counting
    '''
    n_2NN = []
    list_2NN = []
    if G.nodes[i]['color'] == filled:   
        for j in G.neighbors(i):
            if G.nodes[j]['color'] == filled:
                n_2NN.append(num_1NN(G,j)[0])
                list_2NN.append(num_1NN(G,j)[1])
    else:
        print('No atoms detected at this position')            
    return n_2NN, list_2NN

def cal_CN1(G,COsites):
    
    CN1 = []
    sitetype = len(COsites)
    for i in range(sitetype):
        CN1.append(num_1NN(G, COsites[i])[0])
    CN1 = np.mean(np.array(CN1))  
    
    return CN1

def cal_CN2(G,COsites):
    
    list_CN2 = []
    CN2 = []
    
    sitetype = len(COsites)
    
    for i in range(sitetype):
        list_CN2.append(num_2NN(G, COsites[i])[0])

    
    for i in range(sitetype):
        if len(list_CN2[i]) == 0: CN2.append(0)
        else: CN2.append(np.mean(np.array(list_CN2[i])))
    
    CN2 = np.mean(np.array(CN2)) 
    
    return CN2

def cal_GCN(G, COsites):
    
    GCN = []
    
    sitetype = len(COsites)
    list_1NN = []
    
    for i in range(sitetype):
        list_1NN = list_1NN + num_1NN(G, COsites[i])[1]
        
    list_1NN = list(set(list_1NN))
    
    for i in list_1NN:
        GCN.append(num_1NN(G,i)[0])
        
    if len(COsites) == 1: weight = 12
    if len(COsites) == 2: weight = 18
    if len(COsites) == 3: weight = 22
    
    GCN = np.sum(np.array(GCN))/weight

    return GCN



empty = 'grey'
filled = 'r'
occ = [empty, filled]  
    
[Gm, Gsv, Gcv1, Gcv2, Gcv3] = pickle.load(open("clusters_NN1.p", "rb"))

G = Gsv[25] #No.26 Pd10a [0,1,2,6,11,10,14,29,28,27]
COsites =  [2, 14, 27]

sitetype = len(COsites)
CN1 = cal_CN1(G,COsites)
CN2 = cal_CN2(G,COsites)
GCN = cal_GCN(G,COsites)

        
#CN1 = np.mean(np.array(CN1))    