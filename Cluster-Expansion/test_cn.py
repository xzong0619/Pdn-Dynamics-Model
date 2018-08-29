# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:25:35 2018

@author: wangyf
"""
'''
Use No.26 Pd10a to test
'''
from structure_constants import *
import lattice_functions as lf
import pickle
import numpy as np

def num_1NN(G, i):
    
    '''
    G is a networkx graph
    i is the index of node in G
    
    returns number of nearest neighbors of node i 
    and a list of neighbor index number
    '''
    n_1NN = 0   
    list_1NN = []
    if G.nodes[i]['color'] == filled:  # check if the node is occupied 
        for j in list(G.neighbors(i)): #iterate through 1st NN
            if G.nodes[j]['color'] == filled: #check if the node is occupied
                n_1NN = n_1NN + 1
                list_1NN.append(j)
    else:
        print('No atoms detected at this position') 
    return n_1NN, list_1NN

def num_2NN(G,i):
    
    '''
    G is a networkx graph
    i is the index of node in G where CO adsorbs
    
    returns a list of numbers of 2nd nearest neighbors of node i for each 1NN
    and a 2D list of 2NN index numbers
    '''
    n_2NN = []
    list_2NN = []
    if G.nodes[i]['color'] == filled: # check if the node is occupied 
        for j in G.neighbors(i):      # iterate through 1st NN
            if G.nodes[j]['color'] == filled: # check if the node is occupied 
                n_2NN.append(num_1NN(G,j)[0]) # Add number of 2NN for 1NNs
                list_2NN.append(num_1NN(G,j)[1]) # Add neighbor index number
    else:
        print('No atoms detected at this position')            
    return n_2NN, list_2NN

def cal_CN1(G,COsites):
    
    '''
    G is a networkx graph
    COsites are the index of adsorption site
    
    returns 1st CN number
    '''
    
    CN1 = []
    sitetype = len(COsites)
    for i in range(sitetype):
        CN1.append(num_1NN(G, COsites[i])[0])
        
    '''
    take arithmaric mean for CN1 for bridge and hollow sites
    '''
    
    CN1 = np.mean(np.array(CN1))   
    
    return CN1


def cal_CN2(G,COsites):
       
    '''
    G is a networkx graph
    COsites are the index of adsorption site
    
    returns 2nd CN number
    '''
    
    
    list_CN2 = []
    CN2 = []
    
    sitetype = len(COsites)
    
    for i in range(sitetype):
        list_CN2.append(num_2NN(G, COsites[i])[0])
        
    '''
    sum up 2NN numbers for each 1NN
    '''
    
    for i in range(sitetype):
        if len(list_CN2[i]) == 0: CN2.append(0)
        else: CN2.append(np.sum(np.array(list_CN2[i])))
    
    '''
    take arithmaric mean for CN2 for bridge and hollow sites
    '''
    
    CN2 = np.mean(np.array(CN2)) 
    
    return CN2

def cal_GCN(G, COsites):
    
    '''
    G is a networkx graph
    COsites are the index of adsorption site
    
    returns general coordination number
    '''
    
    GCN = []
    
    sitetype = len(COsites)
    list_1NN = []
    '''
    find all avaiable 1NN index 
    '''
    for i in range(sitetype):
        list_1NN = list_1NN + num_1NN(G, COsites[i])[1]
    
    '''
    Use set to avoid double counting
    '''    
    list_1NN = list(set(list_1NN))
    '''
    Get CN for these 1NN nodes
    '''
    for i in list_1NN:
        GCN.append(num_1NN(G,i)[0])

    '''
    Set weight based on Pd(111)
    '''
    
    if len(COsites) == 1: weight = 12
    if len(COsites) == 2: weight = 18
    if len(COsites) == 3: weight = 22
    
    GCN = np.sum(np.array(GCN))/weight

    return GCN


def num_Ce1NN(G, i):
    
    '''
    G is a networkx graph
    i is the index of node in G
    
    returns the flag of whether the atom is next to a Ce atom
    '''
    
    nCe = 0 
    if G.nodes[i]['color'] == filled: # check if the node is occupied
        if G.nodes[i]['z'] == '1':  # check if the node is in base layer
            nCe = 1 # 1 means the atom is in contact with 3 Ce atoms underneath
    return nCe
        
def num_Ce2NN(G, i):
    
    '''
    G is a networkx graph
    i is the index of node in G
    
    returns the number of 2nd nearest Ce neighbors of node i 
    and a list of atom adjacent to Ce base layer
    '''

    n_2NN = 0   
    list_2NN = []
    if G.nodes[i]['color'] == filled: # check if the node is occupied
        for j in list(G.neighbors(i)): # iterate through 1st NN
            n_2NN = n_2NN + num_Ce1NN(G,j)  # check if the node is next to Ce
            if num_Ce1NN(G,j): list_2NN.append(j) # Append the index to a list
    else:
        print('No atoms detected at this position') 
    return n_2NN, list_2NN

def cal_CeCN1(G, COsites):
    
    '''
    G is a networkx graph
    COsites are the index of adsorption site
    
    returns coordination number of Ce
    '''
    
    CN1 = []
    sitetype = len(COsites)
    
    for i in range(sitetype):
        CN1.append(num_Ce1NN(G, COsites[i]))
    
    '''
    take arithmaric mean for CN2 for bridge and hollow sites
    '''
    
    CN1 = np.mean(np.array(CN1)) * 3 # each Pd atom is coordinate by 3 Ce
    
    return CN1
    
def cal_CeCN2(G, COsites):
    
    '''
    G is a networkx graph
    COsites are the index of adsorption site
    
    returns 2nd coordination number of Ce
    '''
    
    CN2 = []
    sitetype = len(COsites)
    
    for i in range(sitetype):
        CN2.append(num_Ce2NN(G, COsites[i])[0])
    
    '''
    take arithmaric mean for CN2 for bridge and hollow sites
    '''
    
    CN2 = np.mean(np.array(CN2)) *3 

    
    return CN2

def cal_CeGCN(G, COsites):

    '''
    G is a networkx graph
    COsites are the index of adsorption site
    
    returns general coordination number of Ce
    '''
    
    GCN = []
    
    sitetype = len(COsites)
    list_1NN = []
    
    '''
    find all avaiable 1NN next to Ce index 
    '''
    
    for i in range(sitetype):
        list_1NN = list_1NN + num_Ce2NN(G, COsites[i])[1]
        
    list_1NN = list(set(list_1NN))
    
    '''
    Check if Ce is around for 1NN nodes
    '''
    
    for i in list_1NN:
        GCN.append(num_Ce1NN(G,i))
    
    
    if len(COsites) == 1: weight = 3
    if len(COsites) == 2: weight = 5
    if len(COsites) == 3: weight = 6
    
    GCN = np.sum(np.array(GCN))/weight * 3

    return GCN
            

empty = 'grey'
filled = 'r'
occ = [empty, filled]  
    
[Gm, Gsv, Gcv1, Gcv2, Gcv3] = pickle.load(open("clusters_NN1.p", "rb"))

G = Gsv[25] #No.26 Pd10a [0,1,2,6,11,10,14,29,28,27]
COsites =  [2, 14, 27]

CN1 = cal_CN1(G,COsites)
CN2 = cal_CN2(G,COsites)
GCN = cal_GCN(G,COsites)

CeCN1 = cal_CeCN1(G,COsites)
CeCN2 = cal_CeCN2(G,COsites)
CeGCN = cal_CeGCN(G,COsites)
        
