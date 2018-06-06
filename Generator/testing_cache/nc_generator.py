# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:58:04 2018

@author: wangyf
"""


from IO_Pdn import *
import Cluster_structure as cs
import networkx as nx


info  = Input()
names = info.surf_spec
dent = info.surf_dent

nc_list = [None,  # neighboring structure
      '1-2',
      '1-2 1-3 2-3',
      '1-2 1-3 2-3 4-2 4-3 ',
      '1-2 1-3 2-3',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 ',
      '1-2 1-3 2-3 4-2 4-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-4 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-4 6-2 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 ',
      '1-2 1-3 2-3 4-2 4-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 7-6 7-3 1-7 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-4 7-6 7-5 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 7-6 7-5 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-4 6-2 7-6 7-2 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-4 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-4 6-2 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 ']

class NCIn():
    
    def __init__(self, str1):
        
        self.l= [str1] 
        self.n = len(self.l)
        
    def addnc(self,str2):
        self.l.append(str2)
        self.n = len(self.l)
        

n_rxn = len(info.rxn) # number of surface reactions
        
types = [] # surface reaction types
for i in range(n_rxn):
    types.append(info.rxn[i].type)

ai= [] # index of a in the neighboring list
for i in range(n_rxn):
    ai.append(info.rxn[i].ai)

bi= []
for i in range(n_rxn):
    if info.rxn[i].bi == []:
        bi.append(None)
    else:
        bi.append(info.rxn[i].bi)

ci= []
for i in range(n_rxn):
    if info.rxn[i].ci == []:
        ci.append(None)
    else:
        ci.append(info.rxn[i].ci)

add_list = []


nc  = []
# Initializing nc object
for i in range(info.n_surf):
        nc.append(NCIn(nc_list[i]))
        
        
def nc_generator(r):


    ni_u = 1
    more_structure_flag  = 0
    new_H_str = None
    
    if types[r] == 'O' or types[r] == 'C': 
        ai_py = ai[r] - 1
        bi_py = bi[r] - 1 
        ci_py = ci[r] - 1
        
        if dent[ci_py] > 5:
    
                bnc = cs.neighbor_string_list([nc_list[bi_py]])
                cnc = cs.neighbor_string_list( [nc_list[ci_py]])
                
                old_edge = bnc
                old_wt_i = cs.neighbor_weight(bnc)
                
                '''
                n = 1 + dent[bi_py]
    
                Pdx = cs.Pdn(n, old_edge, old_wt_i)
                Pdx.new_graph()
                ni_u = Pdx.ni_u
                new_H = Pdx.H
                new_wt_i = Pdx.wt_i
                '''
                new_H  =  [old_edge]
                new_wt_i = [old_wt_i]
                n = dent[bi_py]
                for i in range(dent[ai_py]):
                    
                    n = n+1 
                    Pdx = cs.Parallel_growth(n, ni_u, new_H, new_wt_i)
                    Pdx.unique()
                    ni_u = Pdx.mi_u
                    new_H = Pdx.H_u
                    new_wt_i = Pdx.wt_i_u
                    
                    
                                 
                for i in range(ni_u):
                    A = nx.Graph(new_H[i])
                    B = nx.Graph(cnc)
                    if nx.is_isomorphic(A,B):
                        if not set(new_H[i]) == set(cnc):
                            more_structure_flag  = 1
                            
                            new_H_str = cs.neighbor_list_string(new_H[i])
                            nc[ci_py].addnc(new_H_str) 
        
    return more_structure_flag, new_H_str                    

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
    
                        
for r in range(n_rxn):
    [more_structure_flag, new_H_str] = nc_generator(r)  
    if more_structure_flag:           
        add_list.append(new_H_str)   

for c in range(info.n_surf):
    if nc[c].n > 1:
        H_list = []
        for cn in range(nc[c].n):
            H_list.append(cs.neighbor_string_list([nc[c].l[cn]]))
        i_list = cs.clean_up(H_list)
        ni = len(i_list)
        nc[c].l = []
        for i in range(ni):
            nc[c].addnc(H_list[i_list[i]])