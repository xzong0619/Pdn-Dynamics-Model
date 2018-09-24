# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:07:04 2018

@author: wangyf
"""
import lattice_functions as lf
import numpy as np
from itertools import combinations, product

class reverse():
    
    '''
    Create the reverse graph 
    Input: one configuration grpah Gsv
           the cluster index list Gcv
           the occupancy colors 
    '''
    def __init__(self, occ, Gcv):
        
        self.Gcv = Gcv
        self.occ = occ
        
        niso = []
        for Gi in self.Gcv:
            niso.append(len(Gi)) 
        
        #niso is the vector of number of isomorphic in Gcv
        self.niso = np.array(niso)  
                
    def get_supergraph(self, pi_vector):
        
        # num is number of clusters in the supergraph
        ncluster = np.multiply(pi_vector,self.niso)[0]
        self.ncluster = ncluster.astype(int)
        
        '''
        ncom_list, com_list = self.get_all_com(self.Gcv, self.ncluster)
        n_possible,combo_list = self.possible_config(ncom_list, com_list)
        self.gsuper = self.get_unique_graphs(combo_list)
        '''
    
    @staticmethod
    def get_com(L,n):
        
        com = tuple(combinations(L,n))
        ncom = len(com)
        
        return ncom,com
    

    def get_all_com(self, Gcv, num):   
    
        ncom_list = []
        com_list = []
        
        for gi,ng in enumerate(num):
            
            if not ng==0:
                
                ncom,com =  self.get_com(Gcv[gi], ng)
                ncom_list.append(ncom)
                com_list.append(com)
                
        return ncom_list, com_list    
    
    @staticmethod
    def possible_config(ncom_list, com_list):    
            
        n_possible = np.product(ncom_list)
        nG = len(ncom_list)
        
        combo_list = com_list[0]
        for gi in range(1,nG):
            combo_list = product(combo_list, com_list[gi])
            
        combo_list = list(combo_list)
            
        return n_possible,combo_list
    
    @staticmethod
    def combine_list(ll):
    
        lnew  = ()
        lnodes = []
        for i in range(len(ll)):
            
            if isinstance(ll[i],tuple):
                lnew= lnew + ll[i]
            else:
                lnodes.append(ll[i])
        
        return lnew,lnodes     
    
    
    def get_nodes(self,ll):
        x = ll
        y = []
        while len(x) > 0:
            x,ynew = self.combine_list(x)
            y = y+ynew
        y = tuple(set(y))
        
        return y
    
    def get_unique_graphs(self,combo_list):
    
        y = []
        for gi in combo_list:
            y.append(self.get_nodes(gi))
        
        y = list(set(y))  
        
        return y
   
    def test(self, Gsv):
        
        Cal = lf.calculations(self.occ)
        #pi_vector is the pi vector for this Gsv
        pi_vector = Cal.get_pi_matrix(Gsv, self.Gcv)
        
        self.get_supergraph(pi_vector)
        
        