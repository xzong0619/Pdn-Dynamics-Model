#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 15:23:51 2018

@author: wangyifan
"""

class Neighboring_list(Input):
    
    '''
    #Generating neighboring list for both mechanism_input and enegertics_input
    
    '''
    def __init__(self,  neighboring_list = None):
        
        super().__init__()
        
        if neighboring_list == None:
            self.neighboring_list = []
            for i in range(self.n_surf):
                self.neighboring_list.append('1-2')
        else: 
            self.neighboring_list = neighboring_list
            

        self.n_rxn = len(self.rxn) # number of surface reactions
        
        self.type = [] # surface reaction types
        for i in range(self.n_rxn):
            self.type.append(self.rxn[i].type)

        self.ai= [] # index of a in the neighboring list
        for i in range(self.n_rxn):
            self.ai.append(self.rxn[i].ai)
        
        self.bi= []
        for i in range(self.n_rxn):
            if self.rxn[i].bi == []:
                self.bi.append(None)
            else:
                self.bi.append(self.rxn[i].bi)
        
        self.ci= []
        for i in range(self.n_rxn):
            if self.rxn[i].ci == []:
                self.ci.append(None)
            else:
                self.ci.append(self.rxn[i].ci)
        
        '''
        #Initialize the neighboring list with 1+(n-1) -> n
        '''
        
        nc_edge_list =  []
        old_wt_i = []

                
                
                
                
 
        
        for r in range(self.n_rxn):
            if self.type[r] == 'C' 
                ai = self.ai[r] - 1
                bi = self.bi[r] - 1 
                ci = self.ci[r] - 1