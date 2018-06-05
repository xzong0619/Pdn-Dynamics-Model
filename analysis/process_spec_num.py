# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 14:51:12 2018

@author: wangyf
"""

import os
import numpy as np
import matplotlib.pyplot as plt

class read_Sim:
    
    fname = 'simulation_input.dat'
    
    def __init__(self, fldr = None):
        
        if fldr == None:
            self.fldr = os.getcwd()
        else:
            self.fldr = fldr
    
        
        filepath = os.path.join(self.fldr, self.fname)
        fid = open(filepath, 'r')
        file = fid.read()
        self.line = file.splitlines()
        self.lines  = []
        
        for nline in range(len(self.line)):

            self.lines.append(self.line[nline].split())
        
        for nline in range(len(self.lines)):
            
            if not self.lines[nline] == []: 
                
                txt_contain = self.lines[nline] 
                
                if txt_contain[0] == 'temperature': self.T = float(txt_contain[1])
                if txt_contain[0] == 'pressure': self.P = float(txt_contain[1])
                if txt_contain[0] == 'n_surf_species': self.n_spec = int(txt_contain[1])
                if txt_contain[0] == 'surf_specs_names': self.surf_spec_name = txt_contain[1:]
                    
                    
class spec_info:
    
    def __init__(self,  data,  dict):
        
        '''
        Fetch data from simulation input
        '''
        
        Sim = read_Sim()
        self.n_spec = Sim.n_spec
        self.surf_spec_name = Sim.surf_spec_name
        self.surf_spec = {}
        for i in range(self.n_spec):
            self.surf_spec[self.surf_spec_name[i]] = int(data[dict[self.surf_spec_name[i].lower()]])
        '''
        Fill reactions with species name and species data
        '''
        self.entry = int(data[dict['entry']])        
        self.nevents = int(data[dict['nevents']])     
        self.t = float(data[dict['time']])
        self.T = float(data[dict['temperature']])
        self.E = float(data[dict['energy']])
        '''
        self.spec = {}
        self.Pd1 = int(data[dict['pd1*']])
        self.Pd2 = int(data[dict['pd2*']])
        self.Pd3 = int(data[dict['pd3*']])
        self.Pd4 = int(data[dict['pd4*']])
        self.Pd5 = int(data[dict['pd5*']])
        '''
        self.CO = int(data[dict['co']])                        
    
class read_Single_Spec:
    
    fname =  'specnum_output.txt'
    
    def __init__(self, fldr = None):
        
        if fldr == None:
            self.fldr = os.getcwd()
        else:
            self.fldr = fldr
    
        
        self.filepath = os.path.join(self.fldr, self.fname)

        fid = open(self.filepath, 'r')
        file = fid.read()
        lines = file.splitlines()
        dict_array = lines[0].lower().split()
        
        data = lines[1:]
        dict = {}
        for x in range(0, len(dict_array)):
            dict[dict_array[x]] = x
        
        spec = []
        for s in data:
            spec.append(spec_info(s.split(),dict))
        
        n = len(spec)
        self.n = n       
        self.t = np.zeros(n)
        self.n_spec = spec[0].n_spec
        self.surf_spec_name = spec[0].surf_spec_name
        
        self.surf_spec_dic = {}
        
        for i in range(self.n_spec):      
            self.surf_spec_dic[self.surf_spec_name[i]] = np.zeros(n)
        
        
        for i in range(n):
            self.t[i] = spec[i].t
            
            for j in range(self.n_spec):
                
                self.surf_spec_dic.[self.surf_spec_name[i]][j] = spec[i].[self.surf_spec_name[i]][j]
                
            
            
        '''
    def num_to_cov(self, lattice_dim):
        
        sites_per_unitcell = 4
        total_sites = lattice_dim**2 *sites_per_unitcell
        # in percentage coverage
        self.Pd1 = self.Pd1/total_sites *100
        self.Pd2 = self.Pd2/total_sites *100 *2
        self.Pd3 = self.Pd3/total_sites *100 *3
        self.Pd4 = self.Pd4/total_sites *100 *4
        self.Pd5 = self.Pd5/total_sites *100 *5
x = read_Sim()
y = read_Single_Spec()