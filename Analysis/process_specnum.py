# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 14:51:12 2018

@author: wangyf
"""

import os
import sys
import numpy as np
import pandas as pd

sys.path.append(r"C:\Users\wangyf\Documents\GitHub\Pdn-Dynamics-Model\Generator")



def final_spec_on_surf(single_spec_df):
    
    ''' 
    takes in a single species dataframe
    return a list of available species on the surface at final state
    '''
    
    end_state = single_spec_df[-1:]
    
    
    end_state_v = np.array(end_state)[0]
    surf_spec_names = end_state.columns 
    n_surf  = len(surf_spec_names) -1  #-1 to omit t column 
    s_name = [] # specify the initial surface species
    s_n = [] # specify the corresponding number on the surface 
    s_count = 0 # number of types of initial surface species
    
    for i in range(n_surf):
        if not end_state_v[i] == 0:
            
            s_name.append(surf_spec_names[i])
            s_n.append((end_state_v[i]))
            s_count = s_count +1
    
    return s_name, s_n, s_count        
            
def lifetime_spec_on_surf(single_spec_df):

    ''' 
    takes in a single species dataframe
    return a list of available species on the surface at final state
    '''
    
    surf_spec_names = single_spec_df.columns
    n_surf  = len(surf_spec_names) -1 #-1 to omit t column 
    s_name = [] # specify the initial surface species
    s_count = 0 # number of types of initial surface species
    
    for i in range(n_surf):
        spec_i = np.array(single_spec_df[surf_spec_names[i]])
        
        if np.any(spec_i):
            
            s_name.append(surf_spec_names[i])
            s_count = s_count +1
    
    return s_name, s_count        

def final_lt_spec_on_surf(single_spec_df):
    
    s_name, s_count  = lifetime_spec_on_surf(single_spec_df)
    end_state = single_spec_df[-1:]
    end_state = end_state[s_name]
    s_n = np.array(end_state)[0]
    
    return s_n

#%%

class restart:
    
    fname = 'specnum_output.txt'
    
    '''
    
    Input directory needed (where spectrum.dat file located)
    
    '''
    
    
    def __init__(self, input_dir):
                        
        single_spec_df =  read_Single_Spec(input_dir).surf_spec_df
        
        self.s_name, self.s_n, self.s_count  = final_spec_on_surf(single_spec_df)
        
        
#%%
class read_Sim:
    
    fname = 'simulation_input.dat'
    input_dir = 'zacros_inputs'
    
    def __init__(self, fldr = None):
        
        if fldr == None:
            self.fldr =  os.path.join(os.getcwd(), self.input_dir)
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
        self.CO = int(data[dict['co']])                        
    
class read_Single_Spec:
    
    fname =  'specnum_output.txt'
    input_dir = 'zacros_inputs'
    
    def __init__(self, fldr = None):
        
        if fldr == None:
            self.fldr =  os.path.join(os.getcwd(), self.input_dir)
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
        self.surf_spec_dic['t'] = np.zeros(n)
        for i in range(self.n_spec):      
            self.surf_spec_dic[self.surf_spec_name[i]] = np.zeros(n)
        
        
        for i in range(n):
            
            self.surf_spec_dic['t'][i] = spec[i].t
            
            for j in range(self.n_spec):
            
                self.surf_spec_dic[self.surf_spec_name[j]][i] = spec[i].surf_spec[self.surf_spec_name[j]]
        
        self.surf_spec_df = pd.DataFrame.from_dict(self.surf_spec_dic)
        self.lifetime_spec  = lifetime_spec_on_surf(self.surf_spec_df)
        
                
class read_Multiple_Spec:
    
    input_dir = 'zacros_inputs'     
    
    def __init__(self, n_files, fldr = None):
         
        if fldr == None:
            self.fldr =  os.path.join(os.getcwd(), self.input_dir)
        else:
            self.fldr = fldr
        
        self.filepath = []
        
        for f in range(n_files):
             self.filepath.append(os.path.join(self.fldr, 'outputs', str(f+1)))
        
        single_spec =  read_Single_Spec(self.filepath[0])
        single_spec_df = single_spec.surf_spec_df
        single_spec_n = [single_spec.n]
        
        
        
        for f in range(n_files-1): 
            
            single_spec =  read_Single_Spec(self.filepath[f+1])
            
            single_spec_df = single_spec_df + single_spec.surf_spec_df
            
            single_spec_n.append(single_spec.n)
            
        self.n = min(single_spec_n)
        
        self.multi_spec_ave_df = single_spec_df[:self.n]/n_files
        self.t = self.multi_spec_ave_df['t']
        self.lifetime_spec  = lifetime_spec_on_surf(self.multi_spec_ave_df)
        #del self.multi_spec_ave_df['t']
            

#%%
def num_to_cov(lattice_dim, n_spec, surf_dent_vec, spec_vecs):
        
    sites_per_unitcell = 4
    total_sites = lattice_dim**2 *sites_per_unitcell
    # in percentage coverage
    surf_spec_cov = []
    for i in range(n_spec):
        surf_spec_cov.append(spec_vecs[i] * surf_dent_vec[i]/ total_sites)
    
    return surf_spec_cov
    
'''
#%%     
def plot_single_traj(t_vec, spec_vec, *arg):

    fig = plt.figure()
    ax = plt.axes()
    colors = ['green', 'red', 'purple', 'grey']
    atm = 'Pd5'
    i = 0
    plt.plot(t_vec,spec_vec, color = colors[i])
    for arg in argv:
        i = i+1
        plt.plot(t_vec,spec_vec, color = colors[i+1])

    
    plt.xlim(0, 0.3)
    plt.legend(['Initial Coverage = 1%', 'Initial Coverage = 5%',  
                'Initial Coverage = 10%',  'Initial Coverage = 25%'])
    plt.ylabel('Percentage Coverage (%)',  fontsize=15)
    plt.xlabel('Time (s)',  fontsize=15)
    lat_dim = 10
    title = atm + '_lattice' + str(lat_dim) + '*' +  str(lat_dim)
    plt.title(title,  fontsize=15)
    plt.savefig(atm+ '.png') 


'''
#%%
'''
x = read_Sim()
y = read_Single_Spec()
z = read_Multiple_Spec(10)
t = z.t

#%%
pd2 = z.multi_spec_ave_df['Pd2*']
#pd2 = num_to_cov(25, 2, pd2)
plt.plot(t,pd2)
'''