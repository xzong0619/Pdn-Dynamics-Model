# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 11:42:20 2018

@author: wangyf
"""

import os 
import sys
import numpy as np
import random as _random
from shutil import copyfile
import Cluster_structure as cs
import networkx as nx

'''
Pdn Constructing wrapper 

'''

'''
 Class definitions for:
     User input file        Class name
     (in user_inputs)
     -----------------      ------------
     rxn_input.txt          RxnIn
     cluster_input.txt      ClusterIn
     
     Zacros input file      Class name
     -----------------      ------------
     simulation_input.dat   SimOut
     energetics_input.dat   ClusterOut
     mechanism_input.dat    MechanismOut
     state_input.dat        StateOut

'''


'''
============ Classes to read user input file ============
'''
def ReadIn(filename):
    '''
    Read reaction network file input by the user
    '''
    HomePath = os.path.expanduser('~')
    IO_path = os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Generator')
    input_dir = 'user_inputs'
    filepath = os.path.join(IO_path, input_dir, filename)
    fid = open(filepath, 'r')
    file = fid.read()
    lines = file.splitlines()
    dict_array = lines[0].lower().split('\t')
    data = lines[1:]
    dict = {}
    for x in range(0, len(dict_array)):
        dict[dict_array[x]] = x
    fid.close()
    return (data,dict)


def ReadRxn(filename):
    
    
    data =  ReadIn(filename)[0]
    dict =  ReadIn(filename)[1]
    
    rxn = []
    for s in data:
        rxn.append(RxnIn(s.split('\t'), dict))
    return rxn

def ReadCluster(filename):
    
    data =  ReadIn(filename)[0]
    dict =  ReadIn(filename)[1]
    cluster = []
    
    for s in data:
        cluster.append(ClusterIn(s.split('\t'),dict))
    return cluster

class x_spec():

    def __init__(self,  x,  data, dict):
        
        self.n = int(data[dict[x + '_n']]) # total atoms in A
        self.no_layer = int(data[dict[x + '_no_layers']])
        self.n_1 = int(data[dict[x + '_1st']]) # surface species dentate number
        self.n_2 = int(data[dict[x + '_2nd']])
        self.neighboring = int(data[dict[x + 'i']]) #neighboring list index
        
    
class RxnIn():
        
    '''
    Fill reactions with species name and associated kinetic data
    '''
    def __init__(self,  data,  dict):

        self.no = int(data[dict['rxn_no']])          #reaction no
        self.type = str(data[dict['rxn_type']])     
        '''
        reaction type: O(Ostwald rippeing),C (Coalescence),
                        H(2D to 3D hopping), D (diffusion)
        '''                
        
        self.Ea = float(data[dict['ea']])
        self.Af = float(data[dict['afwd']])
        self.Aratio = float(data[dict['aratio']])
        
        self.ai = int(data[dict['ai']]) 

        if not data[dict['bi']] == '':      
            self.bi = int(data[dict['bi']]) 
        else:
            self.bi = []
        
        if not data[dict['ci']] == '':
            self.ci = int(data[dict['ci']]) 
        else:
            self.ci = []
            
            
class ClusterIn():            
    
    
    def __init__(self,  data,  dict):

        self.surf_spec = str(data[dict['species_name']] + '*') #surface species name
        self.ei = float(data[dict['ei']])       
        self.a = x_spec('a',  data, dict)
        self.surf_dent = self.a.n_1
        self.gm = int(data[dict['graph_multiplicity']]) 

class Input():
    
    def __init__(self): 
        
        self.f_rxn = 'rxn_input.txt'
        self.f_cluster = 'cluster_input.txt'          
        self.rxn  = ReadRxn(self.f_rxn)
        self.cluster = ReadCluster(self.f_cluster)
        
        
        self.n_surf = len(self.cluster) # number of surface species
        
        self.surf_n = [] # number of atoms in one cluster
        for i in range(self.n_surf):
            self.surf_n.append(self.cluster[i].a.n)
        self.n_Pdmax = max(self.surf_n)
        
        self.surf_spec = [] # name of surface species
        for i in range(self.n_surf):
            self.surf_spec.append(self.cluster[i].surf_spec)  
            
        self.surf_dent = []
        for i in range(self.n_surf):
            self.surf_dent.append(self.cluster[i].surf_dent)  
        
        self.gm = []
        for i in range(self.n_surf):
            self.gm.append(self.cluster[i].gm) 
             
        self.fldr = os.getcwd() 
            


class NCIn():
    
    def __init__(self, str1):
        
        self.l= [str1] 
        self.n = len(self.l)
        
    def addnc(self,str2):
        self.l.append(str2)
        self.n = len(self.l)           
                


class RxnCal(Input):
    '''
    Use in debugging
    Calculate the propensity, reaction energy and activation energy 
    for forward and reverse reaction steps 
    '''
    kb = 8.61733E-05 #eV/K
    
    def __init__(self):
        
        super().__init__()
        
        '''
        Process reaction information, calculate activation energy and pre-factors
        
        '''

        self.n_rxn = len(self.rxn) # number of surface reactions
        
        self.type = []
        for i in range(self.n_rxn):
            self.type.append(self.rxn[i].type)
            
        self.Ea = []
        for i in range(self.n_rxn):
            self.Ea.append(self.rxn[i].Ea)
            
        self.Af = []
        for i in range(self.n_rxn):
            self.Af.append(self.rxn[i].Af)
            
        self.Aratio = []
        for i in range(self.n_rxn):
            self.Aratio.append(self.rxn[i].Aratio)
        
        self.ai= [] # index of a in the neirghboring list
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
              
        self.Ei = []
        for i in range(self.n_surf):
            self.Ei.append(self.cluster[i].ei) 
        
        
        
        E_initial = []
        E_final = []
        
        self.Prev = []
    
       
        for r in range(self.n_rxn):
            if self.type[r] == 'O' or self.type[r] == 'C':
                
                ai = self.ai[r] - 1
                bi = self.bi[r] - 1 
                ci = self.ci[r] - 1
            
                E_initial.append(self.Ei[ai]  +  self.Ei[bi])
                E_final.append(self.Ei[ci])
                
            if self.type[r] == 'H':
                
                ai = self.ai[r] - 1
                bi = self.bi[r] - 1 
                E_initial.append(self.Ei[ai])
                E_final.append(self.Ei[bi])
                 
            if self.type[r] == 'D':
                
                ai = self.ai[r] - 1
                E_initial.append(self.Ei[ai])
                E_final.append(self.Ei[ai])
        
        self.E_rxn = list(np.array(E_final) - np.array(E_initial)) #forward reaction energy
        self.Ea_rev = list(np.array(self.Ea) - np.array(self.E_rxn)) # reverse reaction activation energy
        
        self.Arev = list(np.divide(np.array(self.Af), np.array(self.Aratio)))
        
        
        
    def Tdependence(self, T):
        
        self.T  = T # temperature
        #forward propensity
        self.Pfwd = list(np.multiply(np.array(self.Af), np.exp(- np.array(self.Ea)/self.kb/self.T)))
        #reverse propnesity
        self.Prev = list(np.multiply(np.array(self.Arev), np.exp(- np.array(self.Ea_rev)/self.kb/self.T)))
                 
        self.Timescale_fwd = list(1/np.array(self.Pfwd))
        self.Timescale_rev = list(1/np.array(self.Prev))





                       
class SimOut(Input):            
    '''
    Write simulation_input.dat
    '''
    fname = 'simulation_input.dat'
    '''
    f_rxn = 'rxn_input.txt'
    f_cluster = 'cluster_input.txt'
    
    fldr = os.getcwd()
    
    rxn  = ReadRxn(f_rxn)
    cluster = ReadCluster(f_cluster)
    '''
    
    def __init__(self,   n_CO,   T,  Pt,  n_hours,  flag):
        
        super().__init__()
        
        self.Seed = None  
        
        self.T  = T # temperature
        self.Pt = Pt # total pressure
        
        # Put flag as an optional input
        self.flg_on_events = flag[0]
        self.flg_debug = flag[1]
        self.flg_no_restart = flag[2]
        

        self.n_gas = 1 # gas species involved (only CO)
        self.E_gas = 0 # energy of the gas species
        self.MW_gas = 28.0102 # MW of CO
        self.Y_gas = 1 # gas fraction of CO
        
        self.max_time = 1e5
        self.time_snap = self.max_time/10
        self.time_stats = self.max_time/10
        self.wall_time = int(n_hours*3600)
        
        self.max_steps = int(1e2)
        self.event_snap = int(self.max_steps/10)
        self.event_stats = int(self.max_steps/100)
                    
        self.n_CO = n_CO # if there is CO
        '''
        self.n_surf = len(self.cluster) # number of surface species
        
        self.surf_n = []
        for i in range(self.n_surf):
            self.surf_n.append(self.cluster[i].a.n)
        self.n_Pdmax = max(self.surf_n)
        
        self.surf_spec = []
        for i in range(self.n_surf):
            self.surf_spec.append(self.cluster[i].surf_spec)  
            
        self.surf_dent = []
        for i in range(self.n_surf):
            self.surf_dent.append(self.cluster[i].surf_dent) 

        '''
    

    def WriteIn(self, fldr = None):
        
        self.output_dir = 'default_outputs'
        
        if not fldr == None:
            
            self.fldr = fldr 
        else:
            self.fldr = os.path.join(os.getcwd, self.output_dir)
            
            
        with open(os.path.join(self.fldr, self.fname), 'w') as txt:
            
            SeedTxt = ''
            if self.Seed is None:
                _random.seed()
                self.Seed = _random.randint(10000, 99999)
                SeedTxt = ' # Random seed from Python wrapper'
                
            # Write file title
            txt.write('# Pd growth model including upto {:2d} Pd atoms with {:2d} CO\n'
                      .format(self.n_Pdmax, self.n_CO))
            txt.write('############################################################################\n')
                      
            # Write the main part
            txt.write('{:20}{:5d}{}\n\n'.format('random_seed',
                      self.Seed, SeedTxt))
            txt.write('{:20}{:5.1f}\n'.format('temperature', self.T))
            txt.write('{:20}{:5.5e}\n\n'.format('pressure', self.Pt*self.n_CO))
            
            txt.write('{:20}{}\n'.format('n_gas_species',
                      str(self.n_gas)))
            txt.write('{:20}{}\n'.format('gas_specs_names','CO'))
            txt.write('{:20}{:5.3f}\n'.format('gas_energies', self.E_gas))
            txt.write('{:20}{:5.3f}\n'.format('gas_molec_weights', self.MW_gas))
            txt.write('{:20}{:5.3f}\n\n'.format('gas_molar_fracs', self.Y_gas))
            
            txt.write('{:20}{}\n'.format('n_surf_species', self.n_surf))
            txt.write('{:20} '.format('surf_specs_names'))
            for ns in self.surf_spec:
                txt.write('{} '.format(ns))
            txt.write('\n{:20}'.format('surf_specs_dent'))
            
            for sd in self.surf_dent:
                txt.write('{:2d}'.format(sd))
            txt.write('\n\n')
            
            if self.flg_on_events:
                txt.write('{:20}{:10}{:5d}\n'
                          .format('snapshots','on event',self.event_snap))
                txt.write('{:20}{:10}{:5d}\n'
                          .format('process_statistics','on event',self.event_stats))
                txt.write('{:20}{:10}{:5d}\n\n'
                          .format('species_numbers','on event',self.event_stats))
            else:
                txt.write('{:20}  {:10} {:5.3f}\n'
                          .format('snapshots','on time',self.time_snap))
                txt.write('{:20}  {:10} {:5.3f}\n'
                          .format('process_statistics','on time',self.time_stats))
                txt.write('{:20}  {:10} {:5.3f}\n\n'
                          .format('species_numbers','on time',self.time_stats))
                
            if self.flg_debug:
                txt.write('{:20}{:5}\n\n'.format('event_report','on'))
                
            if self.flg_on_events:
                txt.write('{:20}{:5d}\n\n'
                          .format('max_steps',self.max_steps))
            else:
                txt.write('{:20}{:5.3f}\n\n'
                          .format('max_time',self.max_time))
            
            txt.write('{:20} {:5d}'.format('wall_time', self.wall_time))
            
            if self.flg_no_restart:
                txt.write('\n{:20}\n'.format('no_restart'))
                
            if self.flg_debug:
                txt.write('debug_report_global_energetics\n')
                txt.write('debug_report_processes\n')
                txt.write('debug_check_processes\n')
                txt.write('# debug_check_lattice\n\n')
                          
            txt.write('\n{:20}\n'.format('finish'))
            
            print('25%   simulation_input.dat  generated\n')


            
class MechanismOut(Input):
    
    fname = 'mechanism_input.dat'
    
    def __init__(self, neighboring_list = None):
        
        super().__init__()
        '''
        Process cluster information
        '''
        
        if neighboring_list == None:
            self.neighboring_list = []
            for i in range(self.n_surf):
                self.neighboring_list.append('1-2')
        else: 
            self.neighboring_list = neighboring_list
        
        # Initialize nc object
        self.nc = []
        for c in range(self.n_surf):
            self.nc.append(NCIn(neighboring_list[c]))

        
        '''
        Process rxn information
        '''

        self.n_rxn = len(self.rxn) # number of surface species
        
        self.type = []
        for i in range(self.n_rxn):
            self.type.append(self.rxn[i].type)
            
        self.Ea = []
        for i in range(self.n_rxn):
            self.Ea.append(self.rxn[i].Ea)
            
        self.Af = []
        for i in range(self.n_rxn):
            self.Af.append(self.rxn[i].Af)
            
        self.Aratio = []
        for i in range(self.n_rxn):
            self.Aratio.append(self.rxn[i].Aratio)
        
        self.ai= [] # index of a in the neirghboring list
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
                
        self.diffusion_i = [] # diffusion reaction index
        self.diffusion_ai= []
        for i in range(self.n_rxn):
            if self.type[i] == 'D':
                self.diffusion_i.append(self.rxn[i].no)
                self.diffusion_ai.append(self.ai[i])
        
        
        self.Cal = RxnCal()
        self.Ea_rev = self.Cal.Ea_rev
        self.Arev = self.Cal.Arev
        	
    def NC_Generator(self, r):
        
        '''
        r is the reaction index in Python
        '''
        
        # Initialize iteration variables in the loop
        
        ni_u = 1
        more_structure_flag  = 0
        new_H_list = None
        if not r == 27:
        
            if self.type[r] == 'O' or self.type[r] == 'C':
                # A(small clusters) + B -> C
                # Convert cluster index into python index by -1
                ai_py = self.ai[r] - 1
                bi_py = self.bi[r] - 1 
                ci_py = self.ci[r] - 1
                
                if self.surf_dent[ci_py] > 5:
                    
                    bnc = cs.neighbor_string_list([self.neighboring_list[bi_py]])
                    cnc = cs.neighbor_string_list([self.neighboring_list[ci_py]])
                    
                    new_H = [bnc]
                    new_wt_i = [cs.neighbor_weight(bnc)]
                    n =  self.surf_dent[bi_py]
                    for i in range(self.surf_dent[ai_py]):
                        
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
                                new_H_list = cs.neighbor_list_string(new_H[i])
                                self.nc[ci_py].addnc(new_H_list) 
                            
        if self.type[r] == 'H':
            # B + * -> A
            # Convert cluster index into python index by -1
            ai_py = self.ai[r] - 1
            bi_py = self.bi[r] - 1 
            if self.surf_dent[ai_py] > 5:
                
                anc = cs.neighbor_string_list([self.neighboring_list[ai_py]])
                bnc = cs.neighbor_string_list([self.neighboring_list[bi_py]])
                new_H = bnc
                new_wt_i = cs.neighbor_weight(bnc)
                n =  self.surf_dent[bi_py]
                    
                n = n+1 
                Pdx = cs.Pdn(n, new_H, new_wt_i)
                Pdx.new_graph()
                ni_u = Pdx.ni_u
                new_H = Pdx.H
                new_wt_i = Pdx.wt_i
                
                for i in range(ni_u):
                    A = nx.Graph(new_H[i])
                    B = nx.Graph(anc)
                    if nx.is_isomorphic(A,B):
                        if not set(new_H[i]) == set(anc):
                            more_structure_flag  = 1
                            new_H_list = cs.neighbor_list_string(new_H[i])
                            self.nc[ai_py].addnc(new_H_list) 
        
        return more_structure_flag, new_H_list        
                       
    def NC_cleanup(nc):
        
        
        '''
        takes in nc object and clean it up
        '''
        
        for c in range(len(nc)):
            if nc[c].n > 1:
                H_list = []
                for cn in range(nc[c].n):
                    H_list.append(cs.neighbor_string_list([nc[c].l[cn]]))
                i_list = cs.clean_up(H_list)
                ni = len(i_list)
                nc[c].l = []
                for i in range(ni):
                    H_str = cs.neighbor_list_string(H_list[i_list[i]])
                    nc[c].addnc(H_str)
        return nc
        
        
        
    def Write_Single_rxn(self, r, txt):    
            
        rxn_tp = {'O':'Ostwald Rippening',
                  'C': 'Coalescence',
                  'H': 'Hopping',
                  'D': 'Diffusion'}
        
            
        if self.type[r] == 'O' or self.type[r] == 'C':
                
            ai = self.ai[r] - 1
            bi = self.bi[r] - 1 
            ci = self.ci[r] - 1
            
            neighboring_list = self.neighboring_list[ci]
            [more_structure_flag, new_H_list] = self.NC_Generator(r)  
            if more_structure_flag: 
                neighboring_list = new_H_list
            
            txt.write('############################################################################\n')
            txt.write('#{:3d} {} {} + {} -> {} \n'
                      .format(r+1, rxn_tp[self.type[r]], self.surf_spec[ai], 
                              self.surf_spec[bi], self.surf_spec[ci]))
            txt.write('{} {:2d}\n'.format('reversible_step', r+1))
            txt.write('   {} {}\n'.format('neighboring',neighboring_list))
            txt.write('   {} {:2d}\n'.format('sites', self.surf_dent[ci]))
            txt.write('   initial\n')
            for bdent_i in range(self.surf_dent[bi]):
                txt.write('      {} {} {}\n'
                          .format('1', self.surf_spec[bi], bdent_i+1))
            for adent_i in range(self.surf_dent[ai]):
                txt.write('      {} {} {}\n'
                          .format('2', self.surf_spec[ai], adent_i+1))
            txt.write('   final\n')
            for cdent_i in range(self.surf_dent[ci]):
                txt.write('      {} {} {}\n'
                          .format('1', self.surf_spec[ci], cdent_i+1))
            txt.write('   {} {}\n'.format('variant','_'))
            txt.write('   site_types')
            for i_dent in range(self.surf_dent[ci]):
                txt.write(' {}'.format('T'))

            
        '''
        2D to 3D hopping reactions
        #B + * - > A
        '''
        
        if self.type[r] == 'H':
            
            
            
            
            ai = self.ai[r] - 1 
            bi = self.bi[r] - 1
            
            neighboring_list = self.neighboring_list[ai]
            [more_structure_flag, new_H_list] = self.NC_Generator(r)  
            if more_structure_flag: 
                neighboring_list = new_H_list
            
            txt.write('############################################################################\n')
            txt.write('#{:3d} {} {} + * -> {}\n'
                      .format(r+1, rxn_tp[self.type[r]],
                              self.surf_spec[bi], 
                              self.surf_spec[ai]))
            txt.write('{} {:2d}\n'.format('reversible_step', r+1))
            txt.write('   {} {}\n'.format('neighboring',neighboring_list))
            txt.write('   {} {:2d}\n'.format('sites', self.surf_dent[ai]))
            txt.write('   initial\n')
            for bdent_i in range(self.surf_dent[bi]):
                txt.write('      {} {} {}\n'
                          .format('1', self.surf_spec[bi], bdent_i+1))
            txt.write('      {} {} {}\n'.format('2', '*', 1))
            txt.write('   final\n')
            for adent_i in range(self.surf_dent[ai]):
                txt.write('      {} {} {}\n'
                          .format('1', self.surf_spec[ai], adent_i+1))
            
            txt.write('   {} {}\n'.format('variant','_'))
            txt.write('   site_types')
            for i_dent in range(self.surf_dent[ai]):
                txt.write(' {}'.format('T'))

        
        '''
        Small clusters Diffusion Pd1-Pd4
        '''
        
        if self.type[r] == 'D':
            
            ai = self.ai[r] - 1
            a_ni = self.surf_dent.index(self.surf_dent[ai]+1) 
            
            
            txt.write('############################################################################\n')
            txt.write('#{:3d} {} {}\n'
                      .format(r+1, rxn_tp[self.type[r]], 
                              self.surf_spec[ai]))
            txt.write('{} {:2d}\n'.format('reversible_step', r+1))
            txt.write('   {} {}\n'.format('neighboring',self.neighboring_list[a_ni]))
            txt.write('   {} {:2d}\n'.format('sites', self.surf_dent[ai] + 1))
            txt.write('   initial\n')
            for adent_i in range(self.surf_dent[ai]):
                txt.write('      {} {} {}\n'
                          .format('1', self.surf_spec[ai], adent_i+1))
            
            txt.write('      {} {} {}\n'.format('2', '*', 1))
            txt.write('   final\n')
            txt.write('      {} {} {}\n'.format('1', '*', 1))
            for adent_i in range(self.surf_dent[ai]):
                txt.write('      {} {} {}\n'
                          .format('2', self.surf_spec[ai], adent_i+1))
            txt.write('   {} {}\n'.format('variant','_'))
            txt.write('   site_types')
            for i_dent in range(self.surf_dent[ai] + 1):
                txt.write(' {}'.format('T'))

        
        else:
            pass
        if self.type[r] == 'H':
            txt.write('\n   {} {:.5e}\n'.format('pre_expon', self.Arev[r]))
            txt.write('   {} {:.5e}\n'.format('pe_ratio', 1/self.Aratio[r]))
            txt.write('   {} {:5.3f}\n'.format('activ_eng', self.Ea_rev[r]))
        else:
            txt.write('\n   {} {:.5e}\n'.format('pre_expon', self.Af[r]))
            txt.write('   {} {:.5e}\n'.format('pe_ratio', self.Aratio[r]))
            txt.write('   {} {:5.3f}\n'.format('activ_eng', self.Ea[r]))
        txt.write('   {} {:5.3f}\n'.format('prox_factor', 0))
        txt.write('   end_variant\n')
        txt.write('end_reversible_step\n\n') 
        
    def WriteIn(self,  fldr = None, r_index = None):
        
        
        
        if not fldr == None:
            
            self.fldr = fldr 
        
        with open(os.path.join(self.fldr, self.fname), 'w') as txt:
            
            # Write file title
            txt.write('# Pd growth model including upto {:2d} Pd atoms\n'
                      .format(self.n_Pdmax))
            txt.write('#{:3d} {}\n\n'.format(self.n_rxn, 'reactions'))
            txt.write('mechanism\n\n')  
            
            if not r_index == None:
                
                r = r_index
                self.Write_Single_rxn(r, txt)
                
                
                if self.type[r] == 'O' or self.type[r] == 'C':
                    
                    ai = self.diffusion_ai.index(self.ai[r])
                    print(self.diffusion_i[ai])
                    self.Write_Single_rxn(self.diffusion_i[ai] -1 ,  txt)
                    
                    
                    if self.bi[r] in self.diffusion_ai:
                        bi = self.diffusion_ai.index(self.bi[r])
                        if not ai == bi:
                            print(self.diffusion_i[bi])
                            self.Write_Single_rxn(self.diffusion_i[bi] -1 ,  txt)
                                      
                    
            else:
                # write the main part
                for r in range(self.n_rxn):
                    self.Write_Single_rxn(r, txt)
                        
                    
            txt.write('############################################################################\n') 
            txt.write('end_mechanism\n')
            print('75%   mechanism_input.dat  generated\n')
            
class ClusterOut(Input):
    
    fname = 'energetics_input.dat'

    def __init__(self, neighboring_list = None):
        
        super().__init__()
        
        self.n_rxn = len(self.rxn) # number of surface reactions
        
        if neighboring_list == None:
            self.neighboring_list = []
            for i in range(self.n_surf):
                self.neighboring_list.append('1-2')
        else: 
            self.neighboring_list = neighboring_list
        
        Mechan_info = MechanismOut(self.neighboring_list)
        for r in range(self.n_rxn):
            Mechan_info.NC_Generator(r)            
        self.nc = MechanismOut.NC_cleanup(Mechan_info.nc)
        
        
        self.nc_total = 0
        for c in range(self.n_surf):
            self.nc_total = self.nc_total + self.nc[c].n
            
 
        self.Ei = []
        for i in range(self.n_surf):
            self.Ei.append(self.cluster[i].ei) 
   
        
    def WriteIn(self, fldr = None):
        
        self.output_dir = 'default_outputs'
        
        if not fldr == None:
            
            self.fldr = fldr 
        else:
            self.fldr = os.path.join(os.getcwd, self.output_dir)
            
        with open(os.path.join(self.fldr, self.fname), 'w') as txt:
             
             # Write file title
            txt.write('# Pd growth model including upto {:2d} Pd atoms\n'
                      .format(self.n_Pdmax))
            txt.write('#{:3d} {}\n\n'.format(self.n_surf, 'clusters'))
            txt.write('#{:3d} {}\n\n'.format(self.nc_total, 'patterns'))
            txt.write('energetics\n\n')  
            
            # Write the main part
            for c in  range(self.n_surf):
                
                for cn in range(self.nc[c].n):
                    
                    neighboring_list = self.nc[c].l[cn]
                    
                    txt.write('############################################################################\n')
                    txt.write('#{:3d} - {}{} {}\n'.
                              format(c+1, 'Pd', self.cluster[c].a.n, 'cluster'))
                    txt.write('{} {:10}\n'.format('cluster', self.surf_spec[c]))
                    txt.write('   {} {:2d}\n'.format('sites',self.surf_dent[c]))                          
                    if self.surf_dent[c] > 1:
                        txt.write('   {} {} {} {:2d}\n'
                                  .format('neighboring', neighboring_list, 
                                          '#Pattern', cn+1))
                    txt.write('   lattice_state\n')
                    for i_dent in range(self.surf_dent[c]):
                        txt.write('     {:2d} {:10}{:2d}\n'.format(1, self.surf_spec[c], i_dent+1))
                    txt.write('   {} {}\n'.format('variant','_'))
                    txt.write('   site_types')
                    for i_dent in range(self.surf_dent[c]):
                        txt.write(' {}'.format('T'))
                    txt.write('\n   {} {:2d}\n'.format('graph_multiplicity', self.gm[c]))
                    txt.write('   {} {:5.3f}\n'.format('cluster_eng', self.Ei[c]*self.gm[c]))
                    txt.write('   end_variant\n')
                    txt.write('end_cluster\n\n')
            txt.write('############################################################################\n') 
            txt.write('end_energetics')
            print('50%   energetics_input.dat  generated\n')            
            


       
class StateOut(Input):
    
    fname = 'state_input.dat'
    rxn = ReadRxn('rxn_input.txt')
    
    def __init__(self,  neighboring_list = None):
        
        super().__init__()
        
        
        
        if neighboring_list == None:
            self.neighboring_list = []
            for i in range(self.n_surf):
                self.neighboring_list.append('1-2')
        else: 
            self.neighboring_list = neighboring_list
        
        
        self.cluster_no = [] # number of the same clusters to be seeded
        
    def Seed_Single_cluster(self,  cluster_i, n, txt):


       
        txt.write('   {} {} {:5d}\n'
                  .format('seed_multiple',self.surf_spec[cluster_i], n))
            
        txt.write('      site_types')
        for i_dent in range(self.surf_dent[cluster_i]):
            txt.write(' {}'.format('T'))
        txt.write('\n')
        if self.surf_dent[cluster_i] == 1:
            pass
        else:
            txt.write('      {} {}\n'
                      .format('neighboring', self.neighboring_list[cluster_i]))
        txt.write('   end_seed_multiple\n\n')
        
     
    def WriteIn(self, fldr, cluster_list,  cluster_no,  n_cluster):
       
        self.fldr = fldr
        self.cluster_list = cluster_list
        self.cluster_no = cluster_no  # number of same cluster
        self.n_cluster = n_cluster # n types of clusters
        

            
        with open(os.path.join(self.fldr, self.fname), 'w') as txt:
            # Write file title
            txt.write('# Pd growth model including upto {:2d} Pd atoms\n'
                      .format(self.n_Pdmax))
            txt.write('initial_state\n\n')  
            
            # write the main part
            # If the input is a list
            if self.n_cluster > 1 : 
            
                for s in range(self.n_cluster):
                    
                    cluster_i = self.surf_spec.index(self.cluster_list[s])
                    self.Seed_Single_cluster(cluster_i,  cluster_no[s],  txt)

            # Only one 1 input
            if self.n_cluster == 1:
                
                cluster_i = self.surf_spec.index(self.cluster_list)
                
                self.Seed_Single_cluster(cluster_i,   cluster_no,  txt)
                
            else:
                pass
            
            txt.write('end_initial_state\n')
            
            print('100%   state_input.dat  generated\n')
                    
            
    def WriteIn_Single_rxn(self, fldr,  r):
        
        self.fldr = fldr
        self.ai = self.rxn[r].ai           

        if self.rxn[r].bi == []:
            self.bi = None
        else:
            self.bi = self.rxn[r].bi
        
        if self.rxn[r].ci == []:
            self.ci = None
        else:
            self.ci = self.rxn[r].ci
        
        
        with open(os.path.join(self.fldr, self.fname), 'w') as txt:
        # Write file title
            txt.write('# Pd growth model including upto {:2d} Pd atoms\n'
                      .format(self.n_Pdmax))
            txt.write('initial_state\n\n')  
        
            if self.rxn[r].type == 'O' or self.rxn[r].type == 'C':
                
                sa = 1
                sb = 1
                ai = self.ai - 1
                bi = self.bi - 1 
                if ai == bi:
                    self.Seed_Single_cluster(ai, sa + sb,  txt)
                else:
                    self.Seed_Single_cluster(bi, sb,  txt)
                
                    self.Seed_Single_cluster(ai, sa,  txt)
               
      
            if self.rxn[r].type == 'H':
                
                bi = self.bi - 1 
               
                
                self.Seed_Single_cluster(bi,  1,  txt)

                
            if self.rxn[r].type == 'D':
            
                ai = self.ai - 1
                self.Seed_Single_cluster(ai,  1,  txt)
                
            txt.write('end_initial_state\n')
            
            print('100%   state_input.dat  generated\n')

class MakingCopy:
    
    def __init__(self,  fname,  output_fldr):
        
        
        self.Base_path = os.getcwd()
        self.input_dir = os.path.join(os.getcwd(), 'zacros_inputs')
        self.fname = fname
        self.output_fldr = output_fldr
        self.src = os.path.join(self.input_dir, self.fname)
        self.dst = os.path.join(self.output_fldr, self.fname)
        
        copyfile(self.src, self.dst)
        
    
    
    
    
             
                
                


        
            
            
            
            
            
            
            