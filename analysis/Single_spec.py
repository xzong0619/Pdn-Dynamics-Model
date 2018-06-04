# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 13:43:16 2018

@author: wangyf
"""
import os
import numpy as np
import matplotlib.pyplot as plt

def plot_bar(label, data1, data2):
    fig, ax = plt.subplots()
    index = np.arange(len(label))    
    bar_width = 0.35
    opacity = 0.8    
    plt.bar(index, data1, bar_width,
                     alpha=opacity,
                     color='b',
                     label='kMC')
     
    plt.bar(index + bar_width, data2, bar_width,
                     alpha=opacity,
                     color='g',
                     label='Boltzmann')
     
    plt.xlabel('Pdn States', fontsize=20)
    plt.ylabel('Frequency', fontsize=20)
    plt.xticks(index + bar_width, label, fontsize=20)
    plt.ylim((0, 0.5))
    plt.legend()    
    plt.tight_layout()
    #plt.show()  
    plt.savefig('Distribution.png')
    
class spec_info:
    
    def __init__(self,  data,  dict):
        
        '''
        Fill reactions with species name and species data
        '''
        self.entry = int(data[dict['entry']])        
        self.nevents = int(data[dict['nevents']])     
        self.t = float(data[dict['time']])
        self.T = float(data[dict['temperature']])
        self.E = float(data[dict['energy']])
        self.Pd1 = int(data[dict['pd1*']])
        self.Pd2 = int(data[dict['pd2*']])
        self.Pd3 = int(data[dict['pd3*']])
        self.Pd4 = int(data[dict['pd4*']])
        self.Pd5 = int(data[dict['pd5*']])
        self.CO = int(data[dict['co']])

class single_spec:
    
    def __init__(self, filepath):
        
        self.filepath = filepath

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
        self.Pd1 = np.zeros(n)
        self.Pd2 = np.zeros(n)
        self.Pd3 = np.zeros(n)
        self.Pd4 = np.zeros(n)
        self.Pd5 = np.zeros(n)
        
        for i in range(n):
            self.t[i] = spec[i].t
            self.Pd1[i] = spec[i].Pd1
            self.Pd2[i] = spec[i].Pd2
            self.Pd3[i] = spec[i].Pd3
            self.Pd4[i] = spec[i].Pd4
            self.Pd5[i] = spec[i].Pd5
        
    def num_to_cov(self, lattice_dim):
        
        sites_per_unitcell = 4
        total_sites = lattice_dim**2 *sites_per_unitcell
        # in percentage coverage
        self.Pd1 = self.Pd1/total_sites *100
        self.Pd2 = self.Pd2/total_sites *100 *2
        self.Pd3 = self.Pd3/total_sites *100 *3
        self.Pd4 = self.Pd4/total_sites *100 *4
        self.Pd5 = self.Pd5/total_sites *100 *5

    


#Plot Time Average Count of Species 
class process_single():
    
    def __init__(self, fname):
        
        lat_dim = 10
        Base_path = os.getcwd()
        filepath = os.path.join(Base_path, fname)
        self.c = single_spec(filepath)
        self.c.num_to_cov(lat_dim)
        
        
fname1 = 'specnum_output_1.txt'        
cov1 = process_single(fname1).c   

fname2 = 'specnum_output_2.txt'        
cov2 = process_single(fname2).c 

fname3 = 'specnum_output_3.txt'        
cov3 = process_single(fname3).c 

fname4 = 'specnum_output_4.txt'        
cov4 = process_single(fname4).c 
#%%     
fig = plt.figure()
ax = plt.axes()

atm = 'Pd5'
plt.plot(cov1.t,cov1.Pd5, 'g')
plt.plot(cov2.t,cov2.Pd5, 'r')
plt.plot(cov3.t,cov3.Pd5,  color = 'purple')
plt.plot(cov4.t,cov4.Pd5, color = 'grey')

plt.xlim(0, 0.3)
plt.legend(['Initial Coverage = 1%', 'Initial Coverage = 5%',  
            'Initial Coverage = 10%',  'Initial Coverage = 25%'])
plt.ylabel('Percentage Coverage (%)',  fontsize=15)
plt.xlabel('Time (s)',  fontsize=15)
lat_dim = 10
title = atm + '_lattice' + str(lat_dim) + '*' +  str(lat_dim)
plt.title(title,  fontsize=15)
plt.savefig(atm+ '.png') 
 
    
    
    
    
    
    