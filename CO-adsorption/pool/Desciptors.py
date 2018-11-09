# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 16:08:27 2018

@author: wangyf
"""

import glob
from ase.io import read
import pandas as pd
import numpy as np
from ase.data import covalent_radii
from Pdbulk import NN1,NN2

'''
reading CONTCAR file and save atoms object in strutures 
'''
structures = []
files = []
for file in glob.glob("*CONTCAR"):
    files.append(file)
    structures.append(read(file))

#%% descriptor functions and class

Pdr = covalent_radii[46]
Or = covalent_radii[8]
Cr = covalent_radii[6]
CO = Cr + Or
PdO = Pdr + Or
PdC = Cr + Pdr


def sort_i_and_d(D,I):
    
    Dsort = np.sort(D)
    Dsort = list(Dsort)
    D = list(D)
    Isort = []
    for d in Dsort:
        Isort.append(I[D.index(d)])
        
    return Dsort,Isort

    
# write a class for it
class PdCO():

    def __init__(self, nPd):
        
        '''
        Initializing descriptor variables
        '''
    
        
        self.NPd = nPd #No1
        self.Nsites = [] #No2         
        self.NN1_wavg = [] #No3
        self.NN2_wavg = [] #No4
        self.Dsupport = [] #No5
        
        self.sitetype  = []
        self.descriptors = []
        
    def get_descriptors(self, atoms, CI_PdC = 0.3):
        '''
        Takes in atoms onject 
        Input tolerance for Pd-C bond 0.3 A for detecting CO ads site
        '''
        
        
        '''
        Count number of atoms
        '''
        # Atom index
        Pdi = []
        Cei = []
        Ci = []
        
        for i, atom in enumerate(atoms):
            if atom.symbol == 'Pd': Pdi.append(i)
            if atom.symbol == 'Ce': Cei.append(i)
            if atom.symbol == 'C':  Ci.append(i)
        
        #No of Pd atoms in the cluster
        self.NPd = len(Pdi)   
        
        '''
        Get the bond length of C-Ce, Pd-C, Pd-Pd
        and NN table for each Pd
        '''

        Pd_C = atoms.get_distances(Ci[0], Pdi) #all Pd-C bond length        
        Pd_C, Pdi = sort_i_and_d(Pd_C, Pdi) #sorted Pd-C bond length
        
        C_Ce = atoms.get_distances(Ci[0], Cei)#all Ce-C bond length
        C_Ce, Cei = sort_i_and_d(C_Ce, Cei) #sorted Ce-C bond length
        
        Pd_Pd = pd.DataFrame() #Pd to Pd bond length table 
        PdNN = pd.DataFrame() #Pd NN table         
        for i in Pdi:
           PdD =  atoms.get_distances(i, Pdi)
           PdD, Pdisort = sort_i_and_d(PdD,Pdi)
           Pd_Pd['Pd'+str(i)] =  PdD
           Pd_Pd['i'+str(i)] = Pdisort
           PdNN['Pd'+str(i)] = [sum(np.logical_and(PdD>=NN1[0], PdD<=NN1[1])),
                sum(np.logical_and(PdD>=NN2[0], PdD<=NN2[1]))]
        PdNN.index = ['NN1','NN2']   
        
        #take the shortest distance of CO to Ce as the distance to support
        self.Dsupport = C_Ce[0] 
        
        '''
        Detect CO adsorption sites
        '''
          
        PdC_range = (PdC - CI_PdC, PdC + CI_PdC)   
              
        #Atom index of CO adsorption sites
        COsites = np.array(Pdi)[np.logical_and(Pd_C>=PdC_range[0], Pd_C<=PdC_range[1])]
        self.Nsites = len(COsites)
        
        if self.Nsites== 3: self.sitetype = 'h'
        if self.Nsites == 2: self.sitetype = 'b'
        if self.Nsites == 1: self.sitetype = 't'
        
        COsites_cols = []
        for s in range(len(COsites)):
            COsites_cols.append('Pd'+str(COsites[s]))
        
        PdNN_CO = PdNN.loc[:, COsites_cols] #NN dataframe 
        Pd_C_CO = Pd_C[:len(COsites)] #CO distance to neighboring Pd atoms
        
        '''
        Weighted average for NN1, NN2
        '''
        
        weights = np.array(Pd_C_CO)/np.sum(Pd_C_CO) #weights based on CO-Pd distance
        self.NN1_wavg = np.dot(weights, PdNN_CO.loc['NN1'].values)
        self.NN2_wavg = np.dot(weights, PdNN_CO.loc['NN2'].values)
        
        '''
        Make a column in data frame
        '''
        
        self.descriptors =  [self.NPd, self.NN1_wavg, self.NN2_wavg, self.Dsupport]
        # take out the number of sites and NPd

#%% Analyse the structures
        
for struct in structures:


    PdCO_ob = PdCO(struct
    PdCO_ob.get_descriptors())
    
dem = np.zeros((Ntot, len(labels)))


for i in range(Ntot):
    NPd = i + 1
    col = str(NPd)
    PdCO_ob = PdCO(NPd)
    PdCO_ob.get_descriptors()
    dv = PdCO_ob.descriptors
    dem[i,:-1] = np.array(dv)
dem[:,-1] = np.array(charge)






        
        