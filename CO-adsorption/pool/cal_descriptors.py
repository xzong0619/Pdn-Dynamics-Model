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
from ase.visualize import view
from sympy import Plane, Point3D
import pickle 
'''
read adsoprtion energy and barder charge from a csv file
'''
data = pd.read_csv('adsorption_constant.csv', header = 0)

'''
reading CONTCAR file and save atoms object in strutures 
'''
structures = []

for file in glob.glob("*CONTCAR"):
    structure = {'filename': file,
                 'atoms': read(file)}
    structures.append(structure)
    

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

    def __init__(self):
        
        '''
        Initializing descriptor variables
        '''
    
        self.Eads = [] # CO adsorption energy - y
        self.NPd = [] #No1
        self.Nsites = [] #No2         
        self.NN1_wavg = [] #No3
        self.NN2_wavg = [] #No4
        self.Dsupport = [] #No5
        self.charge = []
        
        self.sitetype  = []
        self.descriptors = []
        
    def get_descriptors(self, structure, data, CI_PdC = 0.2):
        '''
        Takes in a structure dictionary with an atoms object and filename 
        Input tolerance for Pd-C bond 0.3 A for detecting CO ads site
        '''
        self.atoms = structure['atoms']
        self.filename =structure['filename']
        
        self.Eads = float(data[data['Filename'] == self.filename]['Eads'])
        self.charge = float(data[data['Filename'] == self.filename]['Charge'])
        self.realsite = data[data['Filename'] == self.filename]['RealSite'].values[0]
        
        '''
        Count number of atoms
        '''
        # Atom index
        Pdi = []
        Cei = []
        Ci = []
        
        for i, atom in enumerate(self.atoms):
            if atom.symbol == 'Pd': Pdi.append(i)
            if atom.symbol == 'Ce': Cei.append(i)
            if atom.symbol == 'C':  Ci.append(i)
        
        #No of Pd atoms in the cluster
        self.NPd = int(len(Pdi))   
        
        '''
        Get the bond length of C-Ce, Pd-C, Pd-Pd
        and NN table for each Pd
        '''

        Pd_C = self.atoms.get_distances(Ci[0], Pdi, mic = True) #all Pd-C bond length 
        
        Pd_C, Pdi = sort_i_and_d(Pd_C, Pdi) #sorted Pd-C bond length
        
        #The distance of Pd to the first nearest C
        PdC3 = np.zeros(3)
        #The bond tolerance is 
        bond_tol = 0.7
        if len(Pd_C) == 1: 
            
            PdC3[0] = Pd_C[0]
            COsites = np.array(Pdi)[:1] #top
        
        if len(Pd_C) == 2:    
            
            PdC3[:2] = Pd_C[:2]
            diff = Pd_C[1] - Pd_C[0]
            
            if diff < bond_tol: COsites = np.array(Pdi)[:2] #bridge
            else: COsites = np.array(Pdi)[:1] #top
            
        if len(Pd_C) >= 3: 
            
            PdC3 = Pd_C[:3] 
            diff1 = Pd_C[1] - Pd_C[0]
            diff2 = Pd_C[2] - Pd_C[1]
            
            if diff1 > bond_tol: COsites = np.array(Pdi)[:1] #top
            else:
                if diff2 > bond_tol: COsites = np.array(Pdi)[:2] #bridge
                else: COsites = np.array(Pdi)[:3] #hollow
        
        self.PdC1 = PdC3[0]
        self.PdC2 = PdC3[1]
        self.PdC3 = PdC3[2]
        
        
        C_Ce = self.atoms.get_distances(Ci[0], Cei, mic = True)#all Ce-C bond length
        C_Ce, Cei = sort_i_and_d(C_Ce, Cei) #sorted Ce-C bond length
        
        Pd_Pd = pd.DataFrame() #Pd to Pd bond length table 
        PdNN = pd.DataFrame() #Pd NN table         
        for i in Pdi:
           PdD =  self.atoms.get_distances(i, Pdi)
           PdD, Pdisort = sort_i_and_d(PdD,Pdi)
           Pd_Pd['Pd'+str(i)] =  PdD
           Pd_Pd['i'+str(i)] = Pdisort
           PdNN['Pd'+str(i)] = [sum(np.logical_and(PdD>=NN1[0], PdD<=NN1[1])),
                sum(np.logical_and(PdD>=NN2[0], PdD<=NN2[1]))]
        PdNN.index = ['NN1','NN2']   
        
        #take the  distance of CO to Ce plane (determined by 3 Ce points)
        # as the distance to support
        Ce_plane = Plane(Point3D(self.atoms[Cei[0]].position), 
                         Point3D(self.atoms[Cei[1]].position), 
                         Point3D(self.atoms[Cei[2]].position))
        self.Dsupport = float(Ce_plane.distance(Point3D(self.atoms[Ci[0]].position)))
        #self.Dsupport = C_Ce[0] 
        
        '''
        Detect CO adsorption sites
        '''
          
        #PdC_range = (PdC - CI_PdC, PdC + CI_PdC)   
              
        #Atom index of CO adsorption sites
        #COsites = np.array(Pdi)[np.logical_and(Pd_C>=PdC_range[0], Pd_C<=PdC_range[1])]
        self.Nsites = len(COsites)
        
        if self.Nsites== 3: self.sitetype = 'hollow'
        if self.Nsites == 2: self.sitetype = 'bridge'
        if self.Nsites == 1: self.sitetype = 'top'
        
        
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
        Make a row in dataframe as an ID for each structure including filenames and properties etc
        '''
        self.structureID =  [self.filename, #filename
                             self.atoms, # atoms object
                             self.Eads, #Eads
                             self.NPd, #Npd
                             self.realsite, #real sitetype
                             self.sitetype, #sitetype from calculation
                             self.NN1_wavg, #CN1
                             self.NN2_wavg, #CN2
                             self.Dsupport, #Z
                             self.charge, #Bader charge
                             self.Nsites, #number of sites
                             self.PdC1, 
                             self.PdC2, 
                             self.PdC3]
                            
        

#%% Analyse the structures
Ntot = len(structures)

labels = ['Filename', 'AtomsObject', 'Eads', 'NPd', 'SiteType', 'RealSite', 
          'CN1', 'CN2', 'Z', 'Charge', 'Nsites', 'PdC1', 'PdC2', 'PdC3']
descriptors =  ['NPd', 'CN1', 'CN2', 'Z', 'Charge', 'Nsites', 'PdC1', 'PdC2', 'PdC3']
fdata = pd.DataFrame(columns = labels)


for i,struct in enumerate(structures):
    
    PdCO_ob = PdCO()
    PdCO_ob.get_descriptors(struct, data)
    fdata.loc[i,:] = PdCO_ob.structureID
fdata.to_csv('descriptor_data.csv')

dem =  np.array(fdata.loc[:,descriptors], dtype = float)
Eads = np.array(fdata.loc[:,'Eads'], dtype = float)
filename_list = list(fdata.loc[:,'Filename'])
sitetype_list = list(fdata.loc[:,'SiteType'])

pickle.dump([dem, Eads, descriptors], open('pca_data.p','wb'))

