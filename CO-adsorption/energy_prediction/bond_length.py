# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 17:00:42 2018

@author: yifan
"""
import numpy as np
import pandas as pd

from ase.io import read, write
from ase.visualize import view
from ase import Atom
from sympy import nsolve
from sympy.abc import x,y,z

from Pdbulk import NN1,NN2
from sympy import Plane, Point3D
#%%
'''
Bond length value calculation from the given dataset
'''
PdCs = ['Pd1C', 'Pd2C', 'Pd3C'] 
PdCi = []

for PdC in PdCs: PdCi.append(descriptors.index(PdC))
Xsite = []   
for cnt in PdCi:
    for site in ['top', 'bridge', 'hollow']:
        indices = np.where(np.array(sitetype_list) == site)[0]
        Xsite.append(X[:,cnt][indices])

top_PC1 = np.mean(Xsite[0])
bridge_PC1 = np.mean(Xsite[1])
hollow_PC1 = np.mean(Xsite[2])

top_PC2 = np.mean(Xsite[3][np.nonzero(Xsite[3])])
bridge_PC2 = np.mean(Xsite[4][np.nonzero(Xsite[4])])
hollow_PC2 = np.mean(Xsite[5][np.nonzero(Xsite[5])])

top_PC3 = np.mean(Xsite[6][np.nonzero(Xsite[6])])
bridge_PC3 = np.mean(Xsite[7][np.nonzero(Xsite[7])])
hollow_PC3 = np.mean(Xsite[8][np.nonzero(Xsite[8])])


#%%

def remove_CO(filename):
    '''
    Read the old CONTCAR file with a CO onto it
    remove the CO and save as a new CONTCAR
    '''
    
    old_name = filename #'pd20-ceria-co-CONTCAR'
    atoms = read(old_name)
    view(atoms)
    
    # find number of Pd
    # find C atom index
    nPd = 0
    
    for i, atom in enumerate(atoms):
        if atom.symbol == 'Pd': 
            nPd = nPd + 1
        if atom.symbol == 'C':
            C_in_CO = i
            
    C_O_Dist  = []
    O_in_CO  = []       
     
    for k, atom in enumerate(atoms):
        if atom.symbol == 'O':
            dist = atoms.get_distance(C_in_CO, k)
            C_O_Dist.append(dist)
            O_in_CO.append(k)
    
    O_in_CO = O_in_CO[C_O_Dist.index(min(C_O_Dist))]
    
    del atoms[[O_in_CO, C_in_CO]]
    view(atoms)
    write('pd'+str(nPd)+'-test-CONTCAR', atoms)
    
#remove_CO('pd20-ceria-co-CONTCAR')
    
'''
Load the CONTCAR file with Pdn structure
'''
    
filename = 'pd20-test-CONTCAR'
atoms = read(filename)
view(atoms)

#%%

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
        
        self.sitetype  = []
        self.descriptors = []
        
    def get_descriptors(self, atoms, CI_PdC = 0.2):
        '''
        Takes in a structure dictionary with an atoms object and filename 
        Input tolerance for Pd-C bond 0.3 A for detecting CO ads site
        '''
        self.atoms = atoms
        
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
        Pd1NN= dict() #Pd NN table 
        
        
        for i in Pdi:
           PdD =  self.atoms.get_distances(i, Pdi)
           PdD, Pdisort = sort_i_and_d(PdD,Pdi)
           
           Pd_Pd['Pd'+str(i)] =  PdD
           Pd_Pd['i'+str(i)] = Pdisort
           PdNN['Pd'+str(i)] = [sum(np.logical_and(PdD>=NN1[0], PdD<=NN1[1])),
                sum(np.logical_and(PdD>=NN2[0], PdD<=NN2[1]))]
           Pd1NN['Pd'+str(i)] = np.array(Pdisort)[np.where(np.logical_and(PdD>=NN1[0],PdD<=NN1[1]))[0]]
           
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
        Pd_C_CO = np.array(Pd_C[:len(COsites)]) #CO distance to neighboring Pd atoms
        
        '''
        Weighted average for NN1, NN2
        '''
        norm_weights = (1/Pd_C_CO)/np.sum(1/Pd_C_CO) #weights based on 1 over CO-Pd distance
        
        self.CN1 = np.dot(norm_weights, PdNN_CO.loc['NN1'].values)
        self.CN2 = np.dot(norm_weights, PdNN_CO.loc['NN2'].values)
        
        '''
        GCN calculation
        '''
        cn_max = [12, 18, 22]
        gcn_sum = 0
        for i in COsites:
            #Pd1NN['Pd'+str(i)] is a list of cn number for 1NN
            for j in Pd1NN['Pd'+str(i)]:
                gcn_sum  = gcn_sum + PdNN.loc['NN1']['Pd'+str(j)]
        
        self.GCN = gcn_sum/cn_max[self.Nsites -1]
        
        
        '''
        Make a row in dataframe as an ID for each structure including filenames and properties etc
        '''
        self.structureID =  [self.atoms, # atoms object
                             self.Eads, #Eads
                             self.NPd, #Npd
                             self.sitetype, #sitetype from calculation
                             self.CN1, #CN1
                             self.CN2, #CN2
                             self.GCN, # general cooridination number
                             self.Dsupport, #Z
                             self.Nsites, #number of sites
                             self.PdC1, #1st Pd-C distance 
                             self.PdC2, #2nd Pd-C distance
                             self.PdC3] #3rd Pd-C distance


Pdi = [115, 113, 112]
nsite = len(Pdi)

if nsite == 1: sitetype = 't'
if nsite == 2: sitetype = 'b'
if nsite == 3: sitetype = 'h'

Pdpos = []
for i in Pdi:
    Pdpos.append(atoms[i].position)
    
def sphere_fun(pos, r):
    f = (x-pos[0])**2 +  (y-pos[1])**2 +  (z-pos[2])**2 - r**2
    return f


f1 = sphere_fun(Pdpos[0],bridge_PC1)
f2 = sphere_fun(Pdpos[1],bridge_PC2)
f3 = sphere_fun(Pdpos[2],bridge_PC3)
initial_guess =  tuple(Pdpos[0])
CO_pos = np.array(list(nsolve((f1, f2, f3), (x,y,z), initial_guess))).astype(float)
#atoms[C_in_CO].position = CO_pos

atoms.append(Atom('C', CO_pos))
view(atoms)