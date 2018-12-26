# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 17:00:42 2018

@author: yifan
"""
import sys
import os
import numpy as np
import pandas as pd


from ase.io import read, write
from ase.visualize import view
from ase import Atom
from ase.data import covalent_radii
from sympy import nsolve
from sympy import Plane, Point3D
from sympy.abc import x,y,z
from itertools import combinations 
from math import sin, cos,radians

HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model','CO-adsorption', 'pool'))

from Pdbulk import NN1,NN2

#%% User define function block

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
    
    #Example: remove_CO('pd20-ceria-co-CONTCAR')

def sort_i_and_d(D,I):
    '''
    Sort indices based on the atomic bond length
    return the sorted bond length and indices in ascending order
    '''
    Dsort = np.sort(D)
    Dsort = list(Dsort)
    D = list(D)
    Isort = []
    for d in Dsort:
        Isort.append(I[D.index(d)])
        
    return Dsort,Isort

def find_all_Pd(atoms):
    '''
    Count number of atoms, return all Pd atom indices
    '''
    Pd_indices = []
    
    for i, atom in enumerate(atoms):
        if atom.symbol == 'Pd': Pd_indices.append(i)  
        
    return Pd_indices

def find_nearest_Pd(Pdi, Pd_indices):
    '''
    Input one target Pd atom index and all the indices for Pd atoms
    return the three relevent indices
    n = 3 for top sites
    '''
    #Get the bond length of C-Ce, Pd-C, Pd-Pd
    #and NN table for each Pd

    Pd_Pd = atoms.get_distances(Pdi, Pd_indices, mic = True) #all Pd-Pd bond length 
    Pd_Pd, Pd_indices = sort_i_and_d(Pd_Pd, Pd_indices) #sorted Pd-Pd bond length    
    Pd_sites = Pd_indices
    
    return Pd_sites

def find_nearest_Pd_top(Pdi, Pd_indices):
    
    Pd_sites = find_nearest_Pd(Pdi, Pd_indices)[:3]
    return Pd_sites

def find_nearest_Pd_bridge(Pdi, Pd_indices):
    '''
    Input two target Pd atom index and all the indices for Pd atoms
    return the three relevent indices
    '''
    Pd1 = Pdi[0]
    Pd2 = Pdi[1]
    
    Pd_sites1 = find_nearest_Pd(Pd1, Pd_indices)
    Pd_sites1 = [x for x in Pd_sites1 if x not in Pdi]
    Pd_sites2 = find_nearest_Pd(Pd2, Pd_indices)
    Pd_sites2 = [x for x in Pd_sites2 if x not in Pdi]
    
    Pd3 = [i for i in Pd_sites1 if i in Pd_sites2]
    
    Pd_sites = [Pd1, Pd2, Pd3[0]]
    
    return Pd_sites

def find_bridge_pairs(Pd_pairs, atoms):
    
    bridge_pairs = []
    for pair in Pd_pairs:
        Pd_Pd = atoms.get_distances([pair[0]], [pair[1]])
        if np.logical_and(Pd_Pd>=NN1[0], Pd_Pd<=NN1[1]):
            bridge_pairs.append(pair)
    
    return bridge_pairs
            
def find_hollow_triples(Pd_triples , atoms):
    
    hollow_triples = []
    for triple in Pd_triples:
        Pd_Pd1 = atoms.get_distances(triple[0], [triple[1], triple[2]])
        Pd_Pd2 = atoms.get_distances([triple[1]], [triple[2]])
        flag1 = np.logical_and(Pd_Pd1>=NN1[0], Pd_Pd1<=NN1[1])
        flag2 = np.logical_and(Pd_Pd2>=NN1[0], Pd_Pd2<=NN1[1])
      
        if np.all(list(flag1)+list(flag2)):
            hollow_triples.append(list(triple))
    
    return hollow_triples
    
def sphere_fun(pos, r):
    f = (x-pos[0])**2 +  (y-pos[1])**2 +  (z-pos[2])**2 - r**2
    return f

def find_CO_pos(Pdi, atoms, sitetype):

    Pdpos = []
    for i in Pdi: Pdpos.append(atoms[i].position)
    PdC = site_PdC[sitetype]
    

    left_side_indcies = [96, 98, 99, 106, 107, 110, 111, 112, 113, 115]
    phi = 60
    
    if sitetype == 'hollow':
        flags =[]        
        for i in Pdi:
            flags.append(i in left_side_indcies)
        
    if sitetype == 'bridge':
        flags = [Pdi[0] in left_side_indcies, Pdi[1] in left_side_indcies]
        
    if sitetype == 'top':
        flags = Pdi[0] in left_side_indcies
    
    if np.all(flags): theta = 180 
    else: theta = 300
        
    rotate = [sin(radians(phi))*cos(radians(theta)),  
              sin(radians(phi))*sin(radians(theta)), 
              cos(radians(phi))]

    initial_guess = tuple(Pdpos[0] + np.array(rotate)*PdC)
#    if sitetype == 'top':
#        initial_guess = tuple(Pdpos[0])
    f1 = sphere_fun(Pdpos[0],PdC[0])
    f2 = sphere_fun(Pdpos[1],PdC[1])
    f3 = sphere_fun(Pdpos[2],PdC[2])

    CO_pos = np.array(list(nsolve((f1, f2, f3), (x,y,z), initial_guess))).astype(float)
    
    return CO_pos
    
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


#%%
'''
Bond length value calculation from the given dataset
'''

csv_file  = os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model','CO-adsorption', 'pool', 'descriptor_data.csv')
fdata = pd.read_csv(csv_file)
descriptors =  ['NPd', 'CN1', 'CN2','GCN', 'Z', 'Charge', 'Nsites', 'Pd1C', 'Pd2C', 'Pd3C'] #10 in total
sitetype_list = list(fdata.loc[:,'SiteType'])
X =  np.array(fdata.loc[:,descriptors], dtype = float)

PdCs = ['Pd1C', 'Pd2C', 'Pd3C'] 
PdCi = []

for PdC in PdCs: PdCi.append(descriptors.index(PdC))
Xsite = []   
for cnt in PdCi:
    for site in ['top', 'bridge', 'hollow']:
        indices = np.where(np.array(sitetype_list) == site)[0]
        Xsite.append(X[:,cnt][indices])

top_PdC1 = np.mean(Xsite[0])
bridge_PdC1 = np.mean(Xsite[1])
hollow_PdC1 = np.mean(Xsite[2])

top_PdC2 = np.mean(Xsite[3][np.nonzero(Xsite[3])])
bridge_PdC2 = bridge_PdC1 #np.mean(Xsite[4][np.nonzero(Xsite[4])])
hollow_PdC2 = hollow_PdC1 #np.mean(Xsite[5][np.nonzero(Xsite[5])])

top_PdC3 = np.mean(Xsite[6][np.nonzero(Xsite[6])])
bridge_PdC3 = np.mean(Xsite[7][np.nonzero(Xsite[7])])
hollow_PdC3 = hollow_PdC1 #np.mean(Xsite[8][np.nonzero(Xsite[8])])

site_PdC = dict([('top', [top_PdC1, top_PdC2, top_PdC3]),
                 ('bridge', [bridge_PdC1, bridge_PdC2, bridge_PdC3]),
                 ('hollow', [hollow_PdC1, hollow_PdC2, hollow_PdC3])])
#%%
    
'''
Load the CONTCAR file with Pdn structure
'''
    
filename = 'pd20-test-CONTCAR'
atoms = read(filename)
view(atoms)

#Find all Pd atoms
Pd_indices =  find_all_Pd(atoms)
Pd_no_interest = [97, 100, 102, 103]
Pd_interest = [x for x in Pd_indices if x not in Pd_no_interest]

#Find all CO adsorption sites
top_sites = []
for Pdi in Pd_interest:
     top_sites.append(find_nearest_Pd_top(Pdi, Pd_indices))
top_sites = [[96, 99, 98],
            [98, 99, 110],
            [99, 96, 98],
            [101, 110, 105],
            [104, 108, 105],
            [105, 104, 101],
            [106, 113, 100],
            [107, 106, 100],
            [108, 114, 104],
            [109, 104, 108],
            [110, 112, 101],
            [111, 110, 101],
            [112, 99, 110],
            [113, 99, 112],
            [114, 105, 112],
            [115, 112, 113]]

bridge_sites = []     
Pd_pairs  = list(combinations(Pd_interest,2))
bridge_pairs = find_bridge_pairs(Pd_pairs, atoms)
for Pdi in bridge_pairs:
     bridge_sites.append(find_nearest_Pd_bridge(Pdi, Pd_indices))
     

bridge_sites = [[107,96,106], [96,98,99],[98,111,110], [106,99,113],[99,110,112],[112,113,115],
                [107,106,96], [106,113,99], [113,115,112],
                [96,99,106], [99,112,113],  [98,110,99],
                [96,106,99], [99,98,110], [99,113,112],
                [115,112,113], [112,110,99], [110,111,98],
                [115,114,112], [114,108,105], [108,109,104],
                [112,114,105], [110,105,112], [105,108,114], [111,101,110], [101,104,105], [104,109,108],
                [114,105,108], [105,101,104], [104,108,109],
                [112,105,114], [105,104,108], [110,101,105]]
        
hollow_sites = [] 
Pd_triples  = list(combinations(Pd_interest,3))
hollow_triples = find_hollow_triples(Pd_triples, atoms)
hollow_sites_no_interest = [[98, 101, 111], [98, 101, 110], [99, 105, 110],[99, 105, 112],[112, 113, 114], [113, 114, 115]]
for triple in hollow_triples:
    if triple not in hollow_sites_no_interest: hollow_sites.append(triple)            
            
            
#Find corresponding CO position
CO_pos_hollow = []
for site in hollow_sites:
    CO_pos_i_hollow = find_CO_pos(site, atoms, 'hollow')
    #atoms.append(Atom('C', CO_pos_i_hollow))
    CO_pos_hollow.append(CO_pos_i_hollow)

CO_pos_bridge = []
for site in bridge_sites:
    CO_pos_i_bridge = find_CO_pos(site, atoms, 'bridge')
    #atoms.append(Atom('C', CO_pos_i_bridge))
    CO_pos_bridge.append(CO_pos_i_bridge)
    
CO_pos_top = []
  
for site in top_sites:
    CO_pos_i_top = find_CO_pos(site, atoms, 'top')
    atoms.append(Atom('C', CO_pos_i_top))
    CO_pos_top.append(CO_pos_i_top)

#%%        



view(atoms)