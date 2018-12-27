# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 12:37:09 2018

@author: yifan
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 17:00:42 2018

@author: yifan
"""
import sys
import os
import numpy as np
import pandas as pd
import pickle 

from ase.io import read, write
from ase.visualize import view
from ase import Atom
from sympy import nsolve
from sympy import Plane, Point3D
from sympy.abc import x,y,z
from itertools import combinations 
from math import sin, cos,radians

from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler 

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
    D = list(D)
    Dsort = np.sort(D)
    D_unqiue_sort = np.sort(np.unique(D))
    D_unqiue_sort = list(D_unqiue_sort)
    Isort = []
    for d in D_unqiue_sort:
        indices = [i for i, x in enumerate(D) if x == d]
        for i in indices: Isort.append(I[i])
        
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
            bridge_pairs.append(list(pair))
    
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
        self.site_PdC = site_PdC
        
    def get_descriptors(self, atoms, CO_sites):
        '''
        Takes in a structure dictionary with an atoms object and filename 
        Input tolerance for Pd-C bond 0.3 A for detecting CO ads site
        '''
        self.atoms = atoms.copy()
        self.CO_sites = CO_sites 

        '''
        Determine C position and site type
        '''
        Pd_pos = []
        for i in self.CO_sites: Pd_pos.append(self.atoms[i].position)
        self.C_pos = np.mean(Pd_pos, axis = 0)
        self.atoms.append(Atom('C', self.C_pos))
        self.Nsites = len(self.CO_sites)
        
        if self.Nsites== 3: self.sitetype = 'hollow'                     
        if self.Nsites == 2: self.sitetype = 'bridge'
        if self.Nsites == 1: self.sitetype = 'top'
        
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
        
        PdC3 = self.site_PdC[self.sitetype]
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
        
        '''
        Detect CO adsorption sites
        '''
               
        COsites_cols = []
        for s in range(len(self.CO_sites)):
            COsites_cols.append('Pd'+str(self.CO_sites[s]))
        
        PdNN_CO = PdNN.loc[:, COsites_cols] #NN dataframe 
        #Pd_C_CO = np.array([:len(self.CO_sites)]) #CO distance to neighboring Pd atoms
        
        '''
        Average for NN1, NN2
        '''
        #norm_weights = (1/Pd_C_CO)/np.sum(1/Pd_C_CO) #weights based on 1 over CO-Pd distance
        
        self.CN1 = np.mean(PdNN_CO.loc['NN1'].values)
        self.CN2 = np.mean(PdNN_CO.loc['NN2'].values)
        
        
        '''
        Make a row in dataframe as an ID for each structure including filenames and properties etc
        '''
        self.structureID =  [self.atoms, # atoms object
                             self.Eads, #Eads
                             self.NPd, #Npd
                             self.sitetype, #sitetype from calculation
                             self.CO_sites, #Pd index at the adsoprtion site
                             self.CN1, #CN1
                             self.CN2, #CN2
                             #self.GCN, # general cooridination number
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
for Pdi in Pd_interest: top_sites.append([Pdi])

bridge_sites = []     
Pd_pairs  = list(combinations(Pd_interest,2))
bridge_sites = find_bridge_pairs(Pd_pairs, atoms)

        
hollow_sites = [] 
Pd_triples  = list(combinations(Pd_interest,3))
hollow_triples = find_hollow_triples(Pd_triples, atoms)
hollow_sites_no_interest = [[98, 101, 111], [98, 101, 110], [99, 105, 110],[99, 105, 112],[112, 113, 114], [113, 114, 115]]
for triple in hollow_triples:
    if triple not in hollow_sites_no_interest: hollow_sites.append(triple)            
   
CO_sites = top_sites + bridge_sites + hollow_sites

#
#%% Analyse the structures and save into a csv file
Ntot = len(CO_sites)        
labels = ['AtomsObject', 'Eads', 'NPd', 'SiteType', 'CO sites',
          'CN1', 'CN2', 'Z', 'Nsites', 'Pd1C', 'Pd2C', 'Pd3C']

fdata = pd.DataFrame(columns = labels)
for i,site in enumerate(CO_sites):
    
    PdCO_ob = PdCO()
    PdCO_ob.get_descriptors(atoms, site)
    #if PdCO_ob.filename != 'pd5-ceria-co-CONTCAR':
    fdata.loc[i,:] = PdCO_ob.structureID
fdata.to_csv('g_descriptor_data.csv', index=False, index_label=False)

#%% Predict the energy 

descriptors_g =  ['CN1', 'Z', 'Nsites',   'Pd1C', 'Pd2C', 'Pd3C'] #6 geometric descriptors
dem =  np.array(fdata.loc[:,descriptors_g], dtype = float)
sitetype_list = list(fdata.loc[:,'SiteType'])

pca = PCA()    
X_g = dem.copy()
X_std_g = StandardScaler().fit_transform(X_g)
Xpc_g = pca.fit_transform(X_std_g) 

estimator_file  = os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model','CO-adsorption', 'pool', 'g_estimator')
pcg_estimator = pickle.load(open(estimator_file,'rb'))
y_pcg = pcg_estimator.predict(Xpc_g)


#%% Plot the energy value for top, bridge, hollow sites




















