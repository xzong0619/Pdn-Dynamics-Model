# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 17:00:42 2018

@author: yifan
"""
import numpy as np
from ase.io import read, write
from ase.visualize import view
from ase import Atom
from sympy import nsolve
from sympy.abc import x,y,z
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