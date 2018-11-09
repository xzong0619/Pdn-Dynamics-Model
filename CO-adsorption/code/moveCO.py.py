# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:30:21 2018

@author: wangyf
"""
#%%
from ase.io import read, write
from ase.visualize import view
from ase.data import covalent_radii
import numpy as np
import os
from math import sin, cos,radians


Pdr = covalent_radii[46]
Or = covalent_radii[8]
Cr = covalent_radii[6]
CO = Cr + Or - 0.3
PdO = Pdr+Or -0.3
#%%
old_name = 'pd3-ceria-co-CONTCAR'
atoms = read(old_name)
view(atoms)

nPd = 0

for i, atom in enumerate(atoms):
    if atom.symbol == 'Pd': 
        nPd = nPd + 1
        
for j, atom in enumerate(atoms):
    if atom.symbol == 'C':
        Cj = j
        
Dist  = []
Ok  = []       
 
for k, atom in enumerate(atoms):
    if atom.symbol == 'O':
        dist = atoms.get_distance(Cj, k)
        Dist.append(dist)
        Ok.append(k)

Ok = Ok[Dist.index(min(Dist))]


#%%
Pdi = [98, 102, 103]
nsite = len(Pdi)

if nsite == 1: sitetype = 't'
if nsite == 2: sitetype = 'b'
if nsite == 3: sitetype = 'h'

Pdpos = []
for i in Pdi:
    Pdpos.append(atoms[i].position)

if nsite == 1: #top site
    phi = 0
    theta = 0
    rotate = [sin(radians(phi))*cos(radians(theta)),  
              sin(radians(phi))*sin(radians(theta)), 
              cos(radians(phi))]
    atoms[Cj].position = Pdpos[0] + np.array([rotate])*PdO
    atoms[Ok].position = atoms[Cj].position + np.array([rotate])*CO

if nsite == 2 or nsite == 3: #bridge or hollow site
    phi = 30
    theta = -60
    rotate = [sin(radians(phi))*cos(radians(theta)),  
              sin(radians(phi))*sin(radians(theta)), 
              cos(radians(phi))]
    translate = [0.8, 0.8, 0.5] 
    atoms[Cj].position = np.mean(Pdpos, axis = 0) + (np.array([translate])*np.array([rotate])) *PdO
    atoms[Ok].position = atoms[Cj].position + (np.array([translate])* np.array([rotate])) *CO   
    
view(atoms)


#%% Save the atom object
index = 1
Base_path = os.getcwd()
filename = 'Pd'+str(nPd) + '-' + sitetype + '-' + str(index) + '-CONTCAR'
output_dir = os.path.join(Base_path, 'outputs')
if not os.path.exists(output_dir): os.makedirs(output_dir)

write(os.path.join(output_dir, filename), atoms)
