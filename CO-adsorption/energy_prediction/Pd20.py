# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 16:24:52 2018

@author: yifan
"""

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
PdC = Pdr+Cr -0.3
#%%
'''
Read the old CONTCAR file with a CO onto it
'''

old_name = 'pd20-ceria-co-CONTCAR'
atoms = read(old_name)


nPd = 0

for i, atom in enumerate(atoms):
    if atom.symbol == 'Pd': 
        nPd = nPd + 1
        
for j, atom in enumerate(atoms):
    if atom.symbol == 'C':
        C_in_CO = j
        
C_O_Dist  = []
O_in_CO  = []       
 
for k, atom in enumerate(atoms):
    if atom.symbol == 'O':
        dist = atoms.get_distance(C_in_CO, k)
        C_O_Dist.append(dist)
        O_in_CO.append(k)

O_in_CO = O_in_CO[C_O_Dist.index(min(C_O_Dist))]


#%%
Pdi = [116]
nsite = len(Pdi)

if nsite == 1: sitetype = 't'
if nsite == 2: sitetype = 'b'
if nsite == 3: sitetype = 'h'

Pdpos = []
for i in Pdi:
    Pdpos.append(atoms[i].position)

if nsite == 1: #top site
    phi = 60
    theta = 180 #300
    rotate = [sin(radians(phi))*cos(radians(theta)),  
              sin(radians(phi))*sin(radians(theta)), 
              cos(radians(phi))]
    atoms[C_in_CO].position = Pdpos[0] + np.array([rotate])*PdC
    atoms[O_in_CO].position = atoms[C_in_CO].position + np.array([rotate])*CO

if nsite == 2 or nsite == 3: #bridge or hollow site
    phi = 30
    theta = -60
    rotate = [sin(radians(phi))*cos(radians(theta)),  
              sin(radians(phi))*sin(radians(theta)), 
              cos(radians(phi))]
    translate = [0.8, 0.8, 0.5] 
    atoms[C_in_CO].position = np.mean(Pdpos, axis = 0) + (np.array([translate])*np.array([rotate])) *PdC
    atoms[O_in_CO].position = atoms[C_in_CO].position + (np.array([translate])* np.array([rotate])) *CO   
    
view(atoms)