# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 11:27:11 2018

@author: wangyf
"""

from ase.io import read, write
from ase.visualize import view
from ase.data import covalent_radii
from ase import Atoms, Atom
from ase.build import surface
from structure_constants import mother
from test_connectivity_fitness import occ_nodes

Pdr = covalent_radii[46]
Or = covalent_radii[8]
PdPd = Pdr*2
PdO = Pdr + Or
PdOcc = mother[occ_nodes]
OO = 3.882

def ceria():
    #Lattice constant
    a = 5.49
    CeO2 = Atoms('Ce4O8', scaled_positions =[ (0., 0., 0.),
                  (0., 0.5, 0.5),
                  (0.5, 0., 0.5),
                  (0.5, 0.5, 0.),
                  (0.75, 0.25, 0.25),
                  (0.25, 0.75, 0.75),
                  (0.75, 0.75, 0.75),
                  (0.25, 0.25, 0.25),
                  (0.25, 0.25, 0.75),
                  (0.75, 0.75, 0.25),
                  (0.25, 0.75, 0.25),
                  (0.75, 0.25, 0.75)],
                  cell = [a,a,a],
    			  pbc = True	)
    #Scales the atomic positions with the unit cell
    #CeO2.set_cell(cell, scale_atoms=True)
    #(1,1,1) is the slab type. There are 2 unit cells along 
    #z direction
    slab = surface(CeO2, (1, 1, 1), 2)
    
    #Repeating the slab 2 unit cells in x and 1 unit cell 
    #in y directions
    slab = slab.repeat((3,3,1))
    slab.center(vacuum=10, axis=2)
    del slab[[atom.index for atom in slab if atom.z>15]]
    
    return slab

support = ceria()
#view(support)

#%%
origin = support[89].position
origin[2] = origin[2] +  PdO
Pdm = PdOcc*PdPd + origin
#PdNP = Atoms('Pd36', positions = Pdm)
nPd = len(PdOcc)
for i in range(nPd):
    support.append(Atom('Pd', position = Pdm[i]))
view(support) 
