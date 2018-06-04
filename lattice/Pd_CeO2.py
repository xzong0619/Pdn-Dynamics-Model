#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 18:02:30 2018

@author: wangyifan
"""

#Lattice constant
a = 5.49


#Coordinates of the unit-cell atoms in fractional 
from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.lattice.surface import hcp0001
from ase.io import write,read
from ase.build import surface,mx2, fcc111, add_adsorbate
from ase.visualize import view
import numpy as np
import math
from numpy import *

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
slab = slab.repeat((1,1,1))
slab.center(vacuum=10, axis=2)
del slab[[atom.index for atom in slab if atom.z>15]]
view(slab)
#Get the position of each atom
p = np.array(slab.get_positions())
#finds the coordinates of Rh atom, read from figure
#Pd_1 = [1.941,7.844]
#Pd_2= [0.00,4.483]
#Adds single Rh atoms to the surface
#add_adsorbate(slab,'Pd',0,Pd_1)
#add_adsorbate(slab,'Pd',0,Pd_2)

view(slab)
write('slab.png' , slab)
#Get the position of each atom
p = slab.get_positions()
#print(p)


sites = np.array(np.matrix([(p[-4]+p[-1])/2, #O bridge
							(p[-3]+p[-2])/2,
							(p[-3]+p[-4])/2,
							(p[-1]+p[-2])/2, 							
							p[-4],  #O atop
							p[-3],
							p[-2],
							p[-1]]))
#Get x and y coordinates of sites 
cart_coords_list = np.delete(sites,[2],axis=1) 
print(cart_coords_list)
#Calculates distance between 2O atoms
x = math.sqrt(abs(sum(np.multiply(p[-1]-p[-2],p[-1]-p[2])))/3)
print(x)