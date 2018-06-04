import sys
import os
sys.path.append('/home/vlachos/wangyf/programs/ase')
sys.path.append('/home/vlachos/wangyf/Zacros-Wrapper/zacros_wrapper')
from Lattice import Lattice as lat
#Creats CeO2(111) surface and place a Rh atom on three-fold O hollow sites
#Identifies site types and generates lattice_input.dat for zacros

#CeO2(111) surface

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
slab = slab.repeat((2,2,1))
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
cell = slab.get_cell()[0:2, 0:2]/2
cell_x = abs(cell[0,0] - cell[0,1])
p = np.array(p)
print(cell_x)

sites = np.array(np.matrix([p[56],  #O atop
							p[57],
							p[58],
							p[59]]))
sites[:,0] = sites[:,0] - cell_x
#Get x and y coordinates of sites 
cart_coords_list = np.delete(sites,[2],axis=1) 
print(cart_coords_list)
#Calculates distance between 2O atoms
x = math.sqrt(abs(sum(np.multiply(p[56]-p[59],p[56]-p[59])))/3)
print(x)

# Set up object KMC lattice
KMC_lat = lat()
KMC_lat.text_only = False 
workingdir = '.'
KMC_lat.lattice_matrix = cell
KMC_lat.site_type_names = ['T']
KMC_lat.site_type_inds = [1,1,1,1]
KMC_lat.set_cart_coords(cart_coords_list)

#Select the nearest neighbors
KMC_lat.Build_neighbor_list(cut = x+x)

# write lattice_input.dat
KMC_lat.Write_lattice_input(workingdir)   
# Create a png file with the lattice drawn
plt = KMC_lat.PlotLattice(plot_neighbs = True)
plt.savefig(os.path.join('.', 'Pd_kmc_lattice_b.png'))
plt.close()