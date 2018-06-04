import sys
import os

#Creats CeO2(111) surface and place a Rh atom on three-fold O hollow sites
#Identifies site types and generates lattice_input.dat for zacros

#CeO2(111) surface

#Lattice constant
a = 5.49


#Coordinates of the unit-cell atoms in fractional 
from ase import Atoms
from ase import Atom
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
slab = slab.repeat((4,2,1))
slab.center(vacuum=10, axis=2)
del slab[[atom.index for atom in slab if atom.z>15]]

extra = [89,91,92,93,94,95,96,97,98,99,109,111,112,113,114,115,116,
         117,118,119,125,127,129,131,132,133,134,135,136,137,138,139,
         145,147,149,151,152,153,154,155,156,157,158,159,124,126,144,146, 140]
del slab[extra]


slab.append(Atom('Pd',(27.174,6.724,10.3)))
slab.append(Atom('Pd',(19.41,2.24,11.5)))
slab.append(Atom('Pd',(20,6.724,13.9)))
slab.append(Atom('Pd',(11.646,6.724,14.5)))
slab.append(Atom('Pd',(9.705,10.086,14.5)))
slab.append(Atom('Pd',(13.587,10.086,14.5)))
slab.append(Atom('Pd',(11.646,8.965,15.5)))

#Get the position of each atom
p = np.array(slab.get_positions())
#finds the coordinates of Rh atom, read from figure
#Pd_1 = [1.941,7.844]
#Pd_2= [0.00,4.483]
#Adds single Rh atoms to the surface
#add_adsorbate(slab,'Pd',0,Pd_1)
#add_adsorbate(slab,'Pd',0,Pd_2)

view(slab)

pov_args = {
	'transparent': True, #Makes background transparent. I don't think I've had luck with this option though
    'canvas_width': 7200., #Size of the width. Height will automatically be calculated. This value greatly impacts POV-Ray processing times
    'display': False, #Whether you want to see the image rendering while POV-Ray is running. I've found it annoying
    'rotation': '61.505x,-48.744y,-173.30z', #Position of camera. If you want different angles, the format is 'ax, by, cz' where a, b, and c are angles in degrees
    'celllinewidth': 0.02, #Thickness of cell lines
    'show_unit_cell': 0 #Whether to show unit cell. 1 and 2 enable it (don't quite remember the difference)
    #You can also color atoms by using the color argument. It should be specified by an list of length N_atoms of tuples of length 3 (for R, B, G)
    #e.g. To color H atoms white and O atoms red in H2O, it'll be:
    #colors: [(0, 0, 0), (0, 0, 0), (1, 0, 0)]
    }

#Write to POV-Ray file

write('step_Pt.POV', slab, **pov_args)
