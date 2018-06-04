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
slab = slab.repeat((2,2,1))
slab.center(vacuum=10, axis=2)
del slab[[atom.index for atom in slab if atom.z>15]]


pov_args = {
	'transparent': True, #Makes background transparent. I don't think I've had luck with this option though
    'canvas_width': 1920., #Size of the width. Height will automatically be calculated. This value greatly impacts POV-Ray processing times
    'display': False, #Whether you want to see the image rendering while POV-Ray is running. I've found it annoying
    'rotation': '-50x', #Position of camera. If you want different angles, the format is 'ax, by, cz' where a, b, and c are angles in degrees
    'celllinewidth': 0.02, #Thickness of cell lines
    'show_unit_cell': 0 #Whether to show unit cell. 1 and 2 enable it (don't quite remember the difference)
    #You can also color atoms by using the color argument. It should be specified by an list of length N_atoms of tuples of length 3 (for R, B, G)
    #e.g. To color H atoms white and O atoms red in H2O, it'll be:
    #colors: [(0, 0, 0), (0, 0, 0), (1, 0, 0)]
    }

#%%
counter = 1
Pdz = 15
#initial positions
Pd1i = (slab[17].position + slab[16].position)/2
Pd1i[2] = Pdz
Pd2i = (slab[37].position + slab[39].position)/2
Pd2i[2] = Pdz
Pd3i = (slab[58].position + slab[79].position)/2
Pd3i[2] = Pdz
Pd4i = (slab[57].position + slab[59].position)/2
Pd4i[2] = Pdz
slab.append(Atom('Pd',Pd1i))
slab.append(Atom('Pd',Pd2i))
slab.append(Atom('Pd',Pd3i))
slab.append(Atom('Pd',Pd4i))

#Pd4 before hopping
Pd3b = (slab[37].position + slab[36].position)/2
Pd4b = (slab[78].position + slab[79].position)/2
Pd4c = (slab[77].position + slab[79].position)/2
Pd4d = np.array((15.528, 2.241, Pdz))

#final positions
Pd1f = np.array((11.646,2.241,Pdz))
Pd2f = np.array((9.705, 5.603, Pdz))
Pd3f = np.array((13.587,5.603,Pdz))
Pd4f = np.array((11.646,4.483,Pdz+1.7))
flname = 'slab'+str(counter)
write(flname+'.png', slab, rotation = '-50x')
write(flname+'.POV', slab, **pov_args)

#%% Pd1+ Pd2
n = 20
dPd1 = (Pd1f-Pd1i)/n
dPd2 = (Pd2f-Pd2i)/n
dPd3 = (Pd3b-Pd3i)/n
dPd4 = (Pd4b-Pd4i)/n

Pd1now = Pd1i
Pd2now = Pd2i
Pd3now = Pd3i
Pd4now = Pd4i

for i in range(n):
    counter = counter +1
    del slab[[-1,-2, -3, -4]]
    
    Pd1now = dPd1 + Pd1now
    Pd2now = dPd2 + Pd2now
    Pd3now = dPd3 + Pd3now
    Pd4now = dPd4 + Pd4now
    slab.append(Atom('Pd',Pd1now))
    slab.append(Atom('Pd',Pd2now))
    slab.append(Atom('Pd',Pd3now))
    slab.append(Atom('Pd',Pd4now))
    flname = 'slab'+str(counter)
    write(flname+'.png', slab, rotation = '-50x')
    write(flname+'.POV', slab, **pov_args)

#%% Pd3+ Pd1/2
n2 = 20

dPd3 = (Pd3f-Pd3b)/n2
dPd4 = (Pd4c-Pd4b)/n2


for i in range(n2):
    counter = counter +1
    del slab[[-1,-2]]
    Pd3now = dPd3 + Pd3now
    Pd4now = dPd4 + Pd4now
    slab.append(Atom('Pd',Pd3now))
    slab.append(Atom('Pd',Pd4now))
    flname = 'slab'+str(counter)
    write(flname+'.png', slab, rotation = '-50x')   
    write(flname+'.POV', slab, **pov_args)

#%% Pd4 diffusion
n3 = 20
dPd4 = (Pd4d - Pd4c)/n3
for i in range(n3):
    counter = counter + 1
    del slab[[-1]]
    Pd4now = dPd4 + Pd4now
    slab.append(Atom('Pd',Pd4now))
    flname = 'slab'+str(counter)
    write(flname + '.png', slab,  rotation = '-50x')
    write(flname+'.POV', slab, **pov_args)
      

#%% hopping
nh = 10
dPd4h = (Pd4f - Pd4d)/nh
for i in range(nh):
    counter = counter + 1
    del slab[[-1]]
    Pd4now = dPd4h + Pd4now
    slab.append(Atom('Pd',Pd4now))
    flname = 'slab'+str(counter)
    write(flname + '.png', slab,  rotation = '-50x')
    write(flname+'.POV', slab, **pov_args)
    
    







'''
#Write to POV-Ray file

write('step_Pt.POV', slab, **pov_args)
'''