# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 14:15:41 2018

@author: yifan
"""

from ase.io import read, write
from ase.visualize import view
import glob
files = []
for file in glob.glob("*CONTCAR"):
    files.append(file)


filename = files[0]

atoms = read('pd20-ceria-co-CONTCAR')



view(atoms)

pov_args = {
	'transparent': True, #Makes background transparent. I don't think I've had luck with this option though
    'canvas_width': 3600., #Size of the width. Height will automatically be calculated. This value greatly impacts POV-Ray processing times
    'display': False, #Whether you want to see the image rendering while POV-Ray is running. I've found it annoying
    'rotation': '0x,-60y,-90z', #Position of camera. If you want different angles, the format is 'ax, by, cz' where a, b, and c are angles in degrees
    'celllinewidth': 0.02, #Thickness of cell lines
    'show_unit_cell': 0 #Whether to show unit cell. 1 and 2 enable it (don't quite remember the difference)
    #You can also color atoms by using the color argument. It should be specified by an list of length N_atoms of tuples of length 3 (for R, B, G)
    #e.g. To color H atoms white and O atoms red in H2O, it'll be:
    #colors: [(0, 0, 0), (0, 0, 0), (1, 0, 0)]
    }

#Write to POV-Ray file

#write('Pd1CO.POV', atoms, **pov_args)

#del atoms[98]
#del atoms[96]
#write('Pd1.POV', atoms, **pov_args)

