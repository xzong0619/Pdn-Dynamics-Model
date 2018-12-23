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

atoms = read('pd1-ceria-co-CONTCAR')



view(atoms)


