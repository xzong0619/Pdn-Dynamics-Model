# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 14:15:41 2018

@author: yifan
"""

from ase.io import read, write
from ase.visualize import view


nPd = 20
sitetype = 'h'
index = 1
filename = 'Pd'+str(nPd) + '-' + sitetype + '-' + str(index) + '-CONTCAR'

atoms = read(filename)
view(atoms)