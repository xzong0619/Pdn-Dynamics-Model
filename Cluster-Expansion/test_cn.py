# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:25:35 2018

@author: wangyf
"""
'''
Use No.26 Pd10a to test
'''
from structure_constants import mother, config, Ec
import lattice_functions as lf
import pickle

[Gm, Gsv, Gcv1, Gcv2, Gcv3] = pickle.load(open("clusters.p", "rb"))

G = Gsv[25] #No.26 Pd10a [0,1,2,6,11,10,14,29,28,27]
COsite =  [2, 14, 27]

sitetype = len(COsite)

for i in range(sitetype):
    