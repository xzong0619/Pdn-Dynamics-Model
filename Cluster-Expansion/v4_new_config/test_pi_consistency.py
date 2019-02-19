# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 17:21:13 2019

@author: wangyf
"""

import GA_functions as GA
import lattice_functions as lf

empty = 'grey'
filled = 'r'
occ = [empty, filled]
individual = index_to_one_hot([0])
Gsv = GA.individual_config(individual)
Cal = lf.calculations(occ)
Gcv = Gcv_nonzero

pi_pred =  Cal.get_pi_matrix(Gsv ,Gcv)