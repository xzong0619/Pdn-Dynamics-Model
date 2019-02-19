# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 17:21:13 2019

@author: wangyf
"""

import GA_functions as GA
import lattice_functions as lf
[Gm, Gsv, Gcv1, Gcv2, Gcv3, Gcv4] = pickle.load(open("clusters.p", "rb"))
Gcv = Gcv1+Gcv2+Gcv3


empty = 'grey'
filled = 'r'
occ = [empty, filled]
individual = GA.index_to_one_hot([0])
Gsv = GA.individual_config(individual)
Cal = lf.calculations(occ)


pi_pred =  Cal.get_pi_matrix(Gsv ,Gcv_nonzero)[0]
E_pred = float(np.dot(pi_pred, J_nonzero) + intercept)

pi_pred_all =  Cal.get_pi_matrix(Gsv ,Gcv)[0][J_index[1:]-1]
E_pred_all = float(np.dot(pi_pred_all, J_nonzero) + intercept)