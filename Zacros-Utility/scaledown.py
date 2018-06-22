# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import os
import sys
sys.path.append('/home/vlachos/wangyf/programs/ZW_git/Zacros-Wrapper')
import zacros_wrapper as zw


zacros_exe = '/home/vlachos/wangyf/programs/ZW_git/Zacros-Wrapper/zacros_ZW.x'
KMC_source = os.getcwd()
output = 'sd_outputs'
BatchPath = os.path.join(KMC_source, output)


'''
Scaldown
'''

# Read input files
x = zw.kmc_traj(path = KMC_source)
x.gas_prod = 'CO'
x.exe_file = zacros_exe
x.ReadAllInput()
final_data = zw.ReachSteadyStateAndRescale(x, BatchPath, n_runs = 2, n_batches = 1000, n_samples = 100, max_iterations = 10)