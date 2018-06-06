# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import os
import sys
sys.path.append('/home/vlachos/wangyf/Zacros-Wrapper')
import zacros_wrapper as zw


'''
Scaldown
'''

# Read input files
x = zw.kmc_traj(path = '/home/vlachos/wangyf/Alumina/Network/12')
x.gas_prod = 'CO2'
x.exe_file = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
x.ReadAllInput()
final_data = zw.ReachSteadyStateAndRescale(x, '/home/vlachos/wangyf/Alumina/Network/12/outputs', n_runs = 16, n_batches = 1000, n_samples = 100, max_iterations = 10)