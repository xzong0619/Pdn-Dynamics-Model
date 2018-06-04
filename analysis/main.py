# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 14:31:23 2018

@author: wangyf
"""

import sys
import os
#sys.path.append('/home/vlachos/wangyf/Zacros-Wrapper')
#import zacros_wrapper as zw

################## User input ##################################

KMC_source = os.getcwd()
output = 'outputs'
BatchPath = os.path.join(KMC_source, output)
print(KMC_source)
print(BatchPath)
n_runs = 10
