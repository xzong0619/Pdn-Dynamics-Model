# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 19:24:29 2018

@author: wangyf
"""

'''
store adsorption constants
'''

import pandas as pd

data = pd.read_clipboard()
data.to_csv('adsorption_constant.csv', index = False)



