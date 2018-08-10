# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 11:08:25 2018

@author: wangyf
"""

from config import mother
import numpy as np
import lat_fun as lf
from itertools import combinations


index = np.arange(len(mother))

c1 = list(combinations(index,1))

c2 = list(combinations(index,2))

c3 = list(combinations(index,3))

'''
use c2 as an example
'''
d2 = []
for i in range(len(c2)):
    a = c2[i][0]
    b = c2[i][1]
    d2.append((lf.two_points_D(mother[a], mother[b]), (mother[a][2], mother[b][2])))
list(set(d2))