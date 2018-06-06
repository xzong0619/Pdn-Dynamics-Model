# -*- coding: utf-8 -*-
"""
Created on Tue May 15 17:07:56 2018

@author: wangyf
"""
import itertools as iter

ni = len(H_list)
i_list = [] #  index of possible graphs including isophoric graphs'''
for i in range(ni):
    i_list.append(i)
index_pair = list(iter.combinations(i_list,2)) # combination of possible index pair
is_iso = [] # Array of whether each pair are isophoric, 1 = yes, None = no'''
for i in range(len(index_pair)):
    a = index_pair[i][0]
    b = index_pair[i][1]
    is_iso.append(H_list[a] == H_list[b])
for i in range(len(index_pair)):
    if is_iso[i] == 1:
        if index_pair[i][1] in i_list:
            i_list.remove(index_pair[i][1])
        #return i_list    
    