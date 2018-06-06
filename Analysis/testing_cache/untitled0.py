# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 10:58:06 2018

@author: wangyf
"""

from process_spec_num import *

fname =  'specnum_output.txt'
fldr = os.getcwd()
filepath = os.path.join(fldr, fname)

fid = open(filepath, 'r')
file = fid.read()
lines = file.splitlines()
dict_array = lines[0].lower().split()

data = lines[1:]
dict = {}
for x in range(0, len(dict_array)):
    dict[dict_array[x]] = x

spec = []
for s in data:
    spec.append(spec_info(s.split(),dict))