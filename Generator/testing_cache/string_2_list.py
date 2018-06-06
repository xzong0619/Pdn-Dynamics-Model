#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 16:09:59 2018

@author: wangyifan
"""
import re
A = '1-2 2-3'
p = re.compile(r'\W+')
B  = p.split(A)
ll = []
nb = int(len(B)/2)
for i in range(nb):
    bi = 2* i 
    ll.append((int(B[bi]), int(B[bi+1])))

#def neighbor_string_list(string):
