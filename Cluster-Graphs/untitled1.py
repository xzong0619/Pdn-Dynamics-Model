#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 13:53:54 2018

@author: wangyifan
"""
import numpy as np
a = np.array([[1, 2], [3, 4]])
b = np.array([[5, 6]])
m = np.concatenate((a, b), axis=0)

n = np.concatenate((a, b.T), axis=1)
