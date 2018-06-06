#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 18:28:09 2018

@author: wangyifan
"""
import Cluster_structure as cs
import networkx as nx


A = '1-2 1-3 2-3 4-2 4-3 5-4 5-3 '
edge_list = cs.neighbor_string_list(A) 
wt_i = cs.neighbor_weight(edge_list)