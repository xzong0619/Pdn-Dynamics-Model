# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 17:33:26 2018

@author: wangyf
"""

import sys
import os

# working on all platforms, using cross-platform home path

HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Cluster-Graphs'))

import Cluster_structure as cs

Descriptors = ['1',
               '1-2',
               '1-2 2-3',
               '1-2 1-3 2-3',
               '1-2 2-3']