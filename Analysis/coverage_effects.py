# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 14:02:48 2018

@author: wangyf
"""
import sys
import os
import numpy as np
import pandas as pd

# working on all platforms, using cross-platform home path

HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Generator'))
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Analysis'))
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Generator','user_inputs'))

import process_specnum as pspec
import plot_specnum as dspec
from IO_Pdn import *
from nc import *


data = pd.read_clipboard() 

################## User input ##################################
T  = 300
xtruncate = 0


output = 'analysis_results'

Parent_dir = os.getcwd()
output_dir = os.path.join(Parent_dir, output)
if not os.path.exists(output_dir): os.makedirs(output_dir)
csv_dir =  os.path.join(output_dir, str(T)+'K_profile.csv')
data.to_csv(csv_dir, sep = '\t')


#%% Plot the results

data = pd.read_csv(csv_dir, sep = '\t')
cov_vecs = np.array([0.005, 	0.01,	0.06,	0.12,	0.18])
s_spec_name = ['Pd1*', 'Pd2*', 'Pd3*', 'Pd4*']
s_spec = np.array(data[s_spec_name].T)

d3_spec_name = ['Pd4_3&1*', 'Pd5_4&1*', 'Pd6_4&2*', 'Pd6_5&1*', 'Pd7_5&2*', 'Pd7_6_1&1*', 'Pd7_6_2&1*','Pd7_6_3&1*']
d3_spec = np.array(data[d3_spec_name].T)

l_spec_name = ['Pd5*', 'Pd6_1*',	'Pd6_2*',	'Pd6_3*', 'Pd7_1*' ,	'Pd7_2*', 'Pd7_3*', 	'Pd7_4*']
l_spec = np.array(data[l_spec_name].T)

xlabel = 'Initial '+ r'$Pd_{1}$'+ ' Coverage (ML)'
ylabel = 'Species count'

if xtruncate == 0:
    xrange = []
else: 
    xrange = xrange

dspec. PlotTimeSeries(cov_vecs, s_spec, 
                      xlab = xlabel, ylab = ylabel, xlimit = xrange,
                      series_labels = s_spec_name,
                      fname = os.path.join(output_dir, 'cov_small_clusters.png'))

dspec. PlotTimeSeries(cov_vecs, l_spec, 
                      xlab = xlabel, ylab = ylabel, xlimit = xrange,
                      series_labels = l_spec_name,
                      fname = os.path.join(output_dir, 'cov_large_clusters.png'))

dspec. PlotTimeSeries(cov_vecs, d3_spec, 
                      xlab = xlabel, ylab = ylabel, xlimit = xrange,
                      series_labels = d3_spec_name,
                      fname = os.path.join(output_dir, 'cov_3d_clusters.png'))