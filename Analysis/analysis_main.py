# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 14:31:23 2018

@author: wangyf
"""

import sys
import os
import numpy as np

# working on all platforms, using cross-platform home path

HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Generator'))
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Analysis'))

import process_specnum as pspec
import plot_specnum as dspec

#%%
################## User input ##################################

output = 'analysis_results'
n_runs = 10
Count_to_Coverage = 1
lattice_dim = 25

Parent_dir = os.getcwd()
output_dir = os.path.join(Parent_dir, output)
if not os.path.exists(output_dir): os.makedirs(output_dir)

# Spec Object
s = pspec.read_Multiple_Spec(n_runs)
# Spec dataframe
s_df = s. multi_spec_ave_df
# Save dataframe as a CSV file
s_df.to_csv(os.path.join(output_dir, 'surf_spec.csv'), sep = '\t')

#%%
'''
# plot the evolution of species number/coverages with time
'''

time_vecs = np.array(s.t)
spec_name, n_spec = s.lifetime_spec

print(spec_name)
surf_spec_vecs = []
surf_spec_dent = [1,2,3,4,3]

for i in range(n_spec):
    surf_spec_vecs.append(np.array(s_df[spec_name[i]]))

if Count_to_Coverage == 0:
    ylabel = 'Species count'
else:
    ylabel = 'Coverage'
    # Convert to surface coverages    
    surf_spec_vecs = pspec.num_to_cov(lattice_dim, n_spec, surf_spec_dent, surf_spec_vecs)



xrange = (0,100)
dspec. PlotTimeSeries(time_vecs, surf_spec_vecs, 
                      xlab = 'Time (s)', ylab = ylabel, #xlimit = xrange,
                      series_labels = spec_name, 
                      fname = os.path.join(output_dir, 'surf_spec_vs_time.png'))

#%%
'''
# plot the distribution of surface species at the end state
'''
# end state species info
es_n = pspec.final_lt_spec_on_surf(s_df)

spec_cov = np.multiply(np.array(surf_spec_dent),np.array(es_n))

dspec.PlotPie(es_n, spec_name,  fname = os.path.join(output_dir, 'surf_spec_pie.png'))


