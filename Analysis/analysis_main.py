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
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Generator','user_inputs'))

#%%
import process_specnum as pspec
import plot_specnum as dspec
from IO_Pdn import *
from nc import *
#%%
################## User input ##################################
output = 'analysis_results'
lattice_dim = 25

'''
Boolean parameters
'''
Count_to_Coverage = 0 #convert species count to coverage
ss_flag = 1 # analyze the data at steady state
xtruncate = 0 # truncate the data when plotting

ss_cut = 10**5# set steady range as the last ss_range seconds
xmax = 10**6 # set plotting range
xrange = (0,xmax) 

################## User input ##################################


Cluster = ClusterOut(nc)
Parent_dir = os.getcwd()
output_dir = os.path.join(Parent_dir, output)
if not os.path.exists(output_dir): os.makedirs(output_dir)

# Spec Object
s = pspec.read_Multiple_Spec()
# Spec dataframe
s_df = s. multi_spec_ave_df
# Save dataframe as a CSV file
s_df.to_csv(os.path.join(output_dir, 'surf_spec.csv'), sep = '\t')

#%%
'''
# plot the evolution of species number/coverages with time
'''

time_vecs = np.array(s_df['t'])
if xmax >= time_vecs[-1]: xtruncate = 0

# Analysis for steady state
if ss_cut >= time_vecs[-1]: ss_flag = 0
if ss_flag == 1:
    s_df = pspec.ss_search(ss_cut, s_df)
    time_vecs = np.array(s_df['t'])

spec_name, n_spec = pspec.lifetime_spec_on_surf(s_df)

print(spec_name)
surf_spec_vecs = []
surf_spec_dent = []
surf_spec_n = []

for i in range(n_spec):
    surf_spec_vecs.append(np.array(s_df[spec_name[i]]))
    surf_spec_dent.append(Cluster.surf_dent[Cluster.surf_spec.index(spec_name[i])])
    surf_spec_n.append(Cluster.surf_n[Cluster.surf_spec.index(spec_name[i])])

if Count_to_Coverage == 0:
    ylabel = 'Species count'
else:
    ylabel = 'Coverage'
    # Convert to surface coverages    
    surf_spec_vecs = pspec.num_to_cov(lattice_dim, n_spec, surf_spec_dent, surf_spec_vecs)

if xtruncate == 0:
    xrange = []
else: 
    xrange = xrange
    
dspec. PlotTimeSeries(time_vecs, surf_spec_vecs, 
                      xlab = 'Time (s)', ylab = ylabel, xlimit = xrange,
                      series_labels = spec_name, 
                      fname = os.path.join(output_dir, 'surf_spec_vs_time.png'))

#%% 
'''
# plot the live version of species number/ coverages plot over time
'''
if Count_to_Coverage == 0:
    ylabel = 'Species count'
    s_df_live = s_df
else:
    ylabel = 'Coverage'
    # Convert to surface coverages    
    s_df_live = pspec.num_to_cov_df(lattice_dim, s_df)

dspec.PlotLiveSpec(s_df_live, xlab = 'Time (s)', ylab = ylabel, xlimit = xrange,
                   fname =  output + '/' + 'live_surf_vs_time.html')



#%%
'''
# plot the distribution of surface species at the end state
'''
# end state species info
es_n = pspec.final_lt_spec_on_surf(s_df)

spec_cov = np.multiply(np.array(surf_spec_n),np.array(es_n))

dspec.PlotPie(spec_cov, spec_name,  fname = os.path.join(output_dir, 'surf_spec_pie.png'))

#%%
'''
# Plot a bar graph of elementary step frequencies versus time - output in elem_step_freqs.png in the directory with the Zacros run
'''
# frequency object
f = pspec.read_Multiple_Procstat()
freq_vecs = f.ave_procstat_freqs
dspec.PlotFreqs(freq_vecs, fname = os.path.join(output_dir, 'elem_step_freqs.png'))
