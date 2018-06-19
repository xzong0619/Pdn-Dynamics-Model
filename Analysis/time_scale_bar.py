# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 13:32:30 2018

@author: wangyf
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib


# working on all platforms, using cross-platform home path

HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Generator'))
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Generator','user_inputs'))

from nc import *
from IO_Pdn import *




#%% Buliding simulation modules
T = 700 # Temperature (K) 

P = 1 # Pressure (atm)
n_CO = 0 # 1 if there is CO
n_hours = 1/60 # simulation max time, number of hours


flg_on_events = 0 # if sampling on events, sampling on time by default
flg_debug = 1 # if debugging mode is on
flg_no_restart = 1  # if restart
single_rxn_mode = 0 # if creat simulation for each rxn

flag = [flg_on_events, flg_debug, flg_no_restart]

# Initial state
d1 = 'Pd1*' # specify the initial surface species, add species name with a *
d2 = 40 # specify the corresponding number on the surface 
di = 1 # number of types of initial surface species

Base_path = os.getcwd()

Output_fldr = 'outputs'
Parent_dir = os.getcwd()
output_dir = os.path.join(Parent_dir, Output_fldr)

Sim = SimOut(n_CO,  T,  P,  n_hours,  flag)
Cluster = ClusterOut(nc)
Mech = MechanismOut(nc)
State =  StateOut(nc)
fldr = os.path.join(Base_path, Output_fldr)

if not os.path.exists(fldr):
    os.makedirs(fldr)
    
Cal = RxnCal()     
Cal.Tdependence(T)

#%% Plot the bar graph
fig = plt.figure()
fs = 25
x = np.arange(54)
s_factor = 10**10

T_fwd = Cal.Timescale_fwd
T_rev = Cal.Timescale_rev

for i in np.arange(-5,0):
    T_fwd[i] = T_fwd[i]/s_factor
    T_rev[i] = T_rev[i]/s_factor
    
yfwd = np.log10(T_fwd)
baseline = -13
yrev = np.log10(T_rev)

width = 0.4
fig, ax = plt.subplots()
fig.set_figwidth(30)
fig.set_figheight(10)
rects1 = ax.bar(x, yfwd - baseline, width, bottom = baseline,  color='b')

for i in x:
        
    if i in np.arange(20):
        rects1[i].set_color('palegreen')
    if i in np.arange(20, 30):
        rects1[i].set_color('lime')
    if i in np.arange(30, 34):
        rects1[i].set_color('yellowgreen')
    if i in np.arange(34, 49):
        rects1[i].set_color('olivedrab')
    if i in np.arange(49, 54):
        rects1[i].set_color('seagreen')
            
rects2 = ax.bar(x + width, yrev - baseline,  width, bottom = baseline, color='g')

for i in x:
        
    if i in np.arange(20):
        rects2[i].set_color('steelblue')
    if i in np.arange(20, 30):
        rects2[i].set_color('c')
    if i in np.arange(30, 34):
        rects2[i].set_color('deepskyblue')
    if i in np.arange(34, 49):
        rects2[i].set_color('lightsteelblue')
    if i in np.arange(49, 54):
        rects2[i].set_color('lightblue')
        
# add some text for labels, title and axes ticks
ax.set_ylabel(r'$log_{10}$'+'(Time Scale (s))')
ax.set_title('T = ' + str(T) + 'K')
ax.set_xlabel('Events Index')
ax.set_ylim([-13, 28])
ax.tick_params(labelsize=fs)
matplotlib.rcParams.update({'font.size': fs})

plt.show()
fig.savefig(fname = os.path.join(output_dir, str(T) + 'K_time_scale.png'))