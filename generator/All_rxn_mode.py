# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 16:39:53 2018

@author: wangyf
"""
from IO_Pdn import *
'''
User Input section (reaction conditions):
--------------------------------------------------------------
'''
HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Pdn-Dynamics-Model', 'Generator','user_inputs'))
from nc import *

T = 300 # Temperature (K) 

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

Output_fldr = 'Output_All'

Sim = SimOut(n_CO,  T,  P,  n_hours,  flag)
Cluster = ClusterOut(nc)
Mech = MechanismOut(nc)
State =  StateOut(nc)
fldr = os.path.join(Base_path, Output_fldr)

if not os.path.exists(fldr):
    os.makedirs(fldr)
            
Sim.WriteIn(fldr)
        
Cluster.WriteIn(fldr)
        
Mech.WriteIn(fldr)
       
State.WriteIn(fldr, d1,  d2,  di)

Cal = RxnCal()     
Cal.Tdependence(T)