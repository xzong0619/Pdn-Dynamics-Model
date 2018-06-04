# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 16:26:30 2018

@author: wangyf
"""

from IO_Pdn import *
'''
User Input section (reaction conditions):
--------------------------------------------------------------
'''
nc = [None,  # neighboring structure
      '1-2',
      '1-2 1-3 2-3',
      '1-2 1-3 2-3 4-2 4-3 ',
      '1-2 1-3 2-3',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 ',
      '1-2 1-3 2-3 4-2 4-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-4 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-4 6-2 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 ',
      '1-2 1-3 2-3 4-2 4-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 7-6 7-3 1-7 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-4 7-6 7-5 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 7-6 7-5 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-4 6-2 7-6 7-2 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-4 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-4 6-2 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 ']


T = 700 # Temperature (K) 

P = 1 # Pressure (atm)
n_CO = 0 # 1 if there is CO
n_hours = 1/60 # simulation max time, number of hours


flg_on_events = 1 # if sampling on events, sampling on time by default
flg_debug = 1 # if debugging mode is on
flg_no_restart = 1  # if restart
single_rxn_mode = 1 # if creat simulation for each rxn

flag = [flg_on_events, flg_debug, flg_no_restart]



'''
--------------------------------------------------------------
'''
rxn = ReadRxn('rxn_input.txt')
n_rxn = len(rxn)
Base_path = os.getcwd()
'''
Output_fldr = 'Output'

Sim = SimOut(n_CO,  T,  P,  n_hours,  flag)
Cluster = ClusterOut(nc)
Mech = MechanismOut(nc)
State = StateOut(nc)

if single_rxn_mode:
    for r in range(n_rxn):
       
        iter_fldr = os.path.join(Base_path, Output_fldr,  str(r+1))
        if not os.path.exists(iter_fldr):
            os.makedirs(iter_fldr)
        
        Sim.WriteIn(iter_fldr)
        
        Cluster.WriteIn(iter_fldr)
        
        Mech.WriteIn(iter_fldr, r)
       
        State.WriteIn_Single_rxn(iter_fldr, r)
        
        #MakingCopy('zacros.qs', iter_fldr)
        
        MakingCopy('lattice_input.dat', iter_fldr)
        
        MakingCopy('Single_trajectory.py', iter_fldr)
        
        MakingCopy('run_single.qs', iter_fldr)
'''
Cal = RxnCal()     
Cal.Tdependence(T)