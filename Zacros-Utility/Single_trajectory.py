# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import os
import sys

sys.path.append('/home/vlachos/wangyf/programs/ZW_git/Zacros-Wrapper')
import zacros_wrapper as zw

''' ------------ User input section ------------ '''
RunPath = os.getcwd()
''' -------------------------------------------- '''

''' Set up data '''
y = zw.kmc_traj()
y.Path = RunPath
y.ReadAllOutput(build_lattice = True)

''' Analyze '''
n_Pd = 1                                                     # number of Pd atoms for normalization
y.PlotSurfSpecVsTime(site_norm = n_Pd)                              # Plots species coverages
y.PlotGasSpecVsTime()                                               # Plots total counts of gas species
y.PlotElemStepFreqs(time_norm = False)             # Plots elementary step frequencies
y.PlotLattice()                                                     # Plots a picture of the bare lattice
y.LatticeMovie()                                                    # creates a subfolder and makes images of each snapshot