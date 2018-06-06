# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import sys
import os
sys.path.append('/home/vlachos/wangyf/Zacros-Wrapper')
import zacros_wrapper as zw

################## User input ##################################

zacros_exe = '/home/vlachos/wangyf/programs/zacros_2.1/build/zacros.x'
KMC_source = os.getcwd()
output = 'outputs'
BatchPath = os.path.join(KMC_source, output)
print(KMC_source)
print(BatchPath)
n_runs = 10

################################################################

# Prepare template
x = zw.Replicates()

x.runtemplate = zw.kmc_traj()
x.runtemplate.Path = KMC_source
x.runtemplate.exe_file = zacros_exe
x.runtemplate.ReadAllInput()
x.ParentFolder = BatchPath

# Build files and run
x.n_runs = n_runs
x.BuildJobFiles()
x.RunAllJobs_parallel_JobArray()

# Read results
x.ReadMultipleRuns()
x.AverageRuns()

# Plot results
x.runAvg.PlotGasSpecVsTime()
x.runAvg.PlotSurfSpecVsTime()
x.runAvg.PlotElemStepFreqs()
