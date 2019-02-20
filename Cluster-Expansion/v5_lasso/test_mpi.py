# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 14:03:30 2019

@author: yifan
"""

from mpi4py import MPI
comm = MPI.COMM_WORLD
print("%d of %d" % (comm.Get_rank(), comm.Get_size()))