# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 14:03:30 2019

@author: yifan
"""

from mpi4py import MPI
import math
comm = MPI.COMM_WORLD

print("%d of %d" % (comm.Get_rank(), comm.Get_size()))

import numpy as np 

def rbind(comm, x):
    return np.vstack(comm.allgather(x))

def rbind2(comm, x):
    
    size =  comm.Get_size()
    m = np.zeros((size, len(x)), dtype = np.int)
    comm.Allgather([x, MPI.INT], [m, MPI.INT])
    
    return m
    
comm = MPI.COMM_WORLD
x = np.arange(4, dtype = np.int) * comm.Get_rank()
a = rbind(comm, x)
a = rbind2(comm, x)

x = range(20)
m = int(math.ceil(float(len(x))/size))
x_chunk = x[rank*m : (rank+1)*m]
r_chunk = map(math.sqrt, x_chunk)
r = comm.allreduce(r_chunk)