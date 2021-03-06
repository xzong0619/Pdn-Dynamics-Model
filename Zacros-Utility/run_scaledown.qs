#!/bin/bash
#$ -cwd
#$ -j y
#$ -N scaledown
#$ -S /bin/bash
#$ -pe openmpi-smp 1
#$ -o scaledown.out

# Get our environment setup:
vpkg_require "python/2.7.8"
vpkg_require "openmpi/1.6.3-gcc"
vpkg_require "python-numpy"
vpkg_require "python-scipy"
#vpkg_require "python-mpi4py"
#vpkg_rollback all
#vpkg_require python-networkx

# The  executable:
export PYTHON_EXE="python scaledown.py"

# Simple summary:
echo ""
echo "Running on ${HOSTNAME} with job id ${JOB_ID}"
echo ""

time ${PYTHON_EXE}