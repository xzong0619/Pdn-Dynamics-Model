#!/bin/bash
#$ -cwd
#$ -j y
#$ -N Single_trajectory
#$ -S /bin/bash
#$ -pe mpich 1
#$ -o Single_trajectory.out

# The  executable:
export KMC_EXE="/home/vlachos/wangyf/programs/zacros_2.1/build/zacros.x"

# Simple summary:
echo ""
echo "Running on ${HOSTNAME} with job id ${JOB_ID}"
echo ""

${KMC_EXE}

# Get our environment setup:
vpkg_require "python/2.7.8"
vpkg_require "openmpi/1.6.3-gcc"
vpkg_require "python-numpy"
vpkg_require "python-scipy"
#vpkg_require "python-mpi4py"
#vpkg_rollback all
#vpkg_require python-networkx

# The  executable:
export PYTHON_EXE="python Single_trajectory.py"

# Simple summary:
echo ""
echo "Running on ${HOSTNAME} with job id ${JOB_ID}"
echo ""

time ${PYTHON_EXE}