#!/bin/bash

# JobName #
#PBS -N Jacobi

#Which Queue to use #
#PBS -q N10C80

# Max Wall time, Example 1 Minute #
#PBS -l walltime=00:20:00

# How many nodes and tasks per node
#PBS -l select=5:ncpus=5:mpiprocs=5:mem=16400000kb

#Change Working directory to SUBMIT directory
cd $PBS_O_WORKDIR

# Run executable #
mpirun jacobi_mpi.x < input2
