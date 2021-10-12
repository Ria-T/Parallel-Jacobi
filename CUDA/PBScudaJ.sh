#!/bin/bash

# Which Queue to use, DO NOT CHANGE #
#PBS -q GPUq

# Max Wall time, Example 1 Minute #
#PBS -l walltime=00:20:00

# How many nodes and tasks per node,  1 nodes with 8 tasks(threads)#
#PBS -lselect=1:ncpus=1:ompthreads=1:ngpus=1 -lplace=excl

# JobName #
#PBS -N JacobiCUDA

#Change Working directory to SUBMIT directory
cd $PBS_O_WORKDIR

# Run executable #
./jacobi_cuda.x < input

# profile executable #
#nvprof ./jacobi_cuda.x < input

