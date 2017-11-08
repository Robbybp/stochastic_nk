#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=scaling
#SBATCH --time=72:00:00


srun julia main.jl -f pglib_opf_case14_ieee_nk.m --algo Lshaped --batchsize 100 --batchid 1 --numbatches 1 -k 2 --solver cplex

