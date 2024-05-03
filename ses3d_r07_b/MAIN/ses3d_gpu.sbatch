#!/bin/bash -l

#SBATCH --account=ch1
#SBATCH --nodes=6
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00

aprun -d 1 -n 6 -N 1 ./ses3d_gpu

