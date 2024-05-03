#!/bin/bash -l
#SBATCH --partition=fichtner_compute
#SBATCH --ntasks=8
#SBATCH --time=01:00:00

mpirun -n 8 ./project_kernel_batch ../../DATA/OUTPUT/1.8s ../../DATA/OUTPUT/1.8s 0 1 4

