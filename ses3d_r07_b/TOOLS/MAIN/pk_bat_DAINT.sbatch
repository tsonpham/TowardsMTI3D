#!/bin/bash -l
#SBATCH --account=ch1
#SBATCH --ntasks=8
#SBATCH --time=01:00:00

aprun -n 8 ./project_kernel_batch ../../DATA/OUTPUT/1.8s ../../DATA/OUTPUT/1.8s 0 1 4

