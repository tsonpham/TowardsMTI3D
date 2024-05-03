#!/bin/bash -l
#SBATCH --partition=fichtner_compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00

srun -c 1 -n 1 --ntasks-per-node 1 ./project_kernel ../../DATA/OUTPUT/1.8s ../../DATA/OUTPUT/1.8s 0

