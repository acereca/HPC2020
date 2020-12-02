#! /bin/bash

#SBATCH --ntasks-per-node=2 --ntasks=2
#SBATCH -o out/%A_%a_1n
#SBATCH --array 0-10

module load mpi
mpirun -np 2 -npernode 2 ./bandwidth.out -n 1000 -s $((2 ** ${SLURM_ARRAY_TASK_ID})) -N 1 -r 100 -b 0
mpirun -np 2 -npernode 2 ./bandwidth.out -n 1000 -s $((2 ** ${SLURM_ARRAY_TASK_ID})) -N 1 -r 100 -b 1
