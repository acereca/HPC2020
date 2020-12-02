#! /bin/bash

#SBATCH --ntasks-per-node=1 --ntasks=2
#SBATCH -o out/%A_%a_2n
#SBATCH --array 0-10

module load mpi
mpirun -np 2 -npernode 1 ./latency.out -n 100000 -s $((2 ** ${SLURM_ARRAY_TASK_ID})) -N 2
