#! /bin/bash

#SBATCH --ntasks-per-node=1 --ntasks=2
#SBATCH -o out/%A_%a_1n
#SBATCH --array 0-10

module load mpi
mpirun -host creek01,creek05 -np 2 ./latency.out -n 100000 -s $((2 ** ${SLURM_ARRAY_TASK_ID})) -N 2
