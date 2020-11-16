#! /bin/bash

#SBATCH --ntasks-per-node=6 --ntasks=24
#SBATCH -o out/%A_%a
#SBATCH --array 2-24:2

module load mpi
mpirun -host creek01,creek06,creek05,creek04 -np ${SLURM_ARRAY_TASK_ID} ./ring.out -n 10000000
