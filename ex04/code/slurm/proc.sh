#! /bin/bash

#SBATCH --ntasks-per-node=4 --ntasks=16
#SBATCH -o out/%A

module load mpi
mpirun -np 1 -npernode 4 ./mmul_para.out -n 2000
for proc in {2..16..2}
    do
        mpirun -np ${proc} -npernode 4 ./mmul_para.out -n 2000
done
