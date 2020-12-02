#! /bin/bash

#SBATCH --ntasks-per-node=4 --ntasks=16
#SBATCH -o out/%A

module load mpi
for proc in 10 16
do
    for size in {7..13}
    do
        mpirun -np ${proc} -npernode 4 ./mmul_para.out -n $((2 ** ${prob}))
    done
done
