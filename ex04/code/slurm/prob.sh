#! /bin/bash

#SBATCH --ntasks-per-node=4 --ntasks=16
#SBATCH -o out/%A

module load mpi
for exp in {1..10}
do
    for proc in 10 16
    do
        for prob in {7..12}
        do
            mpirun -np ${proc} -npernode 4 ./mmul_para.out -n $((2 ** ${prob}))
        done
    done
done
