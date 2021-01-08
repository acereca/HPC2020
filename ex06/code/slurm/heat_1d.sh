#! /bin/bash

#SBATCH --ntasks-per-node=4 --ntasks=12
#SBATCH -o out/%A

module load mpi

mkdir -p out
echo "#size;nodes;iterations;t/us\n" > out/grid.csv  # clear file

for nodes in 1 2 4 6 8 10 12
do
  for size in 128 512 1024 2048 4096
  do
    mpirun -np ${nodes} ./heatgrid_1d.out -s ${size} -i 100
  done
done
