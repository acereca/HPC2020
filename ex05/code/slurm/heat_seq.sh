#! /bin/bash

#SBATCH -o out/%A

for s in 128 512 1024 2048 4096
do
    ./heatgrid_seq.out -s $s -i 100
done
