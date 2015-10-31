#!/bin/bash
#BSUB -n 8
#BSUB -R "span[ptile=2]"
#BSUB -o lsf.out
#BSUB -e lsf.err
#BSUB -q "windfall"
#BSUB -J hello_world

module load openmpi
mpirun -np 4 ./mpi_hello_world > output.txt
