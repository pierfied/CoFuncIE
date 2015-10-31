#!/bin/bash
#BSUB -n 32
#BSUB -R "span[ptile=16]"
#BSUB -o lsf.out
#BSUB -e lsf.err
#BSUB -q "windfall"
#BSUB -J num_core_test

module load openmpi
mpirun -cpus-per-proc 16 -pernode ./num_core_test
