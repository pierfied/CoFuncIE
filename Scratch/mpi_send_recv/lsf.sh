#!/bin/bash
#BSUB -n 32
#BSUB -R "span[ptile=16]"
#BSUB -o lsf.out
#BSUB -e lsf.err
#BSUB -q "windfall"
#BSUB -J send_recv_test

module load openmpi
mpirun -np 32 ./send_recv
