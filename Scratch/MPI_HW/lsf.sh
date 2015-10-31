#!/bin/bash
#BSUB -n 128
#BSUB -R "span[ptile=16]"
#BSUB -o lsf.out
#BSUB -e lsf.err
#BSUB -q "windfall"
#BSUB -J cart_catalog

module load openmpi

# Get time as a UNIX timestamp (seconds elapsed since Jan 1, 1970 0:00 UTC)
T="$(date +%s)"

mpirun -n 128 ./mpi_convert_cart full_data.txt

T="$(($(date +%s)-T))"
echo "Runtime in seconds: ${T}"
