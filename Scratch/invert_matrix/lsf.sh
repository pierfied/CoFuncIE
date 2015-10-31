###========================================
#!/bin/bash
#BSUB -n 1
#BSUB -o lsf.out
#BSUB -e lsf.err
#BSUB -q "windfall"
#BSUB -J invert_test
#---------------------------------------------------------------------


time ./invert_test
