#!/bin/bash

# Number of processes - Adapt with 'param' file
declare -i nproc=64
# Number of output files
declare -i n=32
# Index of output file
declare -i i

# Main loop for batch execution
for ((i=1;i<=$n;i=2*i))
do
   # Execution of parallel version
   mpirun -np $i ./explicitPar_reduced  $i : -np $nproc ./explicitPar  $i
   # Save output data
   mv outputPar.dat outputPar_regen$i.dat
done
