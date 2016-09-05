#!/bin/bash
export OMP_NUM_THREADS=2        
export KMP_STACKSIZE=100m

for cut in  0 3 4 5 7 8 9
do
cencalc_omp  -cutoff $cut -table ccgla_cut${cut}.tab  > ccgla_cut${cut}.out
done

