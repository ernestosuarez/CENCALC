#!/bin/bash
export OMP_NUM_THREADS=2        #Number of cores
export KMP_STACKSIZE=100m

cencalc_omp  -ns 1000 5000 1000 -cutoff 7 -table ccmla_cut7_evol.tab  > ccmla_cut7_evol.out

