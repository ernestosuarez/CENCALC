#!/bin/bash
export OMP_NUM_THREADS=2  
export KMP_STACKSIZE=100m

cencalc_omp  -ns 1000 5000 1000 -o 1 -t mie_o1_evol.tab > mie_o1_evol.out

