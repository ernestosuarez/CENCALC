INSTALLATION
============

The CENCALC software consists mainly of two standalone codes written in the 
FORTRAN90 language, cencalc_prep.f90 and cencalc_omp.f90. The first program 
carries out various preparatory tasks prior to the main entropy calculations 
that are performed by cencalc_omp.f90. As the bulk of the effort to obtain 
the conformational entropy values is expended by cencalc_omp.f90, this 
program takes advantage of shared-memory parallel computers 
through the OpenMP Application Program Interface. 

In order to "install" the program, the first step would be to compile the 
FORTRAN90 code using Intel FORTRAN Compilers (recommended):

>>  ifort  cencalc_prep.f90 -o cencalc_prep
>>  ifort  -openmp cencalc_omp.f90 -o cencalc_omp

In addition CENCALC also includes a Python program, get_tor.py, that reads 
the topology information from AMBER parm files and identifies all the torsions
about single bonds that are required to characterize the conformational state 
of the molecule of interest. 

The second step would be to move the cencalc_prep and cencalc binaries 
together with the get_tor.py script to a run directory included in the PATH
environment variable. 

