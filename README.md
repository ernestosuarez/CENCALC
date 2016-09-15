# Conformational Entropy Calculation

## _CENCALC Users’ Manual_

##Preface 
This software has been designed to estimate the conformational entropy of single molecules from extended computer  simulations,  especially  Molecular  Dynamics  (MD) simulations. On  input CENCALC needs both  trajectory  coordinates  and  topology  information  in  order  to  characterize the conformational states of the molecule of interest. The molecular conformers are identified by discretizing  the  time  evolution  of  internal  rotations.  After  this  transformation,  CENCALC determines the  probability  mass  functions  of  the  individual  torsions  and  uses  them for conformational  entropy  estimations.  CENCALC  can  use up  to four  different  methodologies for approaching  to  the  full  conformational  entropy: the  classical  Mutual  Information  Expansion (MIE),  the Approximate  MIE  (AMIE),  the  so-called  Multibody  Local  Approximation  (MLA), and  the  default  method that corresponds  to  the  correlation  corrected  MLA  (CC-MLA).  All  of these  techniques  can  also  be  combined  with  a  distance-based  cutoff  criterion. In this case, CENCALC  requires as additional input an inter-atomic  distance  matrix  containing  the  mean distance values derived from the MD trajectory in order to include only correlation effects among torsion  angles  whose  mean  separation  is  below a predefined  cutoff.  The  best  cutoff  for  a given amount of sampling can be determined using the CC-MLA method.

All the assumptions and equations defining the various techniques available in CENCALC have been  discussed  in  the  literature as  well  as in  the  accompanying  paper. Users  are  therefore encouraged to consult those references cited below before using the software.  

##0-License and Citation 
Copyright (C) 2012 Ernesto Suárez Álvarez CENCALC is free software: you can redistribute it and/or modify it under the terms  of  the  GNU  General  Public License  as  published  by  the  Free  Software Foundation, either version 3 of the License, or any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>. Any use of the CENCALC software or derivative should include at least the following citation: 1) E.  Suárez,  N.  Díaz,  J.  Méndez  and  D.  Suárez. CENCALC:  A  Computational  Tool  for Conformational Entropy Calculations from Molecular Simulations.  Submitted. The methods implemented in CENCALC are fully described in the following references:  2) E.  Suárez,  N.  Díaz  and  D.  Suárez. Entropy  Calculations  of  Single  Molecules  by Combining  the  RigidRotor  and  Harmonic-Oscillator  Approximations  with Conformational  Entropy  Estimations  from  Molecular  Dynamics Simulations J.  Chem. Theor. Comput. 2011, 7, 2638-2653. 3) E.  Suárez,  D.  Suárez. Multibody  Local  Approximation:  Application to Conformational Entropy Calculations on Biomolecules. J. Chem. Phys. 2012, 137, 084115. All  questions  regarding  the  usage  and distribution  of  CENCALC  or  bug  reports  should  be addressed to Ernesto Suárez (esuarez@pitt.edu  ; ernesto.suarez.a@gmail.com ).  The  software,  which  is  distributed  under  the  GNU  public  license,  together  with  numerical examples  and  this  user’s  manual  are  available  in  the  Supporting  Information  of  reference  1.

##1-Compliling the code 
The  CENCALC  software  consists  mainly  of  two standalone codes  written  in  the  FORTRAN90 language, cencalc_prep.f90 and cencalc_omp.f90.  The  first  program  carries  out  various preparatory tasks prior to the main entropy calculations that are performed by cencalc_omp.f90. As  the  bulk  of  the  effort  to  obtain  the  conformational  entropy  values  is  expended  by cencalc_omp.f90, this  program  takes advantage  of  shared-memory  parallel  computers  through the OpenMP Application Program Interface.

Important  notice: The FORTRAN90 codes in CENCALC have been compiled and fully tested using the Linux version of the gfortran compiler (v. 4.4.5):

```bash
>>  gfortran  cencalc_prep_0.2.2.f90 -o cencalc_prep 
>>  gfortran  -fopenmp cencalc_omp_0.2.2.f90 -o cencalc_omp 
```

##2-Running the code 

##2.1-Files needed to perform conformational entropy calculations  
The following Table describes the data files, in plain text format, that are needed by CENCALC in order to carry out conformational entropy calculations.

| Content       | Filenames*     | Observations  |
| ------------- |:-------------:| -------------:|
| Time series of torsion angles     | d0001.dat, d0002.dat ...  |  Mandatory          |
| Distance matrix      | distance_matrix.dat (Format: F9.3)      |   Only for calculations with cutoff        |
| ID numbers of the central atoms involved in the M torsions  | atoms_in_tor.info | Only for calculations with cutoff |
(*) Default filenames are indicated, but users can choose other filenames using command line options (see below).

In  principle  the d0001.dat, d0002.dat,  ... , d000M.dat files have  each two  columns corresponding  to  the  time  counter  and  the  torsion  angle  value,  respectively. In  any  case CENCALC reads  only  one  column  per  file  whose  column  index  can  be  specified  using  a command line option (see below). Torsion angles must be specified as dihedral angles measured in degrees. Of course the d0001.dat, d0002.dat, ... files should all have the same number N of lines corresponding to the number of MD snapshots being considered. It may be noted that  a time  space  of 1ps between  the  MD  snapshots seems  appropriate  for  most  applications of CENCALC. On the other hand, the value of N will depend on the dimensionality of the problem as  well  as  on  the  conformational  flexibility  of  the  molecular  system: for  medium-sized  systems one should expect that millions of configurations would be required to reach converged entropies.

The file atoms_in_tor.info must have M lines, with M being the number of torsion data files. If, for example, d0001.dat and d0002.dat store the time series of the torsion angles defined by the 1-2-3-4 and 2-3-4-5 atom ID numbers, then the first two lines in atoms_in_tor.info should be  simply 2  3 and 3  4,  respectively.  On  the  other  hand,  the file distance_matrix.dat should contain, in format F9.3, the inter-atomic mean distance matrix of the solute, or at least, the inter-atomic mean distance matrix of the first K atoms in the topology, where K is greater or equal than the highest atom ID number in atoms_in_tor.info. Note that coordinates of solvent molecules and  counterions  should  be  better  removed  before  generating distance_matrix.dat, otherwise the resulting distance matrix could be huge and very expensive to compute.

Of  course  CENCALC  assumes  that  a mutually consistent  atom numbering  is  used  in  the construction of the *atoms_in_tor.info* and *distance_matrix.dat* files. 

##2.2-Runing cencalc_prep 
The cencalc_prep code estimates first the probability density functions of the M torsion angles and characterizes their maxima and minima critical points. Basing on these data, cencalc_prep transforms subsequently the initial time series of N real numbers per torsion angle $\theta$ into a set of N integer  numbers labeling  the  conformational  states  populated by $\theta$.  On  output,  all  the information  is  saved in  a  file named MATRIX.dat,  which has N rows (i.e,  the  number  of  MD snapshots) and M columns (the number of rotatable bonds). Thus, the ith row of MATRIX.dat is an  array  of integer  numbers  that  represent the  conformational  state at the ith-snapshot. In principle, MATRIX.dat should have as  many  columns  as  times  series read  on  input (d0001.dat, d0002, ...). However, cencalc_prep removes by default all the frozen torsions because they do not  represent  conformational  changes  and  do  not  affect  the  conformational  entropy  and, consequently, M could be lower on output than on input. In addition cencalc_prep reads the atoms_in_tor.info and distance_matrix.dat files and generates  a  new  distance  matrix file  named reduced_dist_matrix.dat that  contains the mean distance between  every  pair  of torsions. These distance  values  are  derived  by  applying the following rules: a) The distance d(A,A) between every torsion A and itself is zero; b) the distance d(A,B) between two different torsions A1-A2-A3-A4 and B1-B2-B3-B4 is the mean distance d(Ai,Bj) between the central pair of atoms of both torsions.

A typical execution of cencalc_prep assuming default filenames looks like:

```bash
>>  cencalc_prep  d????.dat 

Any other regular expression besides d????.dat can be used. 
The program will also look for atoms_in_tor.info and distance_matrix.dat. 
```

or without cutoff:

```bash
>> cencalc_prep -nocut d????.dat

Calculations with no cutoff. The program will not look for the atoms_in_tor.info 
and distance_matrix.dat files.
```
For the more general case, we provide the full help of cencalc_prep:

___
**Usage of cencalc_prep** 

**SYNOPSIS:** 

cencalc_prep [OPTIONS]  file1.dat file2.dat ... 
(Default names should be d0001.dat, d0002.dat ...) 

**OPTIONS:** 

**-u/-usecol**  COLNUM    Default: 2 
    
    Column number in file1.dat, file2.dat that contains the times series of the torsion variable. This value is normally 2 since first column often corresponds to the time  or the snapshot number variable. 

-dist  DIST_MATRIX_FILE_NAME  Default: distance_matrix.dat 
    
    This variable specifies the filename of the inter-atomic mean distance matrix. 
    
-i/-info TOR_INFO_FILE_NAME  Default: atoms_in_tor.info 
    
    This file specifies the atoms involved in the torsion angles (only the two central atoms).  For example, if the first row in *TOR_INFO_FILE_NAME*  reads as 3 4,  then file1.dat contains the time series for rotation about the 3-4 bond.
    
 nocut        Default: Use cutoff 
Using this option no cutoff will be applied and the options  
dist/info  are thereby not needed.  s/simplify  yes/no       Default: yes Remove all frozen torsions.  k K_VALUE     Default: 1.0 The k_value (ˆ) sets the smoothing parameter vin the vonMises kernel density estimation as proposed in Eq.(7) in Computational Statistics & Data Analysis, 2008, 52, 34933500. By default ˆ1 as this value ensures slightly over-smoothed Probability Density Functions (PDFs) of individual torsion angles, what is convenient for searching critical points.  ag  yes/no      Default: no If yes then analytic gradients are used for locating the minima critical points of the von-Mises PDFs of individual torsion angles. The default option (no) uses instead a sufficiently accurate and fast linear interpolation scheme for PDF gradient evaluation.  step  STEP_SIZE    Default: 5 (degrees) PDF minimization step size:  xn+1=xn STEP_SIZE * GRADIENT  crit  CONVERGENCE_THRESHOLD  Default: 1.0·10 -4Gradient convergence threshold for the PDF minimizations. 
___
