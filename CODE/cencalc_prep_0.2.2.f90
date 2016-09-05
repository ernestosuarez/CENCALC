!      
! Copyright (C) 2011 Ernesto Suarez Alvarez
!                    ernesto@fluor.quimica.uniovi.es
!                    ernesto.suarez.a@gmail.com
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  Any use of the CENCALC software or derivative should include at least the following 
!  citation:
!
!  1)E. Suarez, N. Diaz, J. Mendez and D. Suarez. CENCALC: A Computational Tool for 
!    Conformational Entropy Calculations from Molecular Dynamics Simulations. 
!    J. Chem. Inf. Model. 2011, 
!
!  The methods implemented in CENCALC are fully described in the following references: 
!
!  2)E. Suarez, N. Diaz and D. Suarez. Entropy Calculations of Single Molecules by 
!    Combining the Rigid-Rotor and Harmonic-Oscillator Approximations with Conformational 
!    Entropy Estimations from Molecular Dynamics Simulations 
!    J. Chem. Theor. Comput. 2011 (Accepted).
!
!  3)E. Suarez, D. Suarez. Multibody Local Approximation: Application in the Conformational 
!    Entropy Calculation on Biomolecules.  2011 (submitted).
!
!  All questions regarding the usage and distribution of CENCALC or bug reports should be
!  addressed to Ernesto Suarez (ernesto@fluor.quimica.uniovi.es). 
!-----------------------------------------------------------------------------------------



!*****************************************************************************************
      MODULE parameters
!*****************************************************************************************
      implicit none
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(10)
      INTEGER, PARAMETER :: SMALL = SELECTED_INT_KIND(1)
      real(DP), parameter :: Pi=3.14159265359
      END module parameters
!-----------------------------------------------------------------------------------------



!     _____________________________
!    /|                           |
!    /|       MAIN PROGRAM        |
!    /|___________________________|
!    //////////////////////////////


!*****************************************************************************************
      PROGRAM CENCALC_PREP
!*****************************************************************************************
! cencalc_prep 
!
! This program prepares "MATRIX.dat" and "reduced_dist_matrix.dat" 
! that are required to run centro_omp program.
!
! Files needed:
!
!        1 - All the time series of each dihedral angle (one file by dihedral,
!            for example: file1.dat file2.dat ...). In each "file*.dat" the first 
!            column is normally the time in (ps) that will not be read, and the 
!            second which is the important one, is the dihedral angle value.
!            It also possible to read another column using the option -usecol. 
!            For more info use the option -h/-help.
!
!             
!        The Cut-off option is activated by default, in this case it 
!        is also needed:
!
!        2 - The distance matrix of the whole protein or system in vacuo 
!            using the file distance_matrix.dat
!        3 - The information about which atoms are involved in each dihedral.
!            The file atoms_in_dih.info should specify the two central atoms
!            of each angle. The first line corresponds to the first dihedral whose
!            time series is in file1.dat and so on.
!
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! QUICK HELP:
!
! A quick help can be viewed from the command line by doing: centro_prep -help 
!-----------------------------------------------------------------------------------------

      use parameters 
      IMPLICIT NONE
      !-VARIABLE DEFINITIONS--------------------------------------------------------------
      integer NumFiles                                        !Number of data files 
      integer UseCol                                          !Which column use in the data files
      integer NumAtoms                                        !Number of atoms in the distance matrix (lines or columns)
      integer NumSnap                                         !Number of snapshots
      character*60  DistMatrixName                            !Distance matrix file name 
      character*60  InfoFileName                              !Info file name
      character*60,dimension (:),allocatable :: FileList      !List of data files (time evolution of torsions)
      integer NumArg                                          !Number of arguments in the command line
      character*60,dimension (:),allocatable :: Arguments     !Arguments in the command line
      character*3 Simplify                                    !Eliminate conformationally rigid torsions? (True/False)
      logical UseCut                                          !Use cutoff? (True/False)
      integer(SMALL),dimension (:,:),allocatable :: BigMatrix !Discretized data matrix before removing rigid torsions
      integer(SMALL),dimension (:,:),allocatable :: Matrix    !Discretized data matrix after removing rigid torsions
      integer,dimension (:,:),allocatable :: Atoms_in_Dih     !Two-column matrix where is saved the content of InfoFileName
      real(DP),dimension (:,:),allocatable :: Dist            !Distance matrix
      real(DP),dimension (:,:),allocatable :: NewDistMat      !Output: reduced distance matrix
      integer NewNumCol                                       !New number of columns (the columns in reduced dist. matrix)
      real(DP),dimension (:),allocatable :: ang               !Dihedral angle value in degrees
      integer(SMALL),dimension (:),allocatable :: DiscreteAng !Discretized dihedral angle
      real(DP) k_value                                        !(see explanation in default values)
      logical  AnalyticGrad                                   !Use analytic gradient in the optimization?
      real(DP) Step                                           !Step for the steepest-descendent optimization
      integer  MaxIterations                                  !Max number of iteration in the optimization
      real(DP) ConvCriterion                                  !Convergence criterion
      integer  MaxNumConf                                     !Max Number of Conformers to look for
      integer  FNumMins                                       !final num of minimums
      real(DP),dimension (9) :: minimum                       !Positions of the minimums in the PDF    

      logical,dimension (:),allocatable :: selected           !
      real(DP) dummy                                          !
      character*1 charact !                                    > Auxiliary or Dummy variables
      character*60 cdummy                                     !
      integer idummy,idummy2,iargc,ios                        !
      integer i,j,k,C,A                                       !
      !-----------------------------------------------------------------------------------


      write(*,'(A33)') " _______________________________ "
      write(*,'(A33)') "||  _____                       |"
      write(*,'(A33)') "||  \\   ||                     |"
      write(*,'(A33)') "||   \\  PROGRAM CENCALC_PREP   |"
      write(*,'(A33)') "||   //          v0.2           |"
      write(*,'(A33)') "||  //___||                     |"
      write(*,'(A33)') "||______________________________|"
      print*," "


!---DEFAULT VALUES--------------------------
      UseCol=2                                                !Which column is going to be discretized?
                                                              
      Simplify="yes"                                          !Eliminate columns with null entropy 
                                                              
      InfoFileName="atoms_in_dih.info"                        !Which atoms define de dihedral
                                                              
      DistMatrixName="distance_matrix.dat"                    !All solute atoms distance matrix
                                                              
      UseCut=.true.                                           !Use cut-off
                                                              
      k_value=1.0                                             !The smoothing parameter "v" of
                                                              !the von-Mises kernel density estimation 
                                                              !depends of the k_value see eq.(7) in ref:
                                                              !Computational Statistics & Data Analysis
                                                              !Volume 52, Issue 7, 15 March 2008, Pages 3493-3500.
                                                              !Here we do not estimate k_value, but simply
                                                              !set the value empirically in order to
                                                              !slightly over-smooth the PDFs. In any case, the user
                                                              !can change this value by the option "-k".
                                                              
      MaxIterations=1000                                      !Max number of Iterations in order to
                                                              !find the PDFs minimums
      Step=5                                                  !Step-size during the optimization
                                                              
      ConvCriterion=1E-4                                      !Convergence criterion (gradient)
                                                              
      MaxNumConf=3                                            !Max number of conformers by torsion to look for 
                                                              
      AnalyticGrad=.false.                                    !Use analytical gradient for optimization, 
                                                              !the ".false." means the usage of linear
                                                              !gradient interpolation
!------------------------------------------

      NumArg=iargc()
      allocate(Arguments(NumArg))
      allocate(FileList(NumArg))

!---READING OPTIONS----------
      CALL Read_Options(FileList,UseCol,InfoFileName,DistMatrixName,&
           NumArg,NumFiles,Simplify,UseCut,k_value,AnalyticGrad,Step,&
           ConvCriterion,MaxIterations,MaxNumConf)


      print*," "
      print*," OPTIONS:"
      print*,"--------------------------------------------------"
      if(UseCut) then     
        write(*,'(A32)') "Distance Matrix file name:     "
        write(*,*) DistMatrixName
        write(*,'(A32)') "Dihedral Information file name:"
        write(*,*) InfoFileName
                            write(*,'(A25,A7)') "Using CutOff:           ","YES"
      else
                            write(*,'(A25,A7)') "Using CutOff:           ","NO"
      endif
                            write(*,'(A25,I7)') "Using column:           ",UseCol
                          write(*,'(A25,F7.1)') "k_value:                ",k_value
                            write(*,'(A25,I7)') "Max. num. of Iterations:",MaxIterations
                          write(*,'(A25,F7.1)') "Step:                   ",Step
                          write(*,'(A25,E7.1)') "Convergence Criterion:  ",ConvCriterion
                            write(*,'(A25,I7)') "Max. num. of Conformers:",MaxNumConf
      if(AnalyticGrad)      write(*,'(A25,A7)') "Using Analytic Gradient:","YES"
      if(.not.AnalyticGrad) write(*,'(A25,A7)') "Using Analytic Gradient:","NO"
      print*," "
      print*, " FILES:"

      A=0
      do i=1,NumFiles
         open(unit=i+6,file=FileList(i),status="OLD",&        !Just to check if the file exists
             iostat=ios)                                       
         if(ios.ne.0) then 
           print*," "
           write(*,'(2A22)') "ERROR during opening: ",FileList(i)
           stop
         endif
         ios=0
         if(i.eq.1) then
               do while(ios==0)
                 read(i+6,'(A1)',iostat=ios) charact          !Just to know the number of snapshots
                 if(ios>0) then
                  print*, "Error(0) while reading ",FileList(i)
                  stop
                 endif
                 if(ios==-1) exit
                 A=A+1
               enddo
         endif
         print*,i,FileList(i)
         close(i+6)
      enddo
      NumSnap=A
      print*,"--------------------------------------------------"
      print*," "
      print*,"Number of Snapshots: ",NumSnap
      print*,"Info: The Number of Snapshots is read from the first file: ",FileList(1)
      print*," "
!---------------------------


!--CREATING THE BIG MATRIX---------
      write(*,'(A49)') " Codifying and creating the full data matrix ... " 
      allocate(BigMatrix(NumSnap,NumFiles))
      allocate(ang(NumSnap))
      allocate(DiscreteAng(NumSnap))

      do i=1,NumFiles
         print*,FileList(i)
         open(unit=i+6,file=FileList(i),status="OLD")
               A=1
               ios=0
               do while(ios==0)
                 read(i+6,*,iostat=ios) (dummy, k=1,(UseCol-1)),ang(A)
                 do while((ang(A).lt.0).or.(ang(A).ge.360))
                  if(ang(A).lt.0.) ang(A)=ang(A)+360.            
                  if(ang(A).ge.360.) ang(A)=ang(A)-360.         
                 enddo
                 if(ios>0) then
                  print*, "Error while reading ",FileList(i)
                  stop
                 endif
                 if(ios==-1) exit
                 A=A+1
               enddo
               A=A-1
      
               if(NumSnap.ne.A) then
                print*," "
                print*, "ERROR: The files ",FileList(1)," and ",FileList(i),&
                        " have different number of snapshots or lines" 
                stop
               endif 
      
               close(i+6)
               call Get_Minimums(ang,NumSnap,k_value,&
                AnalyticGrad,Step,ConvCriterion,MaxIterations,&
                MaxNumConf,minimum,FNumMins)
      
               !-codify
               do k=1,NumSnap
                 if((FNumMins).eq.1) then
                   DiscreteAng(k)=1
                 else
                   do j=1,FNumMins
                     if(j.eq.1) then
                       if((ang(k).lt.minimum(j)).or.(ang(k).ge.minimum(FNumMins))) then
                        DiscreteAng(k)=j 
                       endif
                     else
                       if((ang(k).lt.minimum(j)).and.(ang(k).ge.minimum(j-1))) then
                        DiscreteAng(k)=j
                       endif
                     endif
                   enddo
                 endif
               enddo
               
               do k=1,NumSnap
                 BigMatrix(k,i)=DiscreteAng(k)
               enddo
      enddo
      print*,"OK"
      print*,"  "
!----------------------------------


!--SELECTING THOSE COLUMNS OF BigMatrix THAT CHANGE CONFORMATIONALLY
      write(*,'(A63)') " Selecting the columns to remove from the full data matrix ... " 
      allocate(selected(NumFiles))

      if(simplify.eq."yes") then
          C=0
          do j=1,NumFiles
             selected(j)=.false.
             do i=1,NumSnap
              if((i.gt.1).and.(BigMatrix(i,j).ne.BigMatrix(1,j))) then
                C=C+1
                selected(j)=.true.
                exit
              endif
             enddo
          enddo
          NewNumCol=C                                         !New number of columns after the elimination
      elseif(simplify.eq."no") then
          do j=1,NumFiles
             selected(j)=.true.
          enddo
          NewNumCol=NumFiles
      else
          print*,"ERROR: the option simplify must be yes/no"
          stop
      endif
      print*,"OK"
      print*,"  "
!--------------------------------------------------------------


!--CREATING THE "SMALL" MATRIX ELIMINATING THOSE COLUMNS WITH--
!--NO CONFORMATIONAL CHANGES-----------------------------------
      write(*,'(A65)')" Creating the small data matrix by removing constant columns ... " 
      allocate(Matrix(NumSnap,NewNumCol))
      C=0
      do j=1,NumFiles
         if(selected(j)) then
            C=C+1
            write(*,'(A3,A15,A14,I4)') "   ",FileList(j),"now in column ",C
            do k=1,NumSnap
             Matrix(k,C)=BigMatrix(k,j)                       !Copying in Matrix the heterogeneous columns
            enddo
         endif
         if(.not.selected(j)) write(*,'(A3,A15,A7)') "   ",FileList(j),"removed"
      enddo
      deallocate(BigMatrix)
      print*,"OK"
      print*,"  "
!--------------------------------------------------------------

   
!----OPENING AND READING THE DISTANCE MATRIX-----------
      if(UseCut) then
        write(*,'(A33)')" Reading the distance matrix ... "
       open(unit=9,file=DistMatrixName,status="OLD",iostat=ios)
          if(ios.ne.0) then
            print*," "
            write(*,'(2A22)') "ERROR during opening: ",DistMatrixName
            stop
          endif
        NumAtoms=0
        ios=0
        do while(ios==0)
           read(9,*,iostat=ios) dummy                         !Just to know the number of lines (ATOMS)
           if(ios>0) then
            print*, "Error while reading ",DistMatrixName
            stop
           endif
           if(ios==-1) exit
           NumAtoms=NumAtoms+1
        enddo
        rewind(9)
        allocate(Dist(NumAtoms,NumAtoms))
      
        DO i=1,NumAtoms
          READ(9,*)(Dist(i,j),j=1,NumAtoms)
        ENDDO
        close(9)
        print*,"OK"
        print*,"  "
      endif
!-----------------------------------------------------


!----OPENING AN READING THE INFO FILE-----
      if(UseCut) then
        write(*,'(A27)')" Reading the info file ... "
      
        allocate(Atoms_in_Dih(NewNumCol,2))
       open(unit=10,file=InfoFileName,status="OLD",iostat=ios)
          if(ios.ne.0) then 
            print*," "
            write(*,'(2A22)') "ERROR during opening: ",InfoFileName
            stop
          endif
      
        c=0
        do while(ios==0)
           read(10,*,iostat=ios) idummy                       !Just to know the number of lines
           if(ios>0) then
            print*, "Error while reading ",DistMatrixName
            stop
           endif
           if(ios==-1) exit
           c=c+1
        enddo
        rewind(10)
        if(c.gt.NumFiles) then
          print*, "WARNING: The number of lines in the info file is greater"
          print*, "         than the number of dihedrals (files)"
          print*, " "
        elseif(c.lt.NumFiles) then
          print*, "ERROR: The number of dihedrals (files) is greater than"
          print*, "       the number of lines in the info file, the program"
          print*, "       can not countinue using cut-off."
          print*, " "
          stop
        endif
      
       c=0
        do i=1,NumFiles
          if(selected(i)) then
            c=c+1
            read(10,*,iostat=ios) (Atoms_in_Dih(c,j), j=1,2)
          else
            read(10,*,iostat=ios) idummy,idummy2
          endif
          if(ios.ne.0) then
            print*, "Error while reading ",InfoFileName
            stop
          endif
        enddo
      
        close(10)
        print*,"OK"
        print*,"  "
      endif
!---------------------------------------


!---MAKING THE NEW (OR REDUCED) DISTANCE MATRIX--- 
      if(UseCut) then
        write(*,'(A49)')" Making the new (or reduced) distance matrix ... "
        allocate(NewDistMat(NewNumCol,NewNumCol))
        do i=1,NewNumCol
          do j=1,NewNumCol
             if(i.eq.j) then
              NewDistMat(i,j)=0.0
             else
              NewDistMat(i,j)=(1./4.)*( Dist(Atoms_in_Dih(i,1),Atoms_in_Dih(j,1)) &
                                      + Dist(Atoms_in_Dih(i,2),Atoms_in_Dih(j,2)) &
                                      + Dist(Atoms_in_Dih(i,1),Atoms_in_Dih(j,2)) &
                                      + Dist(Atoms_in_Dih(i,2),Atoms_in_Dih(j,1)))
             endif
          enddo
        enddo
        print*,"OK"
        print*,"  "
      endif
!------------------------------------------------- 


!----PRINTING OUTPUT------
      !writing the reduced distance matrix
      if(UseCut) then
        write(*,'(A40)')" Writing in reduced_dist_matrix.dat ... "
        open(11,file="reduced_dist_matrix.dat")
     
        do I=1,NewNumCol
          write(11,'(10000F9.3)')(NewDistMat(i,j),j=1,NewNumCol)
        enddo
        close(11)
        print*,"OK"
        print*,"  "
      endif
     
      !writing the final matrix
      write(*,'(A26)')" Writing in MATRIX.dat ... "
      open(12,file="MATRIX.dat")
      do i=1,NumSnap
        write(12,'(10000I2)') (Matrix(i,j),j=1,NewNumCol)
      enddo 
      close(12)
      print*,"OK"
      print*,"  "
!--------------------------
        
      STOP " DONE!"
!----------------------------------------------------------------------------------------
      END PROGRAM 
!----------------------------------------------------------------------------------------



!      _____________________________
!     /|                           |
!     /| SUBROUTINES AND FUNCTIONS |
!     /|___________________________|
!     //////////////////////////////


!*****************************************************************************************
      SUBROUTINE Read_Options(FileList,UseCol,InfoFileName,DistMatrixName,&
                 NumArg,c,Simplify,UseCut,k_value,AnalyticGrad,Step,&
                 ConvCriterion,MaxIterations,MaxNumConf)
!*****************************************************************************************
!This subroutine reads all the options that can be given to the 
!program through the command line
!-----------------------------------------------------------------------------------------
      use parameters
      implicit none
      !-VARIABLE DEFINITIONS--------------------------------------------------------------
      integer MaxNumConf                                      !Max number of conf. states by torsion allowed
      logical GivenInfoFileName                               !The info filename is given? (True/False)
      logical GivenDistMatrixName                             !The distance matrix filename is given? (T/F)
      logical UseCut                                          !Use cutoff criterion? (True/False)
      character*60 DistMatrixName                             !Distance matrix filename
      character*60 InfoFileName                               !Info-file filename
      integer NumArg                                          !Number of arguments in the command line
      character*60 Arguments(NumArg)                          !Arguments in the command line
      character*60 FileList(NumArg)                           !List of data files (*.dat) 
      character*3 Simplify                                    !Eliminate rigid torsions? (True/False)
      real(DP) k_value                                        !k_value (see default values in the main program)
      character*3 answer                                      !Use analytic gradient? (yes/no)
      logical  AnalyticGrad                                   !The logical analogous to "answer" (True/False)
      real(DP) Step                                           !Step size in the optimization
      integer  MaxIterations                                  !Max number of iterations in the optimization
      real(DP) ConvCriterion                                  !Convergence criterion for the optimization
      integer UseCol                                          !Column to be used for the discretization in the 
                                                              !data files 

      real(DP)  RealVar                                       !    
      integer i,c                                             !
      integer ios !                                            >  Auxiliary or Dummy variables 
      logical EO                                              !
      character*60 arg                                        !
      !-----------------------------------------------------------------------------------


!-----INITIAL VALUES
      GivenInfoFileName=.false.
      GivenDistMatrixName=.false.
      EO=.false.                                              !End Options is False

!--READING OPTIONS
      do i = 1, NumArg
        call getarg(i, arg)
        Arguments(i)=arg
      enddo
       
      i=1
      c=0
      do while (i.le.NumArg)
        if((Arguments(i).eq.'-u').or.(Arguments(i).eq.'-usecol').and.(.not.EO)) then
          read(Arguments(i+1),*,iostat=ios) UseCol
          read(Arguments(i+1),*,iostat=ios) RealVar
          if((ios.ne.0).or.(UseCol.lt.1).or.(UseCol.ne.RealVar).and.(.not.EO)) & 
           stop "ERROR: Check the -u/-usecol option. Use -help option for quick help"
          i=i+1
        elseif((Arguments(i).eq.'-dist').and.(.not.EO)) then
          read(Arguments(i+1),'(A60)',iostat=ios) DistMatrixName 
          if(ios.ne.0) & 
          stop "ERROR: Check the file name of distance matrix. Use -help option for quick help"
          i=i+1
          GivenDistMatrixName=.true.
        elseif((Arguments(i).eq.'-i').or.(Arguments(i).eq.'-info').and.(.not.EO)) then
          read(Arguments(i+1),'(A60)',iostat=ios) InfoFileName
          if(ios.ne.0) stop "ERROR: Check the info file name. Use -help option for quick help"
          i=i+1
          GivenInfoFileName=.true.
        elseif((Arguments(i).eq.'-nocut').and.(.not.EO)) then
          UseCut=.false.
        elseif((Arguments(i).eq.'-s').or.(Arguments(i).eq.'-simplify').and.(.not.EO)) then
          read(Arguments(i+1),*,iostat=ios) Simplify
          if((ios.ne.0).or.(.not.((Simplify.eq."yes").or.(Simplify.eq."no")))) & 
             stop "ERROR: Check the option simplyfy. Use -help option for quick help"
          i=i+1
          GivenInfoFileName=.true.
        elseif((Arguments(i).eq.'-k').and.(.not.EO)) then
          call getarg(i+1, arg)
          read(arg,*,iostat=ios) k_value
          if(ios>0) stop "ERROR: Check k_value. Use -help option for quick help"
          if(k_value.le.0) stop "ERROR: k_value most be a positive real"
          i=i+1
        elseif((Arguments(i).eq.'-ag').and.(.not.EO)) then
          read(Arguments(i+1),*,iostat=ios) answer
          if((ios.ne.0).or.(.not.((answer.eq."yes").or.(answer.eq."no")))) then
             print*,"ERROR: Check the option -ag. "
             print*,"Use -help option for quick help"
             stop
          endif
          if(answer.eq."yes") AnalyticGrad=.true.
          i=i+1
        elseif((Arguments(i).eq.'-step').and.(.not.EO)) then
          call getarg(i+1, arg)
          read(arg,*,iostat=ios) Step
          if(ios>0) stop "ERROR: Check Step value. Use -help option for quick help"
          if(Step.le.0.0) stop "ERROR: Step most be a positive real"
          i=i+1
        elseif((Arguments(i).eq.'-crit').and.(.not.EO)) then
          call getarg(i+1, arg)
          read(arg,*,iostat=ios) ConvCriterion
          if(ios>0) stop "ERROR: Check Convergence Criterion value. Use -help option for quick help"
          if(ConvCriterion.le.0) stop "ERROR: Convergence Criterion most be a positive real"
          i=i+1
        elseif((Arguments(i).eq.'-maxitr').and.(.not.EO)) then
          read(Arguments(i+1),*,iostat=ios) MaxIterations
          read(Arguments(i+1),*,iostat=ios) RealVar
          if((ios.ne.0).or.(MaxIterations.lt.1).or.(MaxIterations.ne.RealVar)) & 
           stop "ERROR: Check the -maxitr option. Use -help option for quick help"
          i=i+1
        elseif((Arguments(i).eq.'-maxconf').and.(.not.EO)) then
          read(Arguments(i+1),*,iostat=ios) MaxNumConf
          read(Arguments(i+1),*,iostat=ios) RealVar
           if((ios.ne.0).or.(MaxNumConf.lt.2).or.&
              (MaxNumConf.ne.RealVar).or.(MaxNumConf.gt.9)) then 
             print*,"ERROR: Check the -maxconf option. "
             print*,"The maximum number of conformers to select from"
             print*,"each dihedral, most be an integer number between 2 and 9."
             print*,"Use -help option for quick help"
             stop
           endif
           i=i+1
        elseif((Arguments(i).eq.'-help').or.(Arguments(i).eq.'-h').and.(.not.EO)) then
          print*,"QUICK HELP:"
          print*,""
          print*,"USAGE: centro_prep [OPTIONS] file1 file2 ..."
          print*," "
          print*,"The files file1, file2 etc, contain the time series of the dihedral angles."
          print*,"It is also valid the usage of regular expressions, for example:"
          print*,"             "
          print*," centro_prep [OPTIONS] file?"
          print*," centro_prep [OPTIONS] fi*"
          print*,"             "
          print*,"OPTIONS:"
          print*,"--------"
          print*," -u/-usecol NUMCOL                             Default:2"
          print*,"            The number NUMCOL specifies which column of file1, file2 etc." 
          print*,"            contains the dihedral angle variable. This value is normally 2"
          print*,"            since the first column is often the time or the snapshot number"
          print*,"      "
          print*," -dist      DISTANCE_MATRIX_FILE_NAME          Default:distance_matrix.dat" 
          print*,"            Specifies the distance matrix file name of the full system in" 
          print*,"            vacuo."
          print*,""
          print*," -i/-info   DIH_INFO_FILE_NAME                 Default:atoms_in_dih.info" 
          print*,"            Specifies the file name that contains which atoms are involved"
          print*,"            in each dihedral(the two central ones). For example, if the "
          print*,"            first row is 3 4, means that the first dihedral (the one whose"
          print*,"            time series is in file1) is defined over the bond 3-4."
          print*,""
          print*," -nocut                                        Default:Use cut-off" 
          print*,"            Using this option no cut-off will be apply and the options" 
          print*,"            (-dist/-info) and files DISTANCE_MATRIX_FILE_NAME and"
          print*,"            DIH_INFO_FILE_NAME  are not needed." 
          print*,""
          print*," -s/-simplify  yes/no                          Default:yes" 
          print*,"            Remove the dihedrals without conformational changes"
          print*,""
          print*," -k         K_VALUE                            Default: 1.0"
          print*,"            The k_value sets the smoothing parameter v " 
          print*,"            in the von-Mises kernel density estimation as"
          print*,"            was proposed in eq.(7) of ref: "
          print*,"            Computational Statistics & Data Analysis"
          print*,"            Volume 52, Issue 7, 15 March 2008, Pages 3493-3500."
          print*,"            We set by default k_value=1.0 because it slightly"
          print*,"            oversmooth the PDFs, which is convenient for minimizations."
          print*,""
          print*," -ag        yes/no                             Default: no"
          print*,"            Use analytic gradient on each step of the minimizations"
          print*,"            that look for the minimums in the Probability Density"
          print*,"            Functions (PDFs). The default option (no) use a quite"
          print*,"            accurate and fast linear interpolation of the gradient."
          print*,""
          print*," -step      STEP_SIZE                          Default: 5 (degrees)"
          print*,"            Step size a minimization of the form:"
          print*,"            X(n+1)=X(n)-STEP_SIZE*GRADIENT"
          print*,""
          print*," -crit      CONVERGENCE_CRITERION              Default: 1.0E-4"
          print*,"            Convergence criterion in the minimization"
          print*,""
          print*," -maxitr    MAX_NUMBER_OF_ITERATIONS           Default: 1000"
          print*,"            Maximum number of iteration on each minimization"
          print*,""
          print*," -maxconf   MAX_NUMBER_OF_CONFORMERS           Default: 3"
          print*,"            Maximum number of conformers to look for on each"
          print*,"            dihedral angle"
          print*,""
          print*," -help"
          print*,"            Open this quick help"
          print*,""
          print*,"OTHER EXAMPLES:"
          print*,"---------"
          print*,"centro_prep -help"
          print*,"centro_prep f*.out"
          print*,"centro_prep -ag yes -k 1.5 f*.out"
          print*,"centro_prep -nocut f*.out"
          print*,"centro_prep -u 1 -s yes file?.out"
          print*,"centro_prep -u 1 -dist distance_matrix.dat -i atoms_in_dih.info file?.out"
          stop
        else
          EO=.true.
          arg=Arguments(i)
          if(arg(1:1).eq."-") &
            stop "ERROR: Check the options, use -help option for quick help"
          c=c+1
          FileList(c)=arg
        endif
        i=i+1
      enddo
      
      return
      END SUBROUTINE Read_Options
!----------------------------------------------------------------------------------------



!*****************************************************************************************
      SUBROUTINE Get_Minimums(ang,NumSnap,k_value,&
                 AnalyticGrad,Step,ConvCriterion,MaxIterations,&
                 MaxNumConf,minimum,FNumMins)
!*****************************************************************************************
! Get the minimums from the probability density function
!----------------------------------------------------------------------------------------
      use parameters
      implicit none
      !-VARIABLE DEFINITIONS--------------------------------------------------------------
      integer ,intent(in) :: NumSnap                          !Number of snapshots
      real(DP),intent(in) :: ang(NumSnap)                     !Dihedral angle
      real(DP),intent(in) :: k_value                          !k_value (see default values in the main program)
      logical ,intent(in) :: AnalyticGrad                     !Use analytic gradient? (true/false) 
      real(DP),intent(in) :: Step                             !Step for the steepest-descendent optimization
      real(DP),intent(in) :: ConvCriterion                    !Convergence criterion for the steepest-descendent opt.
      integer ,intent(in) :: MaxIterations                    !Max number of iterations in the optimization
      integer ,intent(in) :: MaxNumConf                       !Max number of conformers allowed by torsion 
      real(DP),intent(out) :: minimum(9)                      !Position of the minimums in the PDF
      integer ,intent(out) :: FNumMins                        !Final number of minimums founded
      real(DP) CoordMaxs(9)                                   !Positions of the maximums in the PDF
      real(DP) SecDerivMaxs(9)                                !Positions of the maximums of the 2nd derivate of the (PDF)
      real(DP) BessI                                          !Modified Bessel function 
      real(DP) density(0:360)                                 !Probability Density Function (PDF)
      real(DP) gradient(0:360)                                !Gradient of the PDF
      real(DP) secderiv(0:360)                                !Second derivate of the PDF
      real(DP) SmoothingParam                                 !Smoothing parameter (see ref. 2)
      real(DP) vMgradient                                     !Gradient in a general position (not only in {0,...359})
      real(DP) vMdensity                                      !Density in a general position 
      real(DP) vMsec_deriv                                    !Second derivate in a general position 
      real(DP) Func_vMgradient                                !Function 
      real(DP) Step_times_grad                                !Step x gradient
      integer ObsNmax                                         !Observed Number of Maximums
                                                              
      integer contador                                        !
      integer i,j,Itrs,m                                      !
      real(DP) temp                                           !> Auxiliary or dummy variables
      logical Converged                                       !
      real(DP) grad                                           !
      !-----------------------------------------------------------------------------------
      
     secDerivMaxs=0

!----Computing the Bandwith
      SmoothingParam=(3*NumSnap*(k_value**2)*BessI(2,2*k_value)/&
                (4*(Pi**0.5)*(BessI(0,k_value))**2.) )**(2./5.)


!----Calculating densities, gradient, second derivate
!----and looking for local maximums
      ObsNmax=0    ! Observed Number of Maximums
      do i=0,360
        call vMises(NumSnap,SmoothingParam,(i*Pi/180.),ang,vMdensity,vMgradient,vMsec_deriv)
        gradient(i)=vMgradient
         density(i)=vMdensity
        secderiv(i)=vMsec_deriv
        if(i.gt.0) then
          if(((gradient(i-1)*gradient(i)).le.0.).and.(vMsec_deriv.lt.(-1.0E-4)).and.(ObsNmax.lt.9)) then
           ObsNmax=ObsNmax+1     
           CoordMaxs(ObsNmax)=i/2.+(i-1)/2.                   !Approximate coord of the maximum
           SecDerivMaxs(ObsNmax)=&
               secderiv(i)/2. + secderiv(i-1)/2.              !Approximate second derivate
          endif
        endif
      enddo

!----Sorting to give preference to those maximums with
!----lower(more negative) second derivate
      do i=1,ObsNmax-1
         do j=i+1,ObsNmax
         if(SecDerivMaxs(i).gt.SecDerivMaxs(j)) then
          temp=SecDerivMaxs(i)
          SecDerivMaxs(i)=SecDerivMaxs(j)
          SecDerivMaxs(j)=temp
          temp=CoordMaxs(i)
          CoordMaxs(i)=CoordMaxs(j)
          CoordMaxs(j)=temp
         endif
         enddo
      enddo

!----In a circular variable, the Max. number of minimums 
!----most be equal to the Max. number of maximums, we get 
!----the lowest value between the predefined MaxNumConf 
!----and the observed number of maximums/minimums (ObsNmax) 
      if(MaxNumConf.lt. ObsNmax)  FNumMins=MaxNumConf
      if(MaxNumConf.ge. ObsNmax)  FNumMins=ObsNmax

!----Sorting only the first FNumMins using the coordinates
      do i=1,FNumMins-1
        do j=i+1,FNumMins
         if(CoordMaxs(i).gt.CoordMaxs(j)) then
          temp=CoordMaxs(i)
          CoordMaxs(i)=CoordMaxs(j)
          CoordMaxs(j)=temp
        endif 
        enddo
      enddo

!----Guessing where are the minimums
      if (FNumMins.eq.1) then
        minimum(1)=(CoordMaxs(1)-180.)
        if (minimum(1).lt.0.0) minimum(1)=minimum(1)+360.
      else
        do i=1,FNumMins
          if (i.eq.1) then
           minimum(i)=(CoordMaxs(i)+CoordMaxs(FNumMins)-360.)/2.
           if (minimum(i).lt.0.0) minimum(i)=minimum(i)+360.
          else
           minimum(i)=(CoordMaxs(i)+CoordMaxs(i-1))/2.
          endif
        enddo
      endif

!---Optimizing Minimums 
      if(FNumMins.gt.1) then
       do m=1,FNumMins
        Itrs=0
        Converged=.false.
        do while((.not.Converged).and.(Itrs<MaxIterations))
          Itrs=Itrs+1
           if(AnalyticGrad) then
             grad=Func_vMgradient(NumSnap,SmoothingParam,(Pi/180.)*minimum(m),ang)
           else
             grad=(abs(minimum(m)-1-int(minimum(m))))&        !Gradient Linear interpolation
                *gradient(int(minimum(m)))+&
               (abs(minimum(m)-int(minimum(m)))&
                *gradient(1+int(minimum(m))))
           endif
          if ((abs(grad).lt.ConvCriterion)) then
           Converged=.true.
          else
           if(abs(Step*grad).gt.5.0) then
            Step_times_grad=((Step*grad)/abs(Step*grad))*5.0
           else
            Step_times_grad=Step*grad
           endif
           minimum(m)=minimum(m)-Step_times_grad
           if(minimum(m).lt.0.) minimum(m)=minimum(m)+360.
           if(minimum(m).ge.360.) minimum(m)=minimum(m)-360.
          endif
        enddo
       enddo
      endif

      do i=1,FNumMins-1
        do j=i+1,FNumMins
          if(minimum(i).gt.minimum(j)) then
            temp=minimum(i)
            minimum(i)=minimum(j)
            minimum(j)=temp 
          endif
        enddo      
      enddo

      write(*,'(A40)') "minimums coord. in degrees:"
      do i=1,FNumMins
       write(*,'(F40.2)') minimum(i)
      enddo
      
      return


      END SUBROUTINE Get_Minimums
!----------------------------------------------------------------------------------------



!*****************************************************************************************
      SUBROUTINE vMises(NumSnap,SmoothingParam,X,ang,vMdensity,vMgradient,vMsec_deriv)
!*****************************************************************************************
!Computes the density, gradient and second derivate in a position X
!-----------------------------------------------------------------------------------------
      use parameters
      implicit none
      !-VARIABLE DEFINITIONS--------------------------------------------------------------
      integer NumSnap                                         !Number of snapshots
      real(DP) vMdensity                                      !von Mises density
      real(DP) vMgradient                                     !von Mises gradient
      real(DP) vMsec_deriv                                    !von Mises second derivate
      real(DP) BessI                                          !Modified Bessel function
      real(DP) SmoothingParam                                 !Smoothing parameter
      real(DP) ang(*)                                         !Data Angles in degrees
      real(DP) aj                                             !Data Angles in radians
      real(DP) expo                                           !Auxiliary real variable
      integer j                                               !Dummy integer variable
      real(DP), intent(in) :: X                               !The position where we want to
                                                              !evaluate vMdensity,vMgradient,...
      !-----------------------------------------------------------------------------------

      vMgradient=0.0
      vMgradient=0.0
      vMsec_deriv=0.0
      do j=1,NumSnap
       aj=ang(j)*(Pi/180.)                                    !transforming to radians
       expo=exp( SmoothingParam*cos(X-aj))
       vMdensity=vMdensity+expo
       vMgradient=vMgradient-sin(X-aj)*expo
       vMsec_deriv=vMsec_deriv+&
                   (SmoothingParam*(sin(X-aj))**2.-cos(X-aj))*expo
      enddo
       vMdensity = vMdensity/( 2.0*Pi*NumSnap*BessI(0,SmoothingParam))
       vMgradient= vMgradient*SmoothingParam/( 2.0*Pi*NumSnap*BessI(0,SmoothingParam))
       vMsec_deriv = vMsec_deriv*SmoothingParam/( 2.0*Pi*NumSnap*BessI(0,SmoothingParam))
      return
      END SUBROUTINE vMises
!----------------------------------------------------------------------------------------



!*****************************************************************************************
      FUNCTION Func_vMgradient(NumSnap,SmoothingParam,X,ang)
!*****************************************************************************************
!Function specifically designed to compute the von Misses Gradient
!----------------------------------------------------------------------------------------
      use parameters 
      implicit none
      !-VARIABLE DEFINITIONS--------------------------------------------------------------
      real(DP) Func_vMgradient                                ! 
      real(DP) BessI                                          !Modified Bessel function
      real(DP) SmoothingParam                                 !Smoothing parameter
      real(DP) ang(*)                                         !Data Angle in degrees
      real(DP) aj                                             !Data Angle in radians
      integer NumSnap                                         !Number of snapshots
      integer j                                               !
      real(DP), intent(in) :: X                               !The position where we want to
                                                              !evaluate vMdensity,vMgradient,...
      !-----------------------------------------------------------------------------------
      Func_vMgradient=0.0
      do j=1,NumSnap
       aj=ang(j)*(Pi/180.)                                    !transforming to radians
       Func_vMgradient=Func_vMgradient-SmoothingParam*sin(X-aj)*&
       exp( SmoothingParam*cos(X-aj))/&
       ( 2.0*Pi*NumSnap*BessI(0,SmoothingParam))
      enddo      
      return
      END FUNCTION 
!----------------------------------------------------------------------------------------



!*****************************************************************************************
      FUNCTION BessI(N,X)
!*****************************************************************************************
!     Calc. the first kind modified Bessel function
!
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!-----------------------------------------------------------------------------------------
      use parameters 
      parameter (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      real(DP) X,BessI,BessI0,BessI1,TOX,BIM,BI,BIP

      if (N.EQ.0) then
        BessI = BessI0(X)
        return
      endif

      if (N.EQ.1) then
        BessI = BessI1(X)
        return
      endif

      if(X.EQ.0.D0) then
        BessI=0.D0
        return
      endif

      TOX = 2.D0/X
      BIP = 0.D0
      BI  = 1.D0
      BessI = 0.D0
      M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))

      do J = M,1,-1
        BIM = BIP+DFLOAT(J)*TOX*BI
        BIP = BI
        BI  = BIM
        if (ABS(BI).GT.BIGNO) then
          BI  = BI*BIGNI
          BIP = BIP*BIGNI
          BessI = BessI*BIGNI
        endif
        if (J.eq.N) BessI = BIP
      enddo

      BessI = BessI*BessI0(X)/BI

      return

      END FUNCTION BessI
!----------------------------------------------------------------------------------------



!*****************************************************************************************
      FUNCTION BessI0(X)
!*****************************************************************************************
! Auxiliary Bessel functions for N=0
!-----------------------------------------------------------------------------------------
      use parameters 
      real(DP) X,BessI0,Y,P1,P2,P3,P4,P5,P6,P7,  &
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      data P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
      0.2659732D0,0.360768D-1,0.45813D-2/
      data Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
      0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
      0.2635537D-1,-0.1647633D-1,0.392377D-2/

      if(ABS(X).LT.3.75D0) then
        Y=(X/3.75D0)**2
        BessI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      else
        AX=ABS(X)
        Y=3.75D0/AX
        BX=EXP(AX)/SQRT(AX)
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
        BessI0=AX*BX
      endif

      return
      END FUNCTION BessI0
!-----------------------------------------------------------------------------------------



!*****************************************************************************************
      FUNCTION BessI1(X)
!*****************************************************************************************
! Auxiliary Bessel functions for N=1
!-----------------------------------------------------------------------------------------
      use parameters
      real(DP) X,BessI1,Y,P1,P2,P3,P4,P5,P6,P7,  &
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      data P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
      0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      data Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
      -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
      -0.2895312D-1,0.1787654D-1,-0.420059D-2/

      if(ABS(X).LT.3.75D0) then
        Y=(X/3.75D0)**2
        BessI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      else
        AX=ABS(X)
        Y=3.75D0/AX
        BX=EXP(AX)/SQRT(AX)
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
        BessI1=AX*BX
      endif

      return
      END FUNCTION BessI1
!-----------------------------------------------------------------------------------------
     

