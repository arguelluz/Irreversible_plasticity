﻿# Irreversible_plasticity
Irreversible_plasticity project A.Rago and M.Brun-Usan

## Compiling files
Uses f90 FORTRAN compiler, available by default on most LINUX distributions  
All commands for compiling are in grns.sh  
grns.sh also includes routines to check correct syntax  

## Workflow
Change parameters in file "start"  
Run models from file "grns2"  

Requires 3 main files:  

* start: Declares most variables used in simuations
* development: Contains developmental loop for each individual
* grns2: Contains evolutionary loop (fitness evaluation, mutation etc)
Calls start for global parameters and development for each individual
Loads pre-trained networks named in the GRNames.dat file (for now, working to include it as a parameter in "start")

## Output file format
Running GRN2 generates files in the following format:  
GRN_TARGET_THRESHOLD_REPLICATE_TIME_TARGET1_REPLICATE1_TIME1.dat  
Where GRN is the file identifier  
TARGET is a string of integers which identifies the multivariate targets in format E1T1E1T2ENT1ENT2 
THRESHOLD is an integer stating in which cell the target E1T1E1T2 switches to ENT1ENT2
TIME is an integer describing which simulation checkpoint is recored  
REPLICATE is an integer which denotes the unique ID of the random replicates (random seed)  
The latter part of the string contains the name of the file which we load the initial population from, or BLANK in case the simulation is initiated from empty networks  

Content is formatted as a ascii file with commented header which records the following simulation parameters:  
     write(7000,*)'TARGETS (E1T1,E1T2,ENT1,ENT2)',block(1:2,1), block(1:2,n)  ! 1
     write(7000,*)'THRESHOLDS(CELL).............',thresholds(1)               ! 2
     write(7000,*)'POPULATON SIZE...............',p                           ! 3
     write(7000,*)'STRENGHT OF SELECTION........',ss                          ! 4
     write(7000,*)'RECOMBINATION; 1=YES; 0=NO   ',reco                        ! 5
     write(7000,*)'TRAINING (1) vs TEST (0) SET ',training                    ! 6
     write(7000,*)'NUMBER /  TOTAL REPLICATES...',replica,replicas            ! 7
     write(7000,*)'CURRENT VS MAXIMUM GENERATION',et,etmax,lapso              ! 8
     write(7000,*)'ENV. FACTORS/ENVIRONMENTS....',EF,n                        ! 9
     write(7000,*)'NUMBER GENES,PHEN. DIMENSIONS',ng,PD                       ! 10
     write(7000,*)'TMAX,SDEV,SS,RECO,CAPPED.....',tmax,sdev,ss,reco,capped    ! 11   
     write(7000,*)'CONNECTIVITIES WW / MZZ .....',conWW,conMZZ                ! 12   
Following the header the file records, for each individual, the W and WW matrices, encoding the connection weights of the GRN and the binary matrix detremining the GRN non-zero connections.

Connection matrices are stacked across all individuals recorded in the simulation (rbind).  
At the end of the file there are written the (non-mutable) MZ and MZZ matrices, which are the same for all individuals.
Every matrix written is tab separated.

## File storing.
 
The program uses/requires a folder called "files" located in the same directory. It places there automatically the generated files and takes them from there in the test set.
