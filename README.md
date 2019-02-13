# Irreversible_plasticity
Irreversible_plasticity project  
A. Rago and M. Brun-Usan

## Compiling files
Uses f90 FORTRAN compiler, available by default on most LINUX distributions.  
All commands for compiling are in grns.sh  
grns.sh also includes routines to check correct syntax  

## Workflow
Change parameters in file start.mod.90  
Recompile program by running grns.sh  
Run models using program grns.e  

## Program structure
Requires 3 main files:  

* start: Declares most variables used in simulation
* development: Contains developmental loop for each individual
* grns2: Contains evolutionary loop (fitness evaluation, mutation etc)

Calls start for global parameters and development for each individual.

## Manual input options

If option training = 0, loads pre-trained networks specified by name in start file (lines 92+).  
If option mzadhoc = 1, loads manually specified hidden to phenotype matrices from /files/mzadhoc.dat  
Matrix format is rows = genes, columns = traits. The first rows specify the continuous connection weights, followed by the discrete presence/absence connection data.  

## Output file format
The program produces 2 types of files:  
Population files contain the machine-readable simulation parameters and the matrices necessary to re-load populations.  
Phenotype files contain the simulation results for every generation.

### Population files
Running GRN2 generates files in the following format:  
GRN_TARGET_THRESHOLD_REPLICATE_TIME_TARGET1_REPLICATE1_TIME1.dat  
Where GRN is the file identifier  
TARGET is a string of integers which identifies the multivariate targets in format E1T1E1T2ENT1ENT2  
THRESHOLD is an integer stating in which cell the target E1T1E1T2 switches to ENT1ENT2  
TIME is an integer describing which simulation checkpoint is recorded  
REPLICATE is an integer which denotes the unique ID of the random replicates (random seed)  
The latter part of the string contains the name of the file which we load the initial population from, or BLANK in case the simulation is initiated from empty networks  

Content is formatted as a ascii file with commented header which records the following simulation parameters:  

```
     'TARGETS (E1T1,E1T2,ENT1,ENT2)',block(1:2,1), block(1:2,n)  ! 1  
     'THRESHOLDS(CELL).............',thresholds(1)               ! 2  
     'POPULATON SIZE...............',p                           ! 3  
     'STRENGHT OF SELECTION........',ss                          ! 4  
     'RECOMBINATION; 1=YES; 0=NO   ',reco                        ! 5  
     'TRAINING (1) vs TEST (0) SET ',training                    ! 6  
     'NUMBER /  TOTAL REPLICATES...',replica,replicas            ! 7  
     'CURRENT VS MAXIMUM GENERATION',et,etmax,lapso              ! 8  
     'ENV. FACTORS/ENVIRONMENTS....',EF,n                        ! 9  
     'NUMBER GENES,PHEN. DIMENSIONS',ng,PD                       ! 10  
     'TMAX,SDEV,SS,RECO,CAPPED.....',tmax,sdev,ss,reco,capped    ! 11   
     'CONNECTIVITIES WW / MZZ .....',conWW,conMZZ                ! 12   
```

Following the header the file records, for each individual, the W (weights) and WW (discrete connection) matrices for the hidden layer. Connection matrices are stacked across all individuals recorded in the simulation (rbind).  
At the end of the file there are written the (non-mutable) MZ and MZZ matrices, which are the same for all individuals.
Every matrix written is tab separated.

### Phenotype files
Phenotype files are tab separated annotations of the phenotype of every individual of the population, at each time step and in every environment.  
Each row represents an individual of a population at a given time step in a given environment. NOTE: Individuals change between generations: individual IDs are only used to distinguish between individuals of the same generation.   
Columns are annotated as follows:  

```
"Replicate"      Integer, ID of the simulation replicate  
"Generation"     Integer, Number of generations since the simulation start  
"Individual"     Integer, ID of the individual in the population.  
"Environment"    Integer, ID of the environment the organism is exposed to  
"Trait"          Integer, ID of the trait whose value is shown in column "Phenotype"  
"Phenotype"      Continuous real, numeric value of the trait annotated in column "Trait"  
"Fitness"        Continuous real, average fitness value of the individual across ALL environments
```

## File storing.

The program uses/requires a folder called "files" located in the same directory. It places there automatically the generated files and takes them from there in the test set.
