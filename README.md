# Irreversible_plasticity
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
GRN_TARGET_REPLICATE_TIME_TARGET1_REPLICATE1_TIME1.dat  
Where GRN is the file identifier  
TARGET is a string of integers which identifies the multivariate targets in format E1T1E1T2E2T1E2T2  
TIME is an integer describing which simulation checkpoint is recored  
REPLICATE is an integer which denotes the unique ID of the random replicates (random seed)  
The latter part of the string contains the name of the file which we load the initial population from, or BLANK in case the simulation is initiated from empty networks  

Content is formatted as a ascii file with commented header which records the following simulation parameters:  
* Phenotypic targets
* Population size
* Strength of selection
* Recombination or not
* Recurrent network or not
* Linear activation function or not
* Pre-trained network source (NA if naive network)
* Replicates
* Current generation (et)
* Total length of simulation (etmax)

Following the header the file records the connection weights of the (non-mutable) output layer and following that the recurrent layer of each individual in the population as tab separated matrices.  
Connection matrices are stacked across all individuals recorded in the simulation (rbind).  
