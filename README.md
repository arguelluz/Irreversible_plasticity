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

## Output file format
