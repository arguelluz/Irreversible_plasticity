#!/bin/bash

 #gfortran -w -fexceptions -fno-underscoring -fbounds-check inicio.mod.f90 diffusion.mod.f90 reaction.mod.f90 gnuplotter.mod.f90 idps.f90 -o idps.e 

 gfortran -w -fexceptions -fno-underscoring -fcheck=all -Wall -Wtabs start.mod.f90 development.mod.f90 grns2.f90 -o grns.e
