del *.dat
del *.exe
del *.mod
del GRNstatus.txt

gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs start.mod.f90 development.mod.f90 grns2.f90 -o grns.exe
