
# Define input trained GRNs
trained_grns, = glob_wildcards("files/GRN___________________{base}.dat")
grn_path = expand("files/GRN___________________{base}.dat", base = trained_grns)
modules = ("start.mod.f90", "development.mod.f90", "grns2.f90")
# Run initial problem set on naive networks


# Run trained networks on all problems
rule test:
    input:
        modules
    shell:
        '''
        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs {input} -o grns.e &&
        ./grns.e
        '''
