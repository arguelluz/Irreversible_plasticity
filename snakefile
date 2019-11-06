
# Define input trained GRNs
trained_grns, = glob_wildcards("files/GRN___________________{base}.dat")
grn_path = expand("files/GRN___________________{base}.dat", base = trained_grns)

modules = ("development.mod.f90", "grns2.f90")
problems, = glob_wildcards("./start_{problems}.f90")

# Run initial problem set on naive networks
rule train:
    input:
        modules = modules,
        problem_files = expand("start_{problems}.f90", problems = problems)
    params:
        problem_name = problems
    shell:
        '''
        for problem in {params.problem_name}
        do

        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs {input.modules} start_$problem.f90 -o grns.e &&

        ./grns.e &&

        mv ./files/GRN*.dat ../Simulation_results/test/$problem/
        mv ./files/PHE*.dat ../Simulation_results/test/$problem/

        rm *.mod

        done
        '''

# Run trained networks on all problems
rule test:
    input:
        modules = modules,
        grns = grn_path
    shell:
        '''
        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs {input} -o grns.e &&
        ./grns.e &&
        mv ./files/GRN*.dat ../Simulation_results/test/
        mv ./files/PHE*.dat ../Simulation_results/test/
        '''
