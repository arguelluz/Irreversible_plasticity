
# Define input trained GRNs
trained_grns, = glob_wildcards("files/GRN___________________{base}.dat")
grn_path = expand("files/GRN___________________{base}.dat", base = trained_grns)

modules = ("development.mod.f90", "grns2.f90")
problems_train, = glob_wildcards("./start_{problems, [a-z]_train}.f90")
problems_test, = glob_wildcards("./start_{problems, [a-z]_test}.f90")

# Rule all
rule all:
    input:
        expand("../Simulation_results/test/{problems}/done", problems = problems_train),
        expand("../Simulation_results/test/{problems}/done", problems = problems_test)

# Run initial problem set on naive networks
rule train:
    input:
        modules = modules,
        problem_files = expand("start_{problems}.f90", problems = problems_train)
    output:
        expand("../Simulation_results/test/{problems}/done", problems = problems_train)
    params:
        problem_name = problems_train
    shell:
        '''
        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs {input.modules} start.mod.f90 -o grns.e &&

        for problem in {params.problem_name}
        do

        echo start_$problem.f90

        cp start_$problem.f90 start.mod.f90

        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs {input.modules} start.mod.f90 -o grns.e &&

        ./grns.e &&

        mv ./files/GRN*.dat ../Simulation_results/test/$problem/ &&
        mv ./files/PHE*.dat ../Simulation_results/test/$problem/ &&

        rm grns.e &&

        touch ../Simulation_results/test/$problem/done

        done
        '''

# Run trained networks on all problems
rule test:
    input:
        modules = modules,
        problem_files = expand("start_{problems}.f90", problems = problems_test),
        grn_tokens = expand("../Simulation_results/test/{problems}/done", problems = problems_train)
    output:
        expand("../Simulation_results/test/{problems}/done", problems = problems_test)
    params:
        problem_name = problems_test
    shell:
        '''

        cp ../Simulation_results/test/*/GRN_* ./files &&
        ls files/GRN* | grep -o "GRN.*" > GRNfiles.txt &&

        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs {input.modules} start.mod.f90 -o grns.e &&

        for problem in {params.problem_name}
        do

        echo start_$problem.f90

        cp start_$problem.f90 start.mod.f90

        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs {input.modules} start.mod.f90 -o grns.e &&

        ./grns.e &&

        mv ./files/GRN_[0-9]*.dat ../Simulation_results/test/$problem/ &&
        mv ./files/PHE_[0-9]*.dat ../Simulation_results/test/$problem/ &&

        rm grns.e &&

        touch ../Simulation_results/test/$problem/done

        done

        '''
