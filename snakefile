
# Define input problems and modules
modules = ("development.mod.f90", "grns2.f90")
problems_train, = glob_wildcards("./start_{problems, [a-z]_train}.f90")
problems_test, = glob_wildcards("./start_{problems, [a-z]_test}.f90")

problems_all_timepoints = ("a", "b", "n")
# Rule all
rule all:
    input:
        expand("../Simulation_results/{problems}/done", problems = problems_train),
        expand("../Simulation_results/{problems}_test/done", problems = problems_all_timepoints)

# Run initial problem set on naive networks
rule train:
    input:
        modules = modules,
        problem_files = expand("start_{problems}.f90", problems = problems_train)
    output:
        expand("../Simulation_results/{problems}/done", problems = problems_train)
    params:
        problem_name = problems_train
    shell:
        '''
        for problem in {params.problem_name}
        do

        echo start_$problem.f90

        cp start_$problem.f90 start.mod.f90 &&

        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs start.mod.f90 {input.modules} -o grns.e &&

        ./grns.e &&

        mv ./files/GRN*.dat ../Simulation_results/$problem/ &&
        mv ./files/PHE*.dat ../Simulation_results/$problem/ &&

        rm grns.e start.mod.f90 development.mod start.mod &&

        touch ../Simulation_results/$problem/done

        done
        '''

# Run test simulations initiated from all timepoints of the test set
rule test_all_timepoints:
    input:
        modules = modules,
        problem_files = expand("start_{problems}_test.f90", problems = problems_all_timepoints),
        grn_tokens = expand("../Simulation_results/{problems}_train/done", problems = problems_all_timepoints)
    output:
        expand("../Simulation_results/{problems}_test/done", problems = problems_all_timepoints)
    params:
        problem_name = expand("{problems}_test", problems = problems_all_timepoints)
    shell:
        '''

        rm -f files/GRN* &&

        for grn in {params.problem_name}
        do

        cp -u ../Simulation_results/$grn/GRN_* ./files &&

        done

        ls files/GRN* | grep -o "GRN.*" > GRNfiles.txt

        for problem in {params.problem_name}
        do

        echo start_$problem.f90

        cp start_$problem.f90 start.mod.f90 &&

        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs start.mod.f90 {input.modules} -o grns.e &&

        ./grns.e &&

        mv ./files/GRN_[0-9]*.dat ../Simulation_results/$problem/ &&
        mv ./files/PHE_[0-9]*.dat ../Simulation_results/$problem/ &&
        mv ./GRNstatus.txt ../Simulation_results/$problem/ &&

        rm grns.e start.mod.f90 development.mod start.mod &&

        touch ../Simulation_results/$problem/done

        done

        rm ./files/GRN*

        '''
