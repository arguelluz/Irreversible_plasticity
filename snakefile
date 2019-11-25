
# Define input problems and modules
modules = ("development.mod.f90", "grns2.f90")
problems_train, = glob_wildcards("./start_{problems, [a-z]_train}.f90")
problems_test, = glob_wildcards("./start_{problems, [a-z]_test}.f90")

problems_all_timepoints = ("a", "b", "n")
problems_final_timepoints = ("d", "e", "f")

problem_names = ("n", "a", "b", "d", "e", "f")
problem_codes = ("122345", "312452", "254213", "223344", "443322", "224411")

# Rule all
rule all:
    input:
        expand("../Simulation_results/{problems}/done", problems = problems_train),
        expand("../Simulation_results/{problems}/bomb/done", problems = problems_train),
        expand("../Simulation_results/{problems}_test/done", problems = problems_all_timepoints),
        expand("../Simulation_results/{problems}_test/done", problems = problems_final_timepoints)

# Run initial problem set on naive networks
rule train:
    input:
        modules = modules,
        problem_files = expand("start_{problems}.f90", problems = problems_train)
    output:
        "files/problems_trained"
    params:
        problem_train = problems_train
    shell:
        '''
        rm -f *.e

        # Compile executables for each problem
        for problem in {params.problem_train}
        do
        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs start_$problem.f90 {input.modules} -o $problem.e
        done

        # Run all problems in parallel
        parallel --nice 19 \
            ./{{}}.e \
        ::: {params.problem_train}

        touch files/problems_trained

        '''

# Move results from training into dedicated results folders
rule train_sort:
    input:
        "files/problems_trained"
    output:
        expand("../Simulation_results/{problems}/done", problems = problems_train)
    params:
        problem_names = problem_names,
        problem_codes = problem_codes,
        problems_train = problems_train
    shell:
        '''
        # Clean up binary files
        rm -f development.mod start.mod
        rm -f *.e

        # Transfer all results in respective folders
        parallel --jobs 3 --link \
        mv ./GRN*{{1}}*.dat \
        ../Simulation_results/{{2}}_train/ \
        ::: {params.problem_codes} \
        ::: {params.problem_names}

        parallel --jobs 3 --link \
        mv ./PHE*{{1}}*.dat \
        ../Simulation_results/{{2}}_train/ \
        ::: {params.problem_codes} \
        ::: {params.problem_names}

        for problem in {params.problems_train}
        do
        touch   ../Simulation_results/$problem/done
        done
        '''

# Run test simulations initiated from all timepoints of the test set
rule test_all_timepoints:
    input:
        modules = modules,
        problem_files = expand("start_{problems}_test.f90", problems = problems_all_timepoints),
        grn_tokens = expand("../Simulation_results/{problems}_train/done", problems = problems_all_timepoints)
    output:
        "files/problems_tested_all_timepoints"
    params:
        problem_train = expand("{problems}_train", problems = problems_all_timepoints),
        problem_test = expand("{problems}_test", problems = problems_all_timepoints)
    shell:
        '''
        # Clean any eventual extra GRNs in seed directory
        rm -f files/GRN* &&

        # Copy GRNs to use as source for testing
        for grn in {params.problem_train}
        do
            cp -u ../Simulation_results/$grn/GRN*GRN*[1-9].dat ./files
            cp -u ../Simulation_results/$grn/GRN*GRN*10.dat ./files
        done

        # Create list of GRN sources
        ls files/GRN* | grep -o "GRN.*" > GRNfiles.txt

        # Compile executables for each problem
        for problem in {params.problem_test}
        do
            gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs start_$problem.f90 {input.modules} -o $problem.e
        done

        # Run all testing simulations in parallel
        parallel \
            ./{{}}.e \
        ::: {params.problem_test}

        touch files/problems_tested_all_timepoints
        '''

rule test_final_timepoints:
    input:
        modules = modules,
        problem_files = expand("start_{problems}_test.f90", problems = problems_final_timepoints),
        grn_tokens = expand("../Simulation_results/{problems}_train/done", problems = problems_final_timepoints)
    output:
        "files/problems_tested_fin_timepoints"
    params:
        problem_train = expand("{problems}_train", problems = problems_final_timepoints),
        problem_test = expand("{problems}_test", problems = problems_final_timepoints)
    shell:
        '''
        # Clean any eventual extra GRNs in seed directory
        rm -f files/GRN* &&

        # Copy GRNs to use as source for testing
        for grn in {params.problem_train}
        do
            cp -u ../Simulation_results/$grn/GRN_*_T10.dat ./files
        done

        # Create list of GRN sources
        ls files/GRN* | grep -o "GRN.*" > GRNfiles.txt

        # Compile executables for each problem
        for problem in {params.problem_test}
        do
            gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs start_$problem.f90 {input.modules} -o $problem.e
        done

        # Run all testing simulations in parallel
        parallel \
            ./{{}}.e \
        ::: {params.problem_test}

        touch files/problems_tested_fin_timepoints
        '''

rule test_sort:
    input:
        "files/problems_tested_all_timepoints",
        "files/problems_tested_fin_timepoints"
    output:
        expand("../Simulation_results/{problems}/done", problems = problems_test)
    params:
        problems_test = problems_test,
        problem_codes = problem_codes,
        problem_names = problem_names

    shell:
        '''
        # Clean up binary files
        rm -f development.mod start.mod
        rm -f *.e
        rm files/GRN*


        # Transfer all results in respective folders
        parallel --jobs 3 --link \
        mv ./GRN_{{1}}*.dat \
        ../Simulation_results/{{2}}_test/ \
        ::: {params.problem_codes} \
        ::: {params.problem_names}

        parallel --jobs 3 --link \
        mv ./PHE_{{1}}*.dat \
        ../Simulation_results/{{2}}_test/ \
        ::: {params.problem_codes} \
        ::: {params.problem_names}

        for problem in {params.problems_test}
        do
        touch   ../Simulation_results/$problem/done
        done
        '''

rule train_bomb:
    input:
        expand("../Simulation_results/{problems}/done", problems = problems_train)
    output:
        "files/done_train_bomb"
    params:
        problem_name = problems_train
    shell:
        '''
        # Grep all GRNs and move them into the files folder

        rm -f files/GRN_*
        for problem in {params.problem_name}
        do
            cp -u ../Simulation_results/$problem/GRN_* files
        done

        # Create list of GRN sources
        ls files/GRN* | grep -o "GRN.*" > GRNfiles.txt

        # Compule and run bomb script on all GRNs
        gfortran bomb.f90 -o bomb.e
        ./bomb.e &&
        touch files/done_train_bomb
        '''

rule train_bomb_sort:
    input:
        "files/done_train_bomb"
    output:
        expand("../Simulation_results/{problems}/bomb/done", problems = problems_train)
    params:
        problem_names = problem_names,
        problem_codes = problem_codes,
        problems_train = problems_train
    shell:
        '''
        # Clean up binary files and source GRNs
        rm -f development.mod start.mod
        rm -f *.e
        rm files/GRN*

        # Transfer all results in respective folders
        parallel --jobs 3 --link \
        mv ./GRN*{{1}}*.dat \
        ../Simulation_results/{{2}}_train/bomb/ \
        ::: {params.problem_codes} \
        ::: {params.problem_names}

        for problem in {params.problems_train}
        do
        touch   ../Simulation_results/$problem/bomb/done
        done
        '''
