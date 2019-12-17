
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
        expand("../Simulation_results/{problems}/done", problems = problems_test),
        expand("../Simulation_results/{problems}/bomb", problems = problems_train + problems_test)

# Run initial problem set on naive networks
rule train:
    input:

    output:
        "files/problems_trained"
    params:
        problem_train = problems_train,
        modules = modules
    shell:
        '''
        rm -f *.e

        # Compile executables for each problem
        for problem in {params.problem_train}
        do
        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs start_$problem.f90 {params.modules} -o $problem.e
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
        touch(expand("../Simulation_results/{problems}/done", problems = problems_train))
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

        '''

# Run test simulations initiated from all timepoints of the test set
rule test_all_setup:
    input:
        grn_tokens = expand("../Simulation_results/{problems}/done", problems = problems_train)
    output:
        "{problems, [a,b,n]}_test.e"
    params:
        modules = modules,
        problem_train = expand("{problems}_train", problems = problems_all_timepoints),
        problem_test = "{problems}_test"
    resources:
        GRNfile = 1
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
            gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs start_$problem.f90 {params.modules} -o $problem.e
        done
        '''

rule test_fin_setup:
    input:
        grn_tokens = expand("../Simulation_results/{problems}/done", problems = problems_train)
    output:
        "{problems, [d,e,f]}_test.e"
    params:
        modules = modules,
        problem_train = expand("{problems}_train", problems = problems_final_timepoints),
        problem_test = "{problems}_test"
    resources:
        GRNfile = 1
    shell:
        '''
        # Clean any eventual extra GRNs in seed directory
        rm -f files/GRN* &&

        # Copy GRNs to use as source for testing
        for grn in {params.problem_train}
        do
            cp -u ../Simulation_results/$grn/GRN*GRN*10.dat ./files
        done

        # Create list of GRN sources
        ls files/GRN* | grep -o "GRN.*" > GRNfiles.txt

        # Compile executables for each problem
        for problem in {params.problem_test}
        do
            gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs start_$problem.f90 {params.modules} -o $problem.e
        done
        '''

rule test:
    input:
        '{problem}_test.e'
    output:
        touch('files/done_{problem}_test')
    shell:
        '''
        ./{input}
        '''

rule test_sort:
    input:
        "files/done_{problem_test, [a-z]_test}"
    output:
        touch("../Simulation_results/{problem_test}/done")
    params:
        problem_codes = problem_codes,
        problem_names = problem_names

    shell:
        '''
        # Clean up binary files
        rm -f development.mod start.mod
        rm -f *.e
        rm -f files/GRN*

        # Transfer all results in respective folders
        parallel --jobs 6 --link \
        mv ./GRN_{{1}}*.dat \
        ../Simulation_results/{{2}}_test/ \
        ::: {params.problem_codes} \
        ::: {params.problem_names}

        parallel --jobs 6 --link \
        mv ./PHE_{{1}}*.dat \
        ../Simulation_results/{{2}}_test/ \
        ::: {params.problem_codes} \
        ::: {params.problem_names}

        '''

rule bomb:
    input:
        token = "../Simulation_results/{problem}/done",
        directory = "../Simulation_results/{problem}"
    output:
        touch("files/done_{problem}_bomb")
    resources:
        GRNfile = 1
    shell:
        '''
        # Grep all GRNs and move them into the files folder

        rm -f files/GRN_*
        for problem in {input.directory}
        do
            cp -u ../Simulation_results/$problem/GRN_* files
        done

        # Create list of GRN sources
        ls files/GRN* | grep -o "GRN.*" > GRNfiles.txt

        # Compule and run bomb script on all GRNs
        gfortran bomb.f90 -o bomb.e
        ./bomb.e
        '''

rule bomb_sort:
    input:
        "files/done_{problem}_bomb"
    output:
        directory("../Simulation_results/{problem}/bomb")
    params:
        problem_names = problem_names,
        problem_codes = problem_codes,
    shell:
        '''
        # Clean up binary files and source GRNs
        rm -f development.mod start.mod
        rm -f *.e
        rm -f files/GRN*

        # Transfer all results in respective folders
        problem_codes=({params.problem_codes})
        problem_names=({problem_names})

        for i in $(seq 0 5)
        do
        mkdir -p ../Simulation_results/${{problem_names[$i]}}_train/bomb

        find . -maxdepth 1 -name 'GRN*'${{problem_codes[$i]}}'*.dat' \
        -exec mv -t ../Simulation_results/${{problem_names[$i]}}_train/bomb {{}} \+
        done
        '''
