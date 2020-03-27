
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
        expand("../Simulation_results/{problems}", problems = problems_train),
        expand("../Simulation_results/{problems}", problems = problems_test),
        expand("../Simulation_results/bomb/{problems}", problems = problems_train + problems_test)

# Run initial problem set on naive networks
rule train:
    input:

    output:
        touch("files/problems_trained")
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

        '''

# Move results from training into dedicated results folders
rule train_sort:
    input:
        "files/problems_trained"
    output:
        directory(expand("../Simulation_results/{problems}", problems = problems_train))
    params:
        problem_names = problem_names,
        problem_codes = problem_codes,
        problems_train = problems_train
    shell:
        '''
        # Clean up binary files
        rm -f development.mod start.mod
        rm -f *.e

        # Create results folders
        for problem in {params.problem_names}
        do
        mkdir -p ../Simulation_results/${{problem}}_train
        done

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
        grn_tokens = expand("../Simulation_results/{problems}", problems = problems_train)
    output:
        "{problems, [a,b,n]}_test.e",
    params:
        modules = modules,
        problem_train = expand("{problems}_train", problems = problems_all_timepoints),
        problem_test = "{problems}_test"
    resources:
        GRNfile = 1
    shell:
        '''
        # remove old grnfile
        rm -f GRNfiles.txt

        # Copy GRNs to use as source for testing
        for grn in {params.problem_train}
        do
            cp -u ../Simulation_results/$grn/GRN* ./files
            ls ../Simulation_results/$grn/GRN*GRN* | grep -o "GRN.*" >> GRNfiles.txt
        done

        # Compile problem executable
        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs \
        start_{params.problem_test}.f90 {params.modules} \
        -o {params.problem_test}.e
        '''

rule test_fin_setup:
    input:
        grn_tokens = expand("../Simulation_results/{problems}", problems = problems_train)
    output:
        "{problems, [d,e,f]}_test.e",
    params:
        modules = modules,
        problem_train = expand("{problems}_train", problems = problems_final_timepoints + problems_final_timepoints),
        problem_test = "{problems}_test"
    resources:
        GRNfile = 1
    shell:
        '''
        # remove old grnfile
        rm -f GRNfiles.txt

        # Copy GRNs to use as source for testing and add to GRNfile
        for grn in {params.problem_train}
        do
            cp -u ../Simulation_results/$grn/GRN*GRN*T10.dat ./files
            ls ../Simulation_results/$grn/GRN*GRN*T10.dat | grep -o "GRN.*" >> GRNfiles.txt
        done

        # Compile problem executable
        gfortran -w -fexceptions -fno-underscoring -Wall -Wtabs \
        start_{params.problem_test}.f90 {params.modules} \
        -o {params.problem_test}.e
        '''

rule test:
    input:
        executable = '{problem}_test.e'
    output:
        touch('files/done_{problem}_test')
    shell:
        '''
        ./{input.executable}
        '''

rule test_sort:
    input:
        "files/done_{problem, [a-z]}_test"
    output:
        directory("../Simulation_results/{problem}_test")
    params:
    # This function matches the problem name (as set in the wildcard 'problem')
    # to its problem code (corresponding element in the tuple problem_codes)
        problem_code = lambda wildcards: problem_codes[problem_names.index(wildcards.problem[0])]
    shell:
        '''
        # Create target folder
        mkdir -p {output}

        find . -maxdepth 1 -regextype posix-egrep -regex '.*[0-9]GRN_'{params.problem_code}'.*' \
        -exec mv -t {output} {{}} \;

        find . -maxdepth 1 -regextype posix-egrep -regex '.*PHEN_TE_'{params.problem_code}'*.dat' \
        -exec mv -t {output} {{}} \;

        '''

rule bomb:
    input:
        "../Simulation_results/{problem, ([a-z]_test|[a-z]_train)}"
    output:
        touch("files/done_{problem, ([a-z]_train)|([a-z]_test)}_bomb")
    resources:
        GRNfile = 1
    shell:
        '''
        # Grep all GRNs and move them into the files folder
        rm -f files/GRN_*
        for problem in {wildcards.problem}
        do
            if grep -q 'train' ; then
                find ../Simulation_results/$problem \\
                -regextype posix-extended \\
                -regex '.*/GRN.*(0[0-9]|10).dat' \\
                -exec cp {{}} files \;
            elif grep -q 'test' ; then
                find ../Simulation_results/$problem \\
                -regextype posix-extended \\
                -regex '.*/GRN.*(0[0-9]|10).*GRN.*.dat' \\
                -exec cp {{}} files \;
            fi
        done

        # Create list of GRN sources (grep to remove base path)
        ls files/GRN*.dat | grep -o "GRN.*" > GRNfiles.txt

        # Compule and run bomb script on all GRNs
        gfortran bomb.f90 -o bomb.e
        ./bomb.e
        '''


rule bomb_sort:
# This rule uses the find -exec command combination to move files because we have too many files for the input to mv
# The use of a custom defined lambda function allows the use of the problem wildcard as an input to the param
    input:
        "files/done_{problem}_bomb"
    output:
        directory("../Simulation_results/bomb/{problem, ([a-z]_train)|([a-z]_test)}")
    params:
    # This function matches the problem name (as set in the wildcard 'problem')
    # to its problem code (corresponding element in the tuple problem_codes)
        problem_code = lambda wildcards: problem_codes[problem_names.index(wildcards.problem[0])]
    resources:
        GRNfile = 1
    shell:
        '''
        # Clean up binary files and source GRNs
        rm -f development.mod start.mod
        rm -f *.e
        rm -f files/GRN*

        # Create target folder
        mkdir -p {output}

        # Transfer results from the source problem to respective folder

        find . -maxdepth 1 -name 'GRN*'{params.problem_code}'*.dat' \
        -exec mv -t {output} {{}} \+
        '''
