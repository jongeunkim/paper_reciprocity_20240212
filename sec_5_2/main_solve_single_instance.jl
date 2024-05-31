include("src/mctp_experiment.jl")

function main()
    ##### Parameters
    M, N = 5, 2
    C, D, T = 2, 5, 1
    random_seed = 1 # 1 to 10
    formulation = "TEOM" # "TEOM", "TEOR", "TEOC", "DCC", "DLog", "MC"
    result_file = "result.csv"
    timelimit_sec = 300

    for formulation in ["TEOM", "TEOR", "TEOC", "DCC", "DLog", "MC"]
        ##### Solve
        output = run(formulation, random_seed, M, N, C, D, T; timelimit=timelimit_sec, verbose=0)

        ##### Write
        instance_info = Dict("M"=>M, "N"=>N, "C"=>C, "D"=>D, "T"=>T, "seed"=>random_seed, "formulation"=>formulation, "timelimit"=>timelimit_sec)
        merge!(output, instance_info)
        CSV.write(result_file, DataFrames.DataFrame(output), append=isfile(result_file))
    end
end

main()
