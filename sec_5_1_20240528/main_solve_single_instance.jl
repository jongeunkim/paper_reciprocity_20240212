import Glob
using JuMP

include("src/read_data.jl")
include("src/solve_mip.jl")
include("src/save_result.jl")
include("src/get_output.jl")

function solve_single_instance(result_dir, dataset_dir, num_trees, method, time_limit_sec)
    io = open_global_logger(result_dir * "log.log")

    preprocessing=[]
    data = read_data([dataset_dir], [Array(1:num_trees)], preprocessing=preprocessing)
    save_data(result_dir, data)

    formulation = split(method, "-")[1]
    optimizer_attributes = Dict("LogFile" => "$(result_dir)gurobi.log", "LogToConsole" => 0, "TimeLimit" => time_limit_sec, "MIPGap" => 0.0)
    lazylevel = 2
    model, data = get_model(data, optimizer_attributes, formulation, lazylevel)
    if occursin("-LP", method)
        model = relax_integrality(model)
        model = covert_LP_to_MIP(model)
    end

    optimize!(model)

    output, data = get_output(model, data)
    merge!(output, Dict(:dataset=>dataset_dir, :num_trees=>num_trees, :method=>method, :time_limit_sec=>time_limit_sec))
    CSV.write(result_dir * "result.csv", output)
    CSV.write(result_dir * "df_sol.csv", data["df_sol"])

    close(io)
end

function main()
    result_dir = "result/"
    dataset_dir = "dataset/concrete_bt/"
    num_trees = 10
    method = "misic"    # Ms = ["misic", "misic-fixed", "misic-LP", "liftx", "liftx-fixed", "liftx-LP", "liftall", "liftall-fixed", "liftall-LP"]
    time_limit_sec = 300
    for method in ["misic", "liftx", "liftall"]
        solve_single_instance(result_dir, dataset_dir, num_trees, method, time_limit_sec)
    end
end

main()

