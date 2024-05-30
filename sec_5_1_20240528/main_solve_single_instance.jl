import CSV, DataStructures, DelimitedFiles, Formatting, Glob, Gurobi, Logging, MathOptInterface
using DataFrames, JuMP

include("src/read_data.jl")
include("src/get_model.jl")
include("src/get_output.jl")

const MOI = MathOptInterface
const GRB_ENV = Gurobi.Env()


function open_global_logger(logfile)
    ##### Open a logger
    io = open(logfile, "w+")
    logger = Logging.SimpleLogger(io, Logging.Info)
    Logging.global_logger(logger)
    return io
end

function solve_single_instance(result_dir, dataset_dir, num_trees, formulation, is_LP, time_limit_sec)
    io = open_global_logger(result_dir * "log.log")

    data = read_data(dataset_dir, Array(1:num_trees), dir_to_save_data=result_dir)
 
    optimizer_attributes = Dict("LogFile" => "$(result_dir)gurobi.log", "LogToConsole" => 0, "TimeLimit" => time_limit_sec, "MIPGap" => 0.0)
    lazylevel = 2
    model, data = get_model(data, optimizer_attributes, formulation, lazylevel, is_LP)

    optimize!(model)

    output, data = get_output(model, data)
    merge!(output, Dict(:dataset=>dataset_dir, :num_trees=>num_trees, :formulation=>formulation, :time_limit_sec=>time_limit_sec))
    CSV.write(result_dir * "result.csv", output)
    CSV.write(result_dir * "df_sol.csv", data["df_sol"])

    println("formulation = $(output[:formulation]), objval = $(output[:objval])")

    close(io)
end

function main()
    result_dir = "result/"
    dataset_dir = "dataset/concrete_bt/" # "concrete_bt", "concrete_rf", "redwine_bt", "redwine_rf"
    num_trees = 10
    is_LP = false
    time_limit_sec = 300
    for formulation in ["TEOM", "TEOR", "TEOC"]
        solve_single_instance(result_dir * formulation * "/", dataset_dir, num_trees, formulation, is_LP, time_limit_sec)
    end
end

main()

