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

function solve_single_instance(dataset_dir, num_trees, formulation, is_LP, time_limit_sec)
    io = open_global_logger("log.log")

    # Set Parameters
    optimizer_attributes = Dict("LogFile" => "gurobi.log", "LogToConsole" => 0, "TimeLimit" => time_limit_sec, "MIPGap" => 0.0)
    lazylevel = 2

    # Solve LP
    data = read_data(dataset_dir, Array(1:num_trees))
    model, data = get_model(data, optimizer_attributes, formulation, lazylevel, true)
    optimize!(model)
    output, data = get_output(model, data)
    objval_LP = output[:objval]

    # Solve MIP
    data = read_data(dataset_dir, Array(1:num_trees))
    model, data = get_model(data, optimizer_attributes, formulation, lazylevel, false)
    optimize!(model)
    output, data = get_output(model, data)

    # Save result
    merge!(output, Dict(:dataset=>dataset_dir, :num_trees=>num_trees, :formulation=>formulation, :time_limit_sec=>time_limit_sec, :objval_LP=>objval_LP))
    result_file = "result.csv"
    CSV.write(result_file, DataFrame(output), append=isfile(result_file))

    close(io)
end

function main()
    dataset_dir = "dataset/concrete_bt/" # "concrete_bt", "concrete_rf", "redwine_bt", "redwine_rf"
    num_trees = 10
    is_LP = false
    time_limit_sec = 300
    for formulation in ["TEOM", "TEOR", "TEOC"]
        solve_single_instance(dataset_dir, num_trees, formulation, is_LP, time_limit_sec)
    end
end

main()

