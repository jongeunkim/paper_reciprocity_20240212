# using JuMP 
# import Plots, Missings,
import Glob

include("src/read_data.jl")
include("src/solve_mip.jl")
# include("src/utils_JuMP.jl")
# include("src/post_analysis.jl")
include("src/save_result.jl")
# include("src/partition_leaves.jl")
# include("src/splitconstr")


function solve_instance(dir, datadir, num_trees, method)
    if isdir(dir) 
        if isfile(dir * "result.csv")
            return
        else
            rm(dir, force=true, recursive=true)
            mkdir(dir)
        end    
    else
        rm(dir, force=true, recursive=true)
        mkdir(dir)
    end
    io = open_global_logger(dir * "log.log")

    preprocessing=["remove_nodes_with_empty_box", "revalue_splitval"]
    if occursin("redwinepair", datadir)
        data = read_data(["dataset/RedwinePair1_Q_GBT_T5000/", "dataset/RedwinePair1_A_GBT_T5000/"], 
            [Array(1:Int(num_trees / 2)) for i = 1:2]; weights=[2.0, -1.0], preprocessing=preprocessing)
    elseif occursin("redwinefeas", datadir)
        data = read_data(["dataset/redwine_T5000/", "dataset/redwinefeas_T5000_Depth9/"], 
            [Array(1:Int(num_trees / 2)) for i = 1:2]; weights=[1.0, 1.0], preprocessing=preprocessing)
    else
        data = read_data([datadir], [Array(1:num_trees)], preprocessing=preprocessing)
    end
    save_data(dir, data)

    formulation = split(method, "-")[1]
    optimizer_attributes = Dict("LogFile" => "$(dir)gurobi.log", "LogToConsole" => 0, "TimeLimit" => 3600, "MIPGap" => 0.0)

    objval_boxes = nothing
    if occursin("-fixed", method)
        model, data = get_model(data, optimizer_attributes, formulation, 0)
    elseif occursin("-lazy3", method)
        model, data = get_model(data, optimizer_attributes, formulation, 3)
    else
        model, data = get_model(data, optimizer_attributes, formulation, 2)
    end
    if occursin("-LP", method)
        model = relax_integrality(model)
        model = covert_LP_to_MIP(model)
    end

    optimize!(model)

    output, data = get_model_output(model, data)
    merge!(output, Dict(:dataset=>datadir, :num_trees=>num_trees, :seed=>0, :method=>method, :timelimit=>3600, :preprocessing_level=>1))
    save_result(output, dir * "result.csv")

    CSV.write(dir * "df_sol.csv", data["df_sol"])

    close(io)
end

function main()
    dir_global = "result_compare_formulations_20210121_2/"
    if !isdir(dir_global)
        mkdir(dir_global)
    end

    Ts = 20:20:500
    Ms = ["misic", "misic-LP", "liftall", "liftall-LP"]
    # Ms = ["misic", "misic-fixed", "misic-LP", "liftx", "liftx-fixed", "liftx-LP", "liftall", "liftall-fixed", "liftall-LP"]
    Ds = []
    Ds = [String(split(dir, '/')[end-1]) for dir in Glob.glob("dataset/**/")]
    filter!(r -> occursin("concrete", r), Ds)
    # push!(Ds, "redwinefeas")
    # Ds = [datadir for datadir in Ds if occursin("assortment", datadir) || occursin("concrete", datadir) || occursin("redwine", datadir) || occursin("maxcut_N50", datadir)]

    println("Ts = $Ts")
    println("Ms = $Ms")
    println("Ds = $Ds")

    for num_trees in Ts, data in Ds, method in Ms
        dir_instance = "$(dir_global)$(data)-T$(num_trees)/"
        if !isdir(dir_instance)
            mkdir(dir_instance)
        end

        dir = "$(dir_instance)$(method)/"
        datadir = "dataset/$(data)/"
        println("$num_trees\t$data\t$method")
        solve_instance(dir, datadir, num_trees, method)
    end

    println("main() - all done")
end

main()

println("all done")