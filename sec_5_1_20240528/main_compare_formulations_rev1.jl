import Dates, Glob

include("src/read_data.jl")
# include("src/solve_mip.jl")
include("src/solve_mip_scip.jl")
include("src/save_result.jl")
# include("src/set_optimal_start_values.jl")
include("src/get_model_output.jl")

function is_instance_completed(dir)
    return isdir(dir) && isfile(dir * "result.csv")
end

function solve_instance(dir, datadir, num_trees, method)
    if is_instance_completed(dir)
        println("$dir is already done.")
        return
    end

    rm(dir, force=true, recursive=true)
    mkdir(dir)
    io = open_global_logger(dir * "log.log")

    # preprocessing=["remove_nodes_with_empty_box", "revalue_splitval"]
    preprocessing=[]
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
    optimizer_attributes = Dict("LogFile" => "$(dir)gurobi.log", "LogToConsole" => 0, "TimeLimit" => 3600 * 3, "MIPGap" => 0.0)

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
        if !occursin("-fixed", method)
            model = covert_LP_to_MIP(model)
        end
    end

    @info "First solve"
    optimize!(model)

    if occursin("-givenoptsol", method)
        @info "Solve with the optimal solution"
        x = JuMP.all_variables(model)
        x_solution = Int.(round.(value.(x)))
        JuMP.set_start_value.(x, x_solution)
        optimize!(model)
    end

    output, data = get_model_output(model, data)
    merge!(output, Dict(:dataset=>datadir, :num_trees=>num_trees, :seed=>0, :method=>method, :timelimit=>3600, :preprocessing_level=>1))
    save_result(output, dir * "result.csv")

    CSV.write(dir * "df_sol.csv", data["df_sol"])

    close(io)
end

function GenerateDirectoriesAndReturnFinalPath(dir_name_list)
    path = ""
    for dir_name in dir_name_list
        path *= dir_name * "/"
        if !isdir(path)
            mkdir(path)
        end
    end
    return path
end

function main(args)
    TODAY = Dates.format(Dates.today(), "yyyymmdd")
    DATASET_NUMTREES = Dict(
        "concrete_T1000_depth10" => [100, 300, 500],   
        "concrete_T5000_GBT" => [100, 300, 500],   
        "redwine_T1000_depth10" => [100, 200, 300, 500],   
        "redwine_T5000_GBT" => [100, 200, 300, 500],   
    )
    # METHODS = ["minimal-fixed", "minimal-fixed-LP", "misic-fixed", "misic-fixed-LP", "liftall-fixed", "liftall-fixed-LP", "liftx-fixed", "liftx-fixed-LP"]
    # METHODS = ["minimal-fixed", "misic-fixed", "liftall-fixed", "liftx-fixed"]
    # METHODS = ["minimal-LP", "misic-LP", "liftx-LP", "minimal", "misic", "liftx", "minimal-givenoptsol", "misic-givenoptsol", "liftx-givenoptsol"]
    # METHODS = ["misic", "misic-fixed", "misic-LP", "liftx", "liftx-fixed", "liftx-LP", "liftall", "liftall-fixed", "liftall-LP"]

    METHODS = ["minimal-LP", "misic-LP", "liftx-LP", "minimal", "misic", "liftx", "liftall-LP", "liftall"]

    DIR_GLOBAL = ""
    if length(args) > 0
        DIR_GLOBAL = args[1]
    else
        DIR_GLOBAL = "result_compare_formulations_$TODAY"
        
        cnt = 0
        while isdir("$(DIR_GLOBAL)_$cnt")
            cnt += 1
        end
        mkdir("$(DIR_GLOBAL)_$cnt")
        DIR_GLOBAL = "$(DIR_GLOBAL)_$cnt"
    end
    println("DIR_GLOBAL = ", DIR_GLOBAL)

    instances = []
    for (data, num_trees_list) in DATASET_NUMTREES, method in METHODS
        for num_trees in num_trees_list
            push!(instances, (num_trees, data, method))
        end
    end
    sort!(instances)
    # instances = [(20, "concrete_T1000_depth10", "misic-fixed")]

    num_all_instances = length(instances)
    instances = filter(ins -> !is_instance_completed(GenerateDirectoriesAndReturnFinalPath([DIR_GLOBAL, ins[2], "$(ins[1])", ins[3]])), instances)
    num_instances_to_solve = length(instances)
    println("$num_instances_to_solve / $num_all_instances")
    for (num_trees, data, method) in instances
        dir_instance = GenerateDirectoriesAndReturnFinalPath([DIR_GLOBAL, data, "$num_trees", method])
        println("Solve $dir_instance")
        solve_instance(dir_instance, "dataset/$(data)/", num_trees, method)
    end

    println("main() - all done")
end

main(ARGS)

println("all done")