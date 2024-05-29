using CSV, DataFrames, Dates, DelimitedFiles, Formatting, Glob, Logging, Missings, Printf, Random

include("read_data.jl")
include("solve_mip.jl")

function ensembleopt(dir::String, treedir::String, treeids::Array{Int}, method::String, timelimit::Number; 
                    loglevel="Info", df_sol_fixed=DataFrame(), preprocessing_level=0, boxes=[], z_vartype='C', df_forest_folder="")
    ##### Create a directory to save the progress and the result
    try
        mkdir(dir)
    catch e
        nothing
    end
    
    ##### Open a logger
    io = open_global_logger(dir * "log.log", loglevel=loglevel)

    ##### Check the input
    @info @sprintf("ensembleopt() %s\t%s\t%s\t%s", dir, treedir, length(treeids), method)
    @info "boxes" boxes

    ##### Read data
    preprocessing = []
    if preprocessing_level >= 2 && occursin("assortment", treedir)
        preprocessing = ["remove_nodes_with_empty_box", "revalue_splitval"] 
    elseif preprocessing_level >= 1
        preprocessing = ["remove_nodes_with_empty_box"]
    end

    if df_forest_folder == ""
        data = Dict()
        if treedir == "dataset/redwinepair/"
            data = read_data(["dataset/RedwinePair1_Q_GBT_T5000/", "dataset/RedwinePair1_A_GBT_T5000/"], 
                            [Array(1:Int(length(treeids) / 2)) for i = 1:2]; weights=[2.0, -1.0], preprocessing=preprocessing)
        else
            data = read_data([treedir], [treeids]; preprocessing=preprocessing)
        end
    else
        data = read_data_from_saved(df_forest_folder)
    end

    ##### Partition leaves
    for box in boxes
        # println()
        # println(box)
        data = partition_leaves(data, box)
        data = remove_leaves(data, box)
    end
    boxes = []

    ##### Save data
    # println("save_data")
    save_data(dir, data)


    ##### Get optbox if needed
    optbox = DataFrame()
    if method == "liftall-optbox"
        optdir = replace(dir, "liftall-optbox" => "liftall")
        if isdir(optdir) && isfile(dir * "result.csv")
            optbox = get_optbox_from_result_directory(optdir)
            @info optbox
        else
            return Dict()
        end
    end

    ##### Solve a MIP model
    output = solve_mip(dir, data, method, timelimit, df_sol_fixed=df_sol_fixed, optbox=optbox, boxes=boxes, z_vartype=z_vartype)
    @info join(["$index: $value" for (index, value) in output], "\n")

    ##### Close the logger
    close(io)

    output
end

function check_result(resultfile, dataset, num_trees, seed, method)
    df = DataFrame(CSV.File(resultfile))
    
    for r in eachrow(df)
        if r.dataset == dataset && r.num_trees == num_trees && r.method == method
            if seed == 0 && !("seed" in names(df))
                return true
            elseif seed == 0 && r.seed in [missing, 0]
                return true
            elseif seed > 0 && !("seed" in names(df))
                return false
            elseif seed > 0 && !ismissing(r.seed) && r.seed == seed
                return true
            end
        end
    end

    return false
end
