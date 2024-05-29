using DataFrames, CSV, Random, Glob, Formatting, Combinatorics

# include("utils.jl")
# include("callbacks.jl")
include("ensembleopt2.jl")
include("save_result.jl")
include("partition_leaves.jl")
# include("callbacks.jl")
include("post_analysis.jl")
# include("extended_optbox.jl")

# dir = "result_1201/"
# foreach(rm, Glob.glob(dir*"*"))

# dataset = "assortment_T5000_depth10"
# method = "liftall"
# preprocessing = 0
# dir_IP = "result_1201/$(dataset)_$(num_trees)_0_$(method)_3600_$preprocessing/"
# optbox = get_optbox_from_result_directory(dir_IP)
# ensembleopt(dir, "dataset/$(dataset)/", Array(1:num_trees), "liftall-optbox", 1200; loglevel="Info", preprocessing_level=preprocessing, optbox=optbox)

# dataset = "Redwine_T5000_GBT"
# num_trees = 100
# seed = 0
# method = "liftall-optbox"
# timelimit = 3600
# preprocessing = 1

# dir = "result_1201/$(dataset)_$(num_trees)_$(seed)_$(method)_$(timelimit)_$(preprocessing)/"
# args = [dir, "dataset/$(dataset)/", "$num_trees", "$seed", "$method", "$timelimit", "$preprocessing"]
# main(args)



function get_LPbox(dir, level)
    println("get_LPbox")
    ##### (1) Get boxes: Get optimal unitbox in splitvar1, splitval1
    intervals = Dict()
    df_sol = DataFrame(CSV.File(dir * "df_sol.csv"))
    df_sol_splits = filter(row -> row.treeid == 0 && row.splitvar1 > 0, df_sol)
    for df in groupby(df_sol_splits, ["splitvar1"])
        intervals[df.splitvar1[1]] = []

        gaps = []
        push!(gaps, df.sol[1])
        for i = 1:nrow(df) - 1
            push!(gaps, df.sol[i + 1] - df.sol[i])
        end
        push!(gaps, 1 - df.sol[end])
        sorted_indexes = sortperm(gaps, rev=true)
        # println(gaps)

        i = 1
        nextgap = 1
        while i <= level && nextgap >= 0.1 
            k = sorted_indexes[i]
            if k == 1
                push!(intervals[df.splitvar1[1]], [0, df.splitval1[1]])
            elseif k == length(sorted_indexes)
                push!(intervals[df.splitvar1[1]], [df.splitval1[end], typemax(Int)])
            else
                push!(intervals[df.splitvar1[1]], [df.splitval1[k - 1], df.splitval1[k]])
            end

            i += 1
            if i <= length(sorted_indexes)
                nextgap = gaps[sorted_indexes[i]]
            else
                nextgap = 0
            end
        end

        # gap_best = df.sol[1]
        # lb = 0
        # ub = df.splitval1[1]

        # for i=1:nrow(df)-1
        #     gap_current = df.sol[i+1] - df.sol[i]
        #     if gap_current > gap_best
        #         gap_best = gap_current
        #         lb = df.splitval1[i]
        #         ub = df.splitval1[i+1]
        #     end
        # end

        # if 1 - df.sol[end] > gap_best
        #     gap_best = 1 - df.sol[end] 
        #     lb = df.splitval1[end]
        #     ub = typemax(Int)
        # end

        # box[df.splitvar1[1]] = [lb,ub]
    end
    # println(intervals)

    ##### Cartesian product
    @info "num boxes: " prod([length(v) for (k, v) in intervals]) 

    boxes = [Dict()]
    for k in keys(intervals)
        boxes_new = []
        for (box, interval) in Iterators.product(boxes, intervals[k])
            box[k] = interval
            push!(boxes_new, deepcopy(box))
        end
        boxes = boxes_new
    end

    @info "num boxes: " length(boxes)
    # return 
    # for box in boxes
    #     println(box)
    # end
    # return

    ##### (2) Get boxes using pair of leaves
    # For every pair,
    # compute intersection
    # boxes = []
    

    ### splitvar1  to splitvar
    function splitvar1_to_splitvar(box)
        df_splits = DataFrame(CSV.File(dir * "df_splits.csv"))
        box_new = Dict()
        for (k, v) in box
            df = filter(row -> row.splitvar1 == k, df_splits)
            var = df.splitvar[1]
            
            if v[1] == 0
                lb = -Inf
            else
                df_lb = filter(row -> row.splitval1 == v[1], df)
                lb = df_lb.splitval[1]
            end

            if v[2] == typemax(Int)
                ub = Inf
            else
                df_ub = filter(row -> row.splitval1 == v[2], df)
                ub = df_ub.splitval[1]
            end
            box_new[k] = [lb,ub]
        end
        box_new
    end
    boxes = [splitvar1_to_splitvar(box) for box in boxes]
    println("splitvar1_to_splitvar() done")

    function find_intersection_of_leaves(box)
        df_forest = DataFrame(CSV.File(dir * "df_forest.csv"))
        ptr_forest = get_pointers(df_forest, "treeid")
        num_trees = length(ptr_forest) - 1
        box_new = Dict()
        for t = 1:num_trees
            df_tree = df_forest[ptr_forest[t]:ptr_forest[t + 1] - 1, :]
            leaf = find_leaf_containing_box(df_tree, box)
            if leaf > 0
                box_leaf = get_rectangle_in_tree(df_tree, leaf)
                box_new = intersect_rectangles(box_new, box_leaf)
            end
        end
        box_new
    end
    boxes = [find_intersection_of_leaves(box) for box in boxes]
    # println(box)
    @info "num boxes: " length(boxes) 
    println(format("num boxes: {}", length(boxes))) 
    boxes
end
# dir = "result_1209/redwine_T1000_Depth10_10_1_2/LP1/"
# get_LPbox(dir, 2)

function test(dir, datasetdir, num_trees, preprocessing_level, strategy, num_boxes)
    try
        mkdir(dir)
    catch
        nothing
    end

    treeids = Array(1:num_trees)
    timelimit = 3600
    method = "liftx"
    global_objval = -Inf
    global_break = false

    ### Solve IP
    r = 0
    while r <= 2
        r += 1

        ### Get boxes
        boxes = []
        if r > 1 && strategy in ["LPsol-B", "LPsol-C", "LPsol-gap-B", "LPsol-gap-C"]
            # println("Compute boxes - $strategy")
            for k = 1:(r - 1)
                df_forest = DataFrame(CSV.File(dir * "LP$k-$method/df_forest.csv"))
                df_sol = DataFrame(CSV.File(dir * "LP$k-$method/df_sol.csv"))
                df_forest.ysol = [sum(df_sol[j, "sol"] for j in df_forest[i,"leafbegin"]:df_forest[i,"leafend"]) for i = 1:nrow(df_forest)]
                df_forest.box_varindex = get_boxes_df_forest(df_forest, "splitvar1", "varindex", 0, typemax(Int))
                df_forest.nlsol = [compute_nonlinear_term(df_forest.box_varindex[i], df_sol) for i = 1:nrow(df_forest)]
                df_forest.ynlgap = abs.(df_forest.ysol .- df_forest.nlsol)
                df_forest.box = get_boxes_df_forest(df_forest, "splitvar", "splitval", -Inf, Inf)
                
                ### Sort depth 2 boxes by ynlgap
                df_filtered = DataFrame()
                if strategy in ["LPsol-gap-B", "LPsol-gap-C"]
                    df_filtered = sort(df_forest, ["ynlgap", "depth"], rev=[true, false])
                else 
                    min_gap = 0.01
                    df_filtered = sort(filter(row -> row.ynlgap >= min_gap, df_forest), ["depth", "ynlgap"], rev=[false, true])
                    while nrow(df_filtered) < num_boxes && min_gap >= 0.01
                        min_gap /= 2.0
                        df_filtered = sort(filter(row -> row.ynlgap >= min_gap, df_forest), ["depth", "ynlgap"], rev=[false, true])
                    end
                end
                CSV.write(dir * "df_filtered_$(k)-$(method).csv", df_filtered)

                i = 1
                while i <= nrow(df_filtered) && length(boxes) < k * num_boxes
                    box_candidate = deepcopy(df_filtered.box[i])
                    # if all([isempty_box(intersect_boxes([box, box_candidate])) for box in boxes])
                    push!(boxes, box_candidate)
                    boxes = remove_duplicate_boxes(boxes)
                    # end                    
                    i += 1
                end

                if nrow(df_filtered) == 0 || length(boxes) > 12
                    global_break = true
                end
                # println(length(boxes), ", ", global_break)
            end
            # println(boxes)
            if length(boxes) < (r-1) * num_boxes
                global_break = true
            end

        elseif r == 2 && strategy in ["optsol-B", "optsol-C"]
            df_forest = DataFrame(CSV.File(dir * "IP1-$method/df_forest.csv"))
            df_sol = DataFrame(CSV.File(dir * "IP1-$method/df_sol.csv"))
            df_forest.ysol = [sum(df_sol[j, "sol"] for j in df_forest[i,"leafbegin"]:df_forest[i,"leafend"]) for i = 1:nrow(df_forest)]
            df_forest.box = get_boxes_df_forest(df_forest, "splitvar", "splitval", -Inf, Inf)
            optimal_boxes = [df_forest.box[j] for j = 1:nrow(df_forest) if df_forest.lchild[j] == 0 && df_forest.ysol[j] > 0.5]
            minimal_box = minimal_box_containing_boxes(optimal_boxes)
            println("minimal_box ", minimal_box)
            append!(boxes, deepcopy(optimal_boxes[1:num_boxes]))
        elseif r > 2 && strategy in ["optsol-B", "optsol-C"]
            global_break = true
        end

        global_break && break

        boxes = remove_duplicate_boxes(boxes)
        println(format("round $(r): {}", length(boxes)))

        # println("solve LP$r")
        dir1 = dir * "LP$r-$method/"
        rm(dir1, force=true, recursive=true)
        output = ensembleopt(dir1, datasetdir, treeids, method * "-LP", timelimit, loglevel="Info", 
            preprocessing_level=preprocessing_level, boxes=boxes, z_vartype=strategy[end])
        save_result(output, dir1 * "result.csv")
        current_bound = output[:objval]
               
        # println("solve IP$r")
        if r == 1
            dir1 = dir * "IP$r-$method/"
            rm(dir1, force=true, recursive=true)
            output = ensembleopt(dir1, datasetdir, treeids, method, timelimit, loglevel="Info", 
                preprocessing_level=preprocessing_level, boxes=boxes, df_forest_folder=dir * "LP$r-$method/", z_vartype=strategy[end])
            save_result(output, dir1 * "result.csv")
            global_objval = output[:objval]
        end

        if (current_bound - global_objval) / (global_objval) <= 1e-04
            break
        end
    end
end

function get_result(dir)
    function get_objval(file)
        df = DataFrame(CSV.File(file))

        if df[end, "time_solve"] < 3600
            return df[end, "objval"]
        else
            return 123456789
        end
    end

    function compute_gapclosed(objval, LP, target)
        if LP <= objval
            return -1
        end

        return (LP - target) / (LP - objval)
    end

    function get_numleaves_per_tree(file)
        df = DataFrame(CSV.File(file))
        num_trees = maximum(df.treeid)
        return nrow(df) / num_trees
    end

    function get_numindepvars(file)
        df = DataFrame(CSV.File(file))
        return maximum(df.splitvar1)
    end

    result = Dict()
    result["objval"] = get_objval(dir * "IP1/result.csv")
    result["objLP1"] = get_objval(dir * "LP1/result.csv")
    result["objLP2"] = get_objval(dir * "LP2/result.csv")
    result["objLP-misic"] = get_objval(dir * "LP1-misic/result.csv")
    result["objLP-liftall"] = get_objval(dir * "LP1-liftall/result.csv")
    result["gapclosed"] = compute_gapclosed(result["objval"], result["objLP1"], result["objLP2"])
    result["treesize1"] = get_numleaves_per_tree(dir * "LP1/df_forest.csv")
    result["treesize2"] = get_numleaves_per_tree(dir * "LP2/df_forest.csv")
    result["numindepvars"] = get_numindepvars(dir * "IP1/df_forest.csv")
    result
end


function main()
    preprocessing_level = 0
    for getLPbox_level in [2]
        # for num_trees in [20]
        for num_trees in 10:10:20
            # for dataset in ["Concrete_T5000"]
            # for dataset in ["concrete_T1000_Depth10"]
            for dataset in ["concrete_T1000_Depth10", "redwine_T1000_Depth10", "assortment_T5000_depth10"]
            # for dataset in ["Concrete_T5000_GBT", "Redwine_T5000_GBT", "assortment_T5000_GBT"]
                dir = "result_1209/$(dataset)_$(num_trees)_$(preprocessing_level)_$(getLPbox_level)/"
                try
                    mkdir(dir)
                catch
                    nothing
                end

                datasetdir = "dataset/$(dataset)/"

                test(dir, datasetdir, num_trees, preprocessing_level, getLPbox_level)
                
                ##### get result
                result = get_result(dir)
                result["dataset"] = "$(dataset)_pre$(preprocessing_level)"
                result["num_trees"] = num_trees
                result["getLPbox_level"] = getLPbox_level
                CSV.write(dir * "result.csv", DataFrame(result))
            end
        end
    end
end

# main()

function collect_result(dir)
    files = Glob.glob("$(dir)*/result.csv")
    df_result = DataFrame()
    for file in files
        df = DataFrame(CSV.File(file))
        if occursin("_pre0", df[end,"dataset"]) || occursin("assortment", df[end,"dataset"])
            append!(df_result, DataFrame(df[end,:]), cols=:union)
        end
    end
    header = ["dataset", "num_trees", "getLPbox_level", "numindepvars", "treesize1", "treesize2", "objval", "objLP-misic", "objLP1", "objLP-liftall", "objLP2", "gapclosed"]
    select!(df_result, header)
    sort!(df_result, header)
    CSV.write(dir * "result.csv", df_result)
end
# collect_result("result_1209/")
# result = get_result("result_1208/assortment_T5000_depth10_20_1_2/")
# println(result)

# a = [[1,2,3], [4,5,6], [7,8,9]]
# a = (i for i in a)
# println(a)
# for i in Iterators.product(tuple(a))
#     @show i
# end

# boxes = [Dict(1=>[1,3], 2=>[1,2]), Dict(1=>[2,4], 3=>[1,2])]
# minimal_box = minimal_box_containing_boxes(boxes)
# println(minimal_box)

########################################################################################
# global_dir = "result_partition/"
# try
#     mkdir(global_dir)
# catch
#     nothing
# end

# datasets = ["Redwine_T5000_GBT", "redwine_T1000_Depth10", "hard_d5_i10_T1000", "hard_d5_i100_T1000"]
# for dataset in datasets, strategy in ["optsol-B", "optsol-C"], num_boxes in [1,2,3,4]
#     num_trees = occursin("GBT", dataset) ? 50 : 20
#     preprocessing_level = 1
#     dir = global_dir * "$(dataset)_T$(num_trees)_$(strategy)_$(num_boxes)/"
#     test(dir, "dataset/$(dataset)/", num_trees, preprocessing_level, strategy, num_boxes)
# end

######################################################################################
global_dir = "debug_partition_1221/"
try
    mkdir(global_dir)
catch
    nothing
end

# datasets = ["hard_d5_i10_T1000"] 
# datasets = ["redwine_T1000_Depth10", "Concrete_T5000", "assortment_T5000_depth10"]
# datasets = ["hard_d10_i10_T1000", "hard_d5_i100_T1000", "hard_d10_i100_T1000"]
# datasets = ["Concrete_T5000", "Concrete_T5000_GBT", "assortment_T5000_GBT", "assortment_T5000_depth10"]
# datasets = ["assortment_T5000_depth10"]
# datasets = ["Redwine_T5000_GBT"]


# datasets = ["redwinepair"]
# for dataset in datasets, num_trees in [150], strategy in ["LPsol-B"], num_boxes in [1]
#     # num_trees = 10
#     # num_trees = occursin("GBT", dataset) ? 50 : num_trees
#     # num_trees = occursin("hard", dataset) ? 10 : num_trees

#     preprocessing_level = occursin("assort", dataset) ? 2 : 1
#     dir = global_dir * "$(dataset)_T$(num_trees)_$(strategy)_$(num_boxes)/"
#     if preprocessing_level == 2 
#         dir = dir[1:end - 1] * "_pre$(preprocessing_level)/"
#     end

#     println()
#     println(dataset, " ", strategy, " ", num_boxes)
#     test(dir, "dataset/$(dataset)/", num_trees, preprocessing_level, strategy, num_boxes)
# end


## redpair, T=150
# parent = [0, 0, 0, 1, 1, 1, 0, 0, 8, 8, 8, 8]
# boxes_global = [Dict{Any,Any}(2 => [-Inf, 0.435]), 
#     Dict{Any,Any}(2 => [0.435, Inf],8 => [-Inf, 0.99292]), 
#     Dict{Any,Any}(2 => [0.435, Inf],8 => [0.99292, Inf]), 
#     Dict{Any,Any}(4 => [-Inf, 5.15],8 => [-Inf, 0.99264]),
#     Dict{Any,Any}(4 => [-Inf, 5.15],8 => [0.99264, Inf]),
#     Dict{Any,Any}(4 => [5.15, Inf],8 => [-Inf, 0.99264]),
#     Dict{Any,Any}(2 => [-Inf, 0.435],8 => [-Inf, 0.99292]), 
#     Dict{Any,Any}(2 => [-Inf, 0.435],8 => [0.99292, Inf]),
#     Dict{Any,Any}(3 => [-Inf, 0.475],1 => [-Inf, 11.65]),
#     Dict{Any,Any}(3 => [-Inf, 0.475],1 => [11.65, Inf]),
#     Dict{Any,Any}(3 => [0.475, Inf],1 => [-Inf, 11.65]),
#     Dict{Any,Any}(3 => [0.475, Inf],1 => [11.65, Inf])
#     ]

## redpair, T=500
parent = [0]
boxes_global = [
    nothing
    ]




dataset = "redwinepair"
num_trees = 500
dir = global_dir * "$(dataset)_T$(num_trees)_manual/"
treeids = Array(1:num_trees)

try
    mkdir(dir)
catch
    nothing
end

for r = 1:1
    for method in ["liftx-LP", "liftx"]

        dir1 = dir * "_$r-$method/"
        rm(dir1, force=true, recursive=true)

        df_forest_folder = ""
        if parent[r] > 0
            df_forest_folder = dir * "_$(parent[r])-$method/"
        end

        boxes = [boxes_global[r]]
        if isnothing(boxes[1])
            boxes = []
        end
        output = ensembleopt(dir1, "dataset/$(dataset)/", treeids, method, 3600, loglevel="Info", 
            preprocessing_level=1, boxes=boxes, z_vartype='B', df_forest_folder=df_forest_folder)
        save_result(output, dir1 * "result.csv")

        if occursin("-LP", method)
            df_forest = DataFrame(CSV.File(dir1 * "df_forest.csv"))
            df_sol = DataFrame(CSV.File(dir1 * "df_sol.csv"))
            df_forest.ysol = [sum(df_sol[j, "sol"] for j in df_forest[i,"leafbegin"]:df_forest[i,"leafend"]) for i = 1:nrow(df_forest)]
            df_forest.box_varindex = get_boxes_df_forest(df_forest, "splitvar1", "varindex", 0, typemax(Int))
            df_forest.nlsol = [compute_nonlinear_term(df_forest.box_varindex[i], df_sol) for i = 1:nrow(df_forest)]
            df_forest.ynlgap = abs.(df_forest.ysol .- df_forest.nlsol)
            df_forest.box = get_boxes_df_forest(df_forest, "splitvar", "splitval", -Inf, Inf)
            
            ### Sort depth 2 boxes by ynlgap
            df_filtered = DataFrame()
            mingap = 0.01
            df_filtered = sort(filter(row -> row.ynlgap >= mingap, df_forest), ["depth", "ynlgap", "ysol", "treeid"], rev=[false, true, true, false])
            CSV.write(dir1 * "df_filtered.csv", df_filtered)
        end
    end
end