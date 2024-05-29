using CSV, DataFrames, DelimitedFiles, Formatting, Glob


function last_integer_in_string(str::String)
    new_str = ""
    for ch in str
        if occursin(ch, "0123456789")
            new_str = new_str * ch
        else
            new_str = new_str * " "
        end
    end
    parse(Int, split(new_str)[end])
end


function sort_strings_by_last_integer(files::Array{String,1})
    num_files = length(files)
    last_integers = zeros(Int, length(files))
    for i = 1:num_files
        last_integers[i] = last_integer_in_string(files[i])
    end
    sorted_indexes = sort!([i for i = 1:num_files], by=i -> last_integers[i])
    files[sorted_indexes]    
end


function read_tree(treefile; header=nothing)
    if isnothing(header)
        header = [:nodeid, :lchild, :rchild, :splitvar, :splitval, :nodetype, :pred]
    end
    return DataFrame(CSV.File(treefile, header=header))
end


function reorder_nodes_in_DFS(df_tree)
    num_nodes = nrow(df_tree)
    lchild = df_tree.lchild
    rchild = df_tree.rchild
    old2new = zeros(num_nodes)

    ## Compute old2new
    n = 1
    nodecnt = 1
    old2new[n] = 1
    function update(n)
        if lchild[n] > 0 && old2new[lchild[n]] == 0
            nodecnt += 1
            old2new[lchild[n]] = nodecnt
            update(lchild[n])
        end

        if rchild[n] > 0 && old2new[rchild[n]] == 0
            nodecnt += 1
            old2new[rchild[n]] = nodecnt
            update(rchild[n])
        end
    end        
    update(1)

    ## Replace nodeids according to old2new
    for i = 1:num_nodes
        df_tree.nodeid[i] = old2new[i]
        if df_tree.lchild[i] > 0
            df_tree.lchild[i] = old2new[df_tree.lchild[i]]
        end
        if df_tree.rchild[i] > 0
            df_tree.rchild[i] = old2new[df_tree.rchild[i]]
        end
    end

    ## Remove nodes that is not connected to the root 
    df_tree = df_tree[df_tree.nodeid .> 0, :]

    ## Sort tree array by nodeid
    sort!(df_tree, [:nodeid])

    df_tree
end


##### Preprocessing for a tree - begin
function revalue_splitval(df_tree)
    # @assert true "only for assortment instances"

    for row in eachrow(df_tree)
        newval = 0
        if row.splitval > 1
            newval = 1
        else
            newval = round(row.splitval - (row.splitval % 0.025), digits=3)
        end 
        row.splitval = newval
    end
    df_tree
end


function get_splitval_interval_in_df_tree(df_tree, n, var, colname_var, colname_val)
    @assert n > 0
    @assert n <= nrow(df_tree)
    @assert colname_var in ["splitvar", "splitvar1"]
    @assert colname_val in ["splitval", "splitval1"]

    lb = -Inf
    ub = Inf
    if colname_val == "splitval1"
        lb = 0
        ub = typemax(Int)
    end

    parent = df_tree.parent[n]
    while parent > 0
        if df_tree[parent, colname_var] == var
            if df_tree.lchild[parent] == n
                ub = min(ub, df_tree[parent, colname_val])
            else
                lb = max(lb, df_tree[parent, colname_val])
            end
        end

        n = parent
        parent = df_tree.parent[n]
    end

    return [lb, ub]
end


function remove_nodes_with_empty_box(df_tree)
    if !("parent" in names(df_tree))
        df_tree = add_columns_to_df_tree(df_tree)
    end

    for n in 2:nrow(df_tree)
        # println(n)
        parent = df_tree.parent[n]
        islchild = df_tree.lchild[parent] == n
        splitvar = df_tree.splitvar[n]
        splitval = df_tree.splitval[n]
        if splitvar == 0
            continue
        end

        bounds = get_splitval_interval_in_df_tree(df_tree, n, splitvar, "splitvar", "splitval")
        # println(bounds)
        ch = 0
        if splitval <= bounds[1]
            ### lchild is empty 
            ### => connect parent to rchild
            ch = df_tree.rchild[n]
        elseif splitval >= bounds[2]
            ### rchild is empty 
            ### => connect parent to lchild
            ch = df_tree.rchild[n]
        end

        ### If ch > 0, then box[ch] is empty
        if ch > 0
            # @info "node $ch is empty. $bounds $splitval"
            if islchild
                df_tree.lchild[parent] = ch
            else
                df_tree.rchild[parent] = ch
            end
            df_tree.parent[ch] = parent
        end
    end

    df_tree = reorder_nodes_in_DFS(df_tree)
    df_tree = add_columns_to_df_tree(df_tree)
end

##### Preprocessing for a tree - end
function add_columns_to_df_tree(df_tree; add_parent=true, add_depth=true, add_predmax=true)
    df_tree = deepcopy(df_tree)
    
    num_nodes = nrow(df_tree)
    parent = -ones(Int, num_nodes)
    depth = -ones(Int, num_nodes)
    predmax = -Inf * ones(Float64, num_nodes)

    n = 1
    d = 0
    parent[n] = 0
    depth[n] = 0
    while n > 0
        lch = df_tree.lchild[n]
        rch = df_tree.rchild[n]
        if lch > 0 && parent[lch] == -1
            # non-leaf node and lchild is not updated
            parent[lch] = n
            depth[lch] = d + 1
            n = lch
            d += 1
        elseif rch > 0 && parent[rch] == -1
            # non-leaf node and lchild is updated but rchild is not updated
            parent[rch] = n
            depth[rch] = d + 1
            n = rch
            d += 1
        elseif lch > 0
            # non-leaf node and both lchild and rchild are updated 
            predmax[n] = max(predmax[lch], predmax[rch])
            n = parent[n]
            d -= 1
        else
            # leaf node 
            predmax[n] = df_tree.pred[n]
            n = parent[n]
            d -= 1
        end
    end
    
    if add_parent
        df_tree.parent = parent
    end
    if add_depth
        df_tree.depth = depth
    end
    if add_predmax
        df_tree.predmax = predmax
    end
    
    df_tree
end


function arr_df_tree_to_df_forest(arr_df_tree)
    df_forest = DataFrame()
    ptr_forest = [1]

    for t = 1:length(arr_df_tree)
        df_tree = arr_df_tree[t][:,:]
        df_tree.treeid = t * ones(Int, nrow(df_tree))
        append!(df_forest, df_tree, cols=:union)
        push!(ptr_forest, ptr_forest[end] + nrow(df_tree))
    end

    df_forest.index = 1:nrow(df_forest)

    df_forest = df_forest[:, [:index, :treeid, :nodeid, :lchild, :rchild, :parent, :depth, :splitvar, :splitval, :pred, :predmax]]

    df_forest, ptr_forest
end

function df_forest_to_arr_df_tree(df_forest)
    ptr_forest = get_pointers(df_forest, "treeid")
    num_trees = maximum(df_forest.treeid)
    return [df_forest[ptr_forest[t]:ptr_forest[t + 1] - 1, :] for t = 1:num_trees]
end


function df_forest_to_df_splits(df_forest)
    num_leaves = size(df_forest[df_forest.lchild .== 0, :])[1]

    df_splits = sort(unique(df_forest[:,["splitvar", "splitval"]]))
    df_splits = df_splits[2:end,:]
    df_splits.varindex = (num_leaves + 1):(num_leaves + size(df_splits)[1])

    num_splits = size(df_splits)[1]
    df_splits.splitvar1 = zeros(Int, num_splits)
    df_splits.splitval1 = zeros(Int, num_splits)
    df_splits.splitvar1[1] = 1
    df_splits.splitval1[1] = 1
    for i = 2:num_splits
        if df_splits.splitvar[i - 1] == df_splits.splitvar[i]
            df_splits.splitvar1[i] = df_splits.splitvar1[i - 1]  
            df_splits.splitval1[i] = df_splits.splitval1[i - 1] + 1
        else
            df_splits.splitvar1[i] = df_splits.splitvar1[i - 1] + 1
            df_splits.splitval1[i] = 1
        end
    end
    df_splits = df_splits[:, ["splitvar1", "splitval1", "varindex", "splitvar", "splitval"]]
    df_splits
end

function get_dataframes(df_forest)
    ##### nrows and ncols of df_forest
    nrows, ncols = size(df_forest)

    ##### Construct df_splits
    num_leaves = size(df_forest[df_forest.lchild .== 0, :])[1]

    df_splits = df_forest_to_df_splits(df_forest)
    
    ##### Add "splitvar1", "splitval1", "varindex" of df_splits to df_forest
    select!(df_forest, setdiff(names(df_forest), ["splitvar1", "splitval1", "varindex"]))
    df_forest = leftjoin(df_forest, df_splits, on=[:splitvar, :splitval])
    sort!(df_forest, ["treeid", "nodeid"])


    ##### Fill leaf rows in df_forest
    df_forest[df_forest.lchild .== 0, :splitvar1] .= 0
    df_forest[df_forest.lchild .== 0, :splitval1] .= 0
    df_forest[df_forest.lchild .== 0, :varindex] = 1:num_leaves

    ##### Add leafbegin, leafend to df_forest
    df_forest.leafbegin = zeros(Int, nrows)
    df_forest.leafend = zeros(Int, nrows)
    df_forest[df_forest.lchild .== 0, :leafbegin] = 1:num_leaves
    df_forest[df_forest.lchild .== 0, :leafend] = 1:num_leaves
    for i = reverse(1:nrows)
        if df_forest.lchild[i] != 0
            df_forest.leafbegin[i] = df_forest.leafbegin[i + df_forest.lchild[i] - df_forest.nodeid[i]]
            df_forest.leafend[i] = df_forest.leafend[i + df_forest.rchild[i] - df_forest.nodeid[i]]
        end
    end

    ##### Construct df_leaves
    df_leaves = sort(df_forest[df_forest.lchild .== 0, ["varindex", "treeid", "nodeid", "pred"]])

    ##### Construct df_vars
    df_vars = unique(sort(df_forest[:, ["varindex", "treeid", "nodeid", "pred", "splitvar1", "splitval1", "splitvar", "splitval"]]), 1)
    df_vars[num_leaves + 1:end, "treeid"] .= 0 
    df_vars[num_leaves + 1:end, "nodeid"] .= 0
    df_vars[num_leaves + 1:end, "pred"] .= 0
    df_vars.vartype = ["x or y" for i=1:nrow(df_vars)]
    df_vars.vartype[df_vars.treeid .> 0] .= "y"
    df_vars.vartype[df_vars.treeid .== 0] .= "x"

    ##### df_forest_index
    df_forest.index = 1:size(df_forest)[1]

    output = Dict()
    output["df_forest"] = df_forest
    output["df_splits"] = df_splits
    output["df_leaves"] = df_leaves
    output["df_vars"] = df_vars

    output
end


function read_data_in_df_forest(arr_treedir, arr_treeids, weights, preprocessing)
    df_forest = DataFrame()
    ptr_forest = Int[]
    num_trees_global = 0
    intercept_global = 0
    for T = 1:length(arr_treedir)
        treedir = arr_treedir[T]
        treeids = arr_treeids[T]

        intercept = 0
        ensemble_model = "random_forest"
        if occursin("GBT", treedir)
            ensemble_model = "GBT"
            @assert isfile(treedir * "intercept.txt") "there is no intercept.txt in $treedir."
            intercept = DelimitedFiles.readdlm(treedir * "intercept.txt")[1,1]
        end
        intercept *= weights[T]

        treefiles = glob(string(treedir, "*tree*.txt"))
        treefiles = sort_strings_by_last_integer(treefiles)
        treefiles = treefiles[treeids]
        @info treefiles

        arr_df_tree = []
        for treefile in treefiles
            ### read tree and reorder nodeids
            df_tree = read_tree(treefile)
            df_tree = reorder_nodes_in_DFS(df_tree)

            if "revalue_splitval" in preprocessing && occursin("assortment", treedir)
                df_tree = revalue_splitval(df_tree)
            end

            ### update pred
            df_tree.pred .*= weights[T]
            if ensemble_model == "random_forest"
                df_tree.pred ./= length(treeids)
            end
            minvalue = minimum(df_tree.pred)
            df_tree.pred .-= minvalue
            intercept += minvalue

            ### add columns
            df_tree = add_columns_to_df_tree(df_tree)
            if "remove_nodes_with_empty_box" in preprocessing
                df_tree = remove_nodes_with_empty_box(df_tree)
            end
            push!(arr_df_tree, df_tree)
        end

        _df_forest, _ptr_forest = arr_df_tree_to_df_forest(arr_df_tree)
        
        if isempty(df_forest)
            df_forest = deepcopy(_df_forest)
            ptr_forest = deepcopy(_ptr_forest)
        else
            _df_forest.treeid .+= num_trees_global
            append!(df_forest, _df_forest, cols=:union)
            _ptr_forest .+= ptr_forest[end] - 1
            pop!(ptr_forest)
            append!(ptr_forest, _ptr_forest)
        end
        num_trees_global += length(treeids)
        intercept_global += intercept
    end

    data = Dict("df_forest" => df_forest, "ptr_forest" => ptr_forest, "intercept" => intercept_global)
end


function read_data(arr_treedir, arr_treeids; weights=nothing, preprocessing=nothing)
    t_begin = time()

    if isnothing(weights)
        weights = [1.0]
    end
    if isnothing(preprocessing)
        preprocessing = []
    end
    @assert length(arr_treedir) == length(arr_treeids) "length(arr_treedir) != length(arr_treeids)"
    @assert length(arr_treedir) == length(weights) "length(arr_treedir) != length(weights)"

    data = read_data_in_df_forest(arr_treedir, arr_treeids, weights, preprocessing)

    ##### get df_forest (updated), df_splits, df_leaves, df_vars
    merge!(data, get_dataframes(data["df_forest"]))

    @info(format("read_data() (duration: {:.2f} sec)\n", time() - t_begin))

    data
end


function save_data(dir, data)
    try
        mkdir(dir)
    catch
        nothing
    end

    # CSV.write(dir * "df_splits.csv", data["df_splits"])
    # CSV.write(dir * "df_leaves.csv", data["df_leaves"])
    for (name, value) in data
        # if occursin("df_", dfname)
        if name in ["df_forest", "df_splitconstr", "df_splitconstr_boxes"]
            CSV.write("$(dir)$(name).csv", value)
        elseif name in ["df_vars", "df_sol"]
            header = ["varindex", "treeid", "nodeid", "pred", "splitvar1", "splitval1", "splitvar", "splitval", "vartype", "boxindex", "sol"]
            CSV.write("$(dir)$(name).csv", select(value, intersect(header, names(value))))
        elseif name == "intercept"
            io = open("$(dir)$(name).txt", "w")
            write(io, format("{:.2f}\n", value))
            close(io)
        end
    end
    # CSV.write(dir * "df_vars.csv", data["df_vars"])
    # CSV.write(dir * "df_forest.csv", data["df_forest"])
    # if "df_splitconstr" in keys(data)   
    #     CSV.write(dir * "df_splitconstr.csv", data["df_splitconstr"])
    # end
end

function read_data_from_saved(dir)
    df_forest = DataFrame(CSV.File(dir * "df_forest.csv"))
    intercept = readdlm(dir * "intercept.txt")[1,1]
    data = Dict()
    data = get_dataframes(df_forest)
    data["ptr_forest"] = get_pointers(df_forest, "treeid")
    data["intercept"] = intercept
    data
end