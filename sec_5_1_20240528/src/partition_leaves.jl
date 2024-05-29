using CSV, DataFrames, Dates, DelimitedFiles, Formatting, Glob, Logging, Missings, Printf, Random


include("read_data.jl")
include("utils.jl")

function get_rectangle_in_tree(df_tree, node)
    box = Dict()
    while df_tree.parent[node] > 0
        parent = df_tree.parent[node]
        splitvar = df_tree.splitvar[parent]
        splitval = df_tree.splitval[parent]
        islchild = df_tree.lchild[parent] == node
        
        if islchild
            if splitvar in keys(box)
                box[splitvar] = [box[splitvar][1], min(splitval, box[splitvar][2])]
            else
                box[splitvar] = [-Inf, splitval]
            end
        else 
            if splitvar in keys(box)
                box[splitvar] = [max(splitval, box[splitvar][1]), box[splitvar][2]]
            else
                box[splitvar] = [splitval, Inf]
            end
        end

        node = parent
    end
    box
end 

function find_leaf_containing_box(df_tree, box; allow_nonleaf=false)
    # println(box)
    
    node = 1
    box_node = Dict()
    while df_tree.lchild[node] > 0
        splitvar = df_tree.splitvar[node]
        splitval = df_tree.splitval[node]
        # println("$node: X$splitvar <= $splitval")

        if !(splitvar in keys(box))
            return 0
        end

        (splitval_lb, splitval_ub) = (box[splitvar][1], box[splitvar][2])
        # println("[$splitval_lb, $splitval_ub]")
        
        if splitval_ub <= splitval
            # box_node = intersect_rectangles(box_node, Dict(splitvar=>[])
            node = df_tree.lchild[node]
        elseif splitval_lb >= splitval
            # box_node = Dict()
            node = df_tree.rchild[node]
        else
            return 0
        end
    end

    return node
end

function get_needed_splits(box, box_leaf)
    splits = []
    for v in keys(box)
        if box[v][1] > -Inf
            if !(v in keys(box_leaf)) || box_leaf[v][1] < box[v][1]
                push!(splits, [v, box[v][1], "right"])
            end
        end

        if box[v][2] < Inf
            if !(v in keys(box_leaf)) || box_leaf[v][2] > box[v][2]
                push!(splits, [v, box[v][2], "left"])
            end
        end

    end
    splits
end

function partition_leaf(df_tree, leaf, splitvar, splitval)
    num_nodes = nrow(df_tree)
    treeid = df_tree.treeid[leaf]
    pred = df_tree.pred[leaf]

    ### Change leaf
    df_tree.splitvar[leaf] = splitvar
    df_tree.splitval[leaf] = splitval
    df_tree.lchild[leaf] = num_nodes + 1
    df_tree.rchild[leaf] = num_nodes + 2

    ### Add children
    lch = Dict("treeid" => treeid, "nodeid" => num_nodes + 1, "lchild" => 0, "rchild" => 0, "splitvar" => 0, "splitval" => 0, "pred" => pred, "parent" => leaf)
    rch = Dict("treeid" => treeid, "nodeid" => num_nodes + 2, "lchild" => 0, "rchild" => 0, "splitvar" => 0, "splitval" => 0, "pred" => pred, "parent" => leaf)
    append!(df_tree, lch, cols=:union)
    append!(df_tree, rch, cols=:union)

    df_tree
end

function partition_leaves_in_tree(df_tree, box)
    function compute_partition_number(box_leaf, box)
        """
        partition number = required number of additional leaves
        """
        if isempty_box(intersect_rectangles(box_leaf, box))
            return 0
        end

        if issubset_box(box_leaf, box)
            return 0
        end

        partition_number = 0
        for (k, v) in box
            for splitval in v
                if splitval > -Inf && splitval < Inf
                    if k in keys(box_leaf)
                        if box_leaf[k][1] < splitval && splitval < box_leaf[k][2]
                            partition_number += 1
                        end
                    else
                        partition_number += 1
                    end
                end
            end
        end
        return partition_number
    end

    function find_split(box_leaf, box)
        for (k, v) in box
            for splitval in v
                if splitval > -Inf && splitval < Inf
                    if k in keys(box_leaf)
                        if box_leaf[k][1] < splitval && splitval < box_leaf[k][2]
                            return (k, splitval)
                        end
                    else
                        return (k, splitval)
                    end
                end
            end
        end


        return (0, 0.0)
    end

    function create_row()
        Dict{Any,Any}("nodeid"=>0, "lchild"=>0, "rchild"=>0, "splitvar"=>0, "splitval"=>0.0, "pred"=>0.0)
    end

    function split_leaf(df_tree, n, splitvar, splitval)
        lch = nrow(df_tree) + 1
        rch = nrow(df_tree) + 2

        ### edit the current row
        df_tree[n,"lchild"] = lch
        df_tree[n,"rchild"] = rch
        df_tree[n,"splitvar"] = splitvar
        df_tree[n,"splitval"] = splitval

        ### add lchild
        row = create_row()
        row["nodeid"] = lch
        row["pred"] = df_tree[i, "pred"]
        row["box"] = intersect_rectangles(df_tree[i,"box"], Dict(splitvar=>[-Inf, splitval]))
        append!(df_tree, row, cols=:union)

        ### add rchild
        row = create_row()
        row["nodeid"] = rch
        row["pred"] = df_tree[i, "pred"]
        row["box"] = intersect_rectangles(df_tree[i,"box"], Dict(splitvar=>[splitval, Inf]))
        append!(df_tree, row, cols=:union)

        df_tree
    end

    df_tree = deepcopy(df_tree)

    ### Search all leaves that intersect with box
    df_tree.box = get_boxes_df_forest(df_tree, "splitvar", "splitval", -Inf, Inf)
    df_tree.partition_number = [df_tree.lchild[i] == 0 ? compute_partition_number(df_tree.box[i], box) : 0 for i in 1:nrow(df_tree)]
    # println("sum of partition number = ", sum(df_tree.partition_number))
    # println("nrow (before) = ", nrow(df_tree))

    ### For each leaf, partition it
    i = 1
    while i <= nrow(df_tree)
        if df_tree.lchild[i] == 0
            partition_number = compute_partition_number(df_tree.box[i], box)
            # println("$i, $partition_number, ", df_tree.box[i], ", ", box)
            if partition_number > 0
                (splitvar, splitval) = find_split(df_tree.box[i], box)
                df_tree = split_leaf(df_tree, i, splitvar, splitval)
            end
        end
        i += 1
    end

    # println("nrow (after) = ", nrow(df_tree))
    # CSV.write("debug/df_tree.csv", df_tree)
    df_tree
end


function partition_leaves(data, box)
    t_begin = time()

    arr_df_tree = df_forest_to_arr_df_tree(data["df_forest"])
    header = ["nodeid", "lchild", "rchild", "splitvar", "splitval", "pred"]
    arr_df_tree = [select(partition_leaves_in_tree(df_tree, box), header) for df_tree in arr_df_tree]
    arr_df_tree = [reorder_nodes_in_DFS(df_tree) for df_tree in arr_df_tree]
    arr_df_tree = [add_columns_to_df_tree(df_tree) for df_tree in arr_df_tree]
    arr_df_tree = [remove_nodes_with_empty_box(df_tree) for df_tree in arr_df_tree]
    df_forest, ptr_forest = arr_df_tree_to_df_forest(arr_df_tree)


    intercept = data["intercept"]
    
    data = get_dataframes(df_forest)
    data["ptr_forest"] = ptr_forest
    data["intercept"] = intercept
    
    elapsed_time = time() - t_begin
    @info "time to partition leaves:: $elapsed_time seconds"

    data
end

function retain_leaves_by_box_in_tree(df_tree, box, in_or_out)
    @assert in_or_out in ["in", "out"]

    function make_row_empty(row)
        row.lchild = 0
        row.rchild = 0
        row.splitvar = 0
        row.splitval = 0
        row.pred = 0.0
    end

    df_tree.box = get_boxes_df_forest(df_tree, "splitvar", "splitval", -Inf, Inf)

    for row in eachrow(df_tree)
        if !issubset_box(row.box, box) && !isempty_box(intersect_boxes([row.box, box]))
            @assert true "not supposed to"
        end

        if in_or_out == "in"
            # if !issubset_box(row.box, box)
            if isempty_box(intersect_boxes([row.box, box]))
                make_row_empty(row)
            end
        else
            if issubset_box(row.box, box)
                make_row_empty(row)
            end
        end

    end

    df_tree
end


function retain_leaves_by_box(data, box, in_or_out)
    """
    remove leaves that does not intersect the box
    """
    t_begin = time()
    @assert in_or_out in ["in", "out"]
    
    arr_df_tree = df_forest_to_arr_df_tree(data["df_forest"])
    header = ["nodeid", "lchild", "rchild", "splitvar", "splitval", "pred"]
    arr_df_tree = [select(retain_leaves_by_box_in_tree(df_tree, box, in_or_out), header) for df_tree in arr_df_tree]
    arr_df_tree = [reorder_nodes_in_DFS(df_tree) for df_tree in arr_df_tree]
    arr_df_tree = [add_columns_to_df_tree(df_tree) for df_tree in arr_df_tree]
    arr_df_tree = [remove_nodes_with_empty_box(df_tree) for df_tree in arr_df_tree]
    df_forest, ptr_forest = arr_df_tree_to_df_forest(arr_df_tree)

    intercept = data["intercept"]

    # save_data("debug_bnb_debug/", data)
    
    data = get_dataframes(df_forest)
    data["ptr_forest"] = ptr_forest
    data["intercept"] = intercept

    elapsed_time = time() - t_begin
    @info "time to retain_leaves_by_box:: $elapsed_time seconds"
    
    data
end


# data = read_data(["dataset/Redwine_T5000/"], [Array(1:10)])
# CSV.write("debug/df_forest_before.csv", data["df_forest"])
# box = Dict{Any,Any}(7 => [-Inf, 65.6], 10 => [-Inf, 0.675])
# data = partition_leaves(data, box)
# CSV.write("debug/df_forest_after.csv", data["df_forest"])
