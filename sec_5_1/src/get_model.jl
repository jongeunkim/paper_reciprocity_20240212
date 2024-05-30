function minimal_box_containing_boxes(boxes)
    minimal_box = deepcopy(boxes[1])

    for box in boxes[2:end]
        to_be_deleted_keys = []
        for v in keys(minimal_box)
            if !(v in keys(box))
                push!(to_be_deleted_keys, v)
            else
                minimal_box[v] = [min(minimal_box[v][1], box[v][1]), max(minimal_box[v][2], box[v][2])]
                if minimal_box[v][2] >= Inf && minimal_box[v][1] <= -Inf
                    push!(to_be_deleted_keys, v)
                end
            end
        end

        for v in to_be_deleted_keys
            pop!(minimal_box, v)
        end
    end

    return minimal_box
end

function get_interval(df_forest, i, splitvar1)
    lb = 0
    ub = typemax(Int)

    while i > 0
        if df_forest.parent[i] <= 0
            break
        end

        nodeid = df_forest.nodeid[i]
        p = i + df_forest.parent[i] - nodeid
        if df_forest.splitvar1[p] == splitvar1
            if ub == typemax(Int) && nodeid == df_forest.lchild[p]
                # ub = df_forest.varindex[p]
                ub = min(ub, df_forest.varindex[p])
            elseif lb == 0 && nodeid == df_forest.rchild[p]
                # lb = df_forest.varindex[p]
                lb = max(lb, df_forest.varindex[p])
            end
        end
        
        # if lb > 0 && ub < typemax(Int)
        #     break
        # end
        
        i = p
    end

    lb, ub
end

function intersect_rectangles(rec1, rec2)
    rec = deepcopy(rec1)

    for (k, v) in rec2
        if k in keys(rec)
            rec[k] = [max(rec[k][1], rec2[k][1]), min(rec[k][2], rec2[k][2])]
        else
            rec[k] = rec2[k]
        end
    end
    
    rec
end

function get_boxes_df_forest(df_forest, varname, valname, minvalue, maxvalue)
    @assert varname in names(df_forest)
    @assert valname in names(df_forest)

    ### Initialize boxes
    boxes = [Dict() for i = 1:nrow(df_forest)] 

    ### For a nonleaf node, compute a box of a child by intersecting the node's box and the split
    for i = 1:nrow(df_forest)
        ### If it is a leaf node, then continue
        if df_forest[i, "lchild"] == 0
            continue
        end

        ### lch, rch: df_forest indexes of lchild and rchild
        nodeid = df_forest[i, "nodeid"]
        lch = i - nodeid + df_forest[i, "lchild"]
        rch = i - nodeid + df_forest[i, "rchild"]

        ### split information of parent
        var = df_forest[i, varname]
        val = df_forest[i, valname]

        ### Update boxes of lchild and rchild
        boxes[lch] = intersect_rectangles(boxes[i], Dict(var => [minvalue, val]))
        boxes[rch] = intersect_rectangles(boxes[i], Dict(var => [val, maxvalue]))
    end

    boxes
end

function relax_binary_variables(binary_variables)
    for v in binary_variables
        JuMP.unset_binary(v)
        JuMP.set_upper_bound(v, 1)
        JuMP.set_lower_bound(v, 0)
    end
end

function relax_integrality(model)
    binary_variables = [v for v in JuMP.all_variables(model) if JuMP.is_binary(v)] 
    relax_binary_variables(binary_variables)
    model
end

function covert_LP_to_MIP(model)
    """
    Add dummy binary variables and a dummy constraint
    This enables to use `lazy constraint` feature when solving LP
    """

    @variable(model, [i=1:2], Bin)
    var = JuMP.all_variables(model)[end-1:end]
    @constraint(model, var[1] + var[2] <= 1)
    
    model
end

function remove_dominated(df_splitconstr)
    df_splitconstr.dominated = [false for i = 1:nrow(df_splitconstr)]
    for df in groupby(df_splitconstr, "leaves")
        lb_max = maximum(df.lb)
        ub_min = minimum(df.ub)

        if lb_max == 0 && ub_min == typemax(Int)
            df.dominated .= true
        else
            for row in eachrow(df)
                if row.lb < lb_max || row.ub > ub_min
                    row.dominated = true
                end
            end
        end
    end
    filter!(row -> !row.dominated, df_splitconstr)
    select!(df_splitconstr, Not("dominated"))
    df_splitconstr
end

function convert_to_half_interval(df)
    ### lb
    arr_lbs = sort(unique(df.lb))
    num_lbs = length(arr_lbs)
    df_lb = DataFrame(lb=arr_lbs, ub=maximum(df.ub) * ones(Int, num_lbs), leaves=[DataStructures.SortedSet() for i = 1:num_lbs])
    for row in eachrow(df), row_lb in eachrow(df_lb)
        if row.lb >= row_lb.lb
            union!(row_lb.leaves, row.leaves)
        end
    end
    df_lb = remove_dominated(df_lb)

    ### ub
    arr_ubs = sort(unique(df.ub), rev=true)
    num_ubs = length(arr_ubs)
    df_ub = DataFrame(lb=minimum(df.lb) * ones(Int, num_ubs), ub=arr_ubs, leaves=[DataStructures.SortedSet() for i = 1:num_ubs])
    for row in eachrow(df), row_ub in eachrow(df_ub)
        if row.ub <= row_ub.ub
            union!(row_ub.leaves, row.leaves)
        end
    end
    df_ub = remove_dominated(df_ub)

    df_lb, df_ub
end

function convert_to_bounded_interval(df_lb, df_ub)
    df_lb = select(df_lb, ["lb", "leaves"])
    rename!(df_lb, "leaves" => "leaves_lb")
    df_lb.dummy = ones(Int, nrow(df_lb))
    df_ub = select(df_ub, ["ub", "leaves"])
    rename!(df_ub, "leaves" => "leaves_ub")
    df_ub.dummy = ones(Int, nrow(df_ub))
    
    df_bounded = outerjoin(df_lb, df_ub, on="dummy")
    df_bounded.leaves = [intersect(row.leaves_lb, row.leaves_ub) for row in eachrow(df_bounded)]
    select!(df_bounded, ["lb", "ub", "leaves"])
    filter!(row -> length(row.leaves) > 0, df_bounded)

    ### remove dominated ones
    df_bounded = remove_dominated(df_bounded)

    df_bounded
end

function convert_to_liftall(df)
    df_lb, df_ub = convert_to_half_interval(df)
    df_bounded = convert_to_bounded_interval(df_lb, df_ub)

    df_all = DataFrame(lb=Int[], ub=Int[], leaves=DataStructures.SortedSet[])
    append!(df_all, df_lb)
    append!(df_all, df_ub)
    append!(df_all, df_bounded)
    df_all = remove_dominated(df_all)

    df_all.index = 1:nrow(df_all)

    df_all
end

function initialize_model(optimizer_attributes)
    ### Initialize model with Gurobi
    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV))
    for (name, value) in optimizer_attributes
        JuMP.set_optimizer_attribute(model, name, value)
    end
    model
end

function add_basic_formulation(model, data)
    @assert "df_forest" in keys(data)
    @assert "df_vars" in keys(data)
    @assert "intercept" in keys(data)

    df_forest = data["df_forest"]
    df_vars = data["df_vars"]
    intercept = data["intercept"]

    ptr_forest = get_pointers(df_forest, "treeid")
    num_trees = length(ptr_forest) - 1
    num_leaves = sum(df_vars.vartype .== "y")
    num_splits = sum(df_vars.vartype .== "x")
    num_vars = nrow(df_vars)

    ### (1) Add variables [y, x]
    @variable(model, variables[1:num_vars], lower_bound = 0, upper_bound = 1)
    JuMP.set_binary.(variables[(num_leaves + 1):end])
    for row in eachrow(df_vars) 
        if row.vartype == "y" && row.pred <= 0.0
            JuMP.fix(variables[row.varindex], 0.0, force=true)
        end
    end

    ### (2) Set objective
    @objective(model, Max, intercept + sum(df_vars.pred[i] * variables[i] for i = 1:num_leaves))

    ### (3) Set constraints: one node per tree
    leaves_in_tree = Dict(df.treeid[1] => df.varindex for df in groupby(df_vars, ["treeid"]) if df.treeid[1] > 0)
    @constraint(model, constr_one_node[t=1:num_trees], sum(variables[i] for i in leaves_in_tree[t]) <= 1) 

    ### (4) Set constraints: between two split variables
    df = filter(row -> row.vartype == "x" && row.splitval1 > 1, df_vars)
    @constraint(model, constr_splitvar[i=1:nrow(df)], variables[df.varindex[i] - 1] <= variables[df.varindex[i]])

    model
end

function update_splitconstr_info_at_df_forest(df_forest, df_splitconstr)
    df_forest_new = df_forest[:,:]
    df_forest_new.conindex = zeros(Int, nrow(df_forest_new))
    for r in eachrow(df_splitconstr)
        for i in r.forestindex
            df_forest_new.conindex[i] = r.index
        end
    end
    df_forest_new
end

function get_df_splitconstr(method, df_forest, df_vars)
    function get_empty_row()
        Dict(:forestindex=>Set(), :treeid=>0, :splitvar1=>0, :leaves=>DataStructures.SortedSet(), :lb=>0, :ub=>typemax(Int))
    end

    function drop_variables_with_pred_zero(df_splitconstr, df_vars)
        nr = nrow(df_splitconstr)

        yvars_positive = DataStructures.SortedSet(filter(row -> row.vartype == "y" && row.pred > 0, df_vars)[:, "varindex"])
        df_splitconstr.leaves = [intersect(row.leaves, yvars_positive) for row in eachrow(df_splitconstr)]
        filter!(row -> length(row.leaves) > 0, df_splitconstr)

        @info Formatting.format("num constraints from {} -> {}", nr, nrow(df_splitconstr))

        df_splitconstr
    end

    @info "get_df_splitconstr(method=$method)"
    @assert any([occursin(m, method) for m in ["misic", "liftx", "liftall"]])   

    t_begin = time()

    df_splitconstr = DataFrame()
    if occursin("misic", method)
        for r in eachrow(df_forest)
            if r.nodeid > 1
                constr = get_empty_row()
                p = r.index - (r.nodeid - r.parent)
                push!(constr[:forestindex], r.index)
                constr[:treeid] = r.treeid
                constr[:splitvar1] = df_forest.splitvar1[p]
                constr[:leaves] = DataStructures.SortedSet(r.leafbegin:r.leafend)
                if df_forest.lchild[p] == r.nodeid
                    constr[:ub] = df_forest.varindex[p]
                else
                    constr[:lb] = df_forest.varindex[p]
                end
                append!(df_splitconstr, constr)
            end
        end
    elseif occursin("liftx", method)
        df_splitconstr, df_forest = get_df_splitconstr("misic", df_forest, df_vars)

        if !("box_varindex" in names(df_vars))
            df_boxes = deepcopy(df_forest)
            df_boxes.box_varindex = get_boxes_df_forest(df_boxes, "splitvar1", "varindex", 0, typemax(Int))
            df_boxes = select(df_boxes, ["treeid", "nodeid", "box_varindex"])
            df_vars = leftjoin(df_vars, df_boxes, on=["treeid", "nodeid"])
        end

        for row in eachrow(df_splitconstr)
            ### Existing
            row.lb, row.ub = get_interval(df_forest, collect(row.forestindex)[1], row.splitvar1)

            ### Narrow interval by dropping leaves with pred=0
            row.leaves = DataStructures.SortedSet([leaf for leaf in row.leaves if df_vars[leaf, "pred"] > 0.0])

            if length(row.leaves) > 0
                boxes = [df_vars.box_varindex[leaf] for leaf in row.leaves]
                minimal_box = minimal_box_containing_boxes(boxes)
                if minimal_box[row.splitvar1][1] < row.lb
                    @assert true "error"
                end
                if minimal_box[row.splitvar1][2] > row.ub
                    @assert true "error"
                end
                row.lb = minimal_box[row.splitvar1][1]
                row.ub = minimal_box[row.splitvar1][2]
            end
        end
    elseif occursin("liftall", method)
        df_splitconstr, df_forest = get_df_splitconstr("misic", df_forest, df_vars)
        df_splitconstr = drop_variables_with_pred_zero(df_splitconstr, df_vars)
        df_splitconstr = covert_df_splitconstr_liftall(df_splitconstr)

        elapsed_time = time() - t_begin
        @info "get_df_splitconstr($method):: $(elapsed_time) seconds"
        return df_splitconstr, df_forest

        @assert true "not using this part"
        df_liftall = DataFrame()
        for _df1 in groupby(df_splitconstr, [:treeid, :splitvar1])
            treeid = _df1.treeid[1] 
            splitvar1 = _df1.splitvar1[1]
            df = DataFrame()

            ##### Combine splitconstriant with same RHS
            for _df2 in groupby(_df1, [:lb, :ub])
                forestindex_new = DataStructures.SortedSet()
                leaves_new = DataStructures.SortedSet()
                for r in eachrow(_df2)
                    union!(forestindex_new, r.forestindex)
                    union!(leaves_new, r.leaves)
                end

                combined_constr = Dict(:forestindex=>forestindex_new, :leaves=>leaves_new, :lb=>_df2.lb[1], :ub=>_df2.ub[1])
                append!(df, combined_constr, cols=:union)
            end

            ##### Add leaf from other constraints
            sort!(df, [:lb, :ub])
            i = 1
            while (i+1 <= nrow(df) && df.lb[i] == 0 && df.lb[i+1] == 0)
                union!(df.leaves[i+1], df.leaves[i]) 
                i += 1
            end

            sort!(df, [:ub, :lb], rev=true)
            i = 1
            while (i+1 <= nrow(df) && df.ub[i] == typemax(Int) && df.ub[i+1] == typemax(Int))
                union!(df.leaves[i+1], df.leaves[i]) 
                i += 1
            end

            ##### Create constraints with bounded intervals
            df_left = df[df.lb .== 0, :]
            df_right = df[df.ub .== typemax(Int), :]
            for l in eachrow(df_left), r in eachrow(df_right)
                leaves_new = intersect(l.leaves, r.leaves)

                forestindex_new = DataStructures.SortedSet()
                for row in eachrow(_df1)
                    if issubset(row.leaves, leaves_new)
                        union!(forestindex_new, row.forestindex)
                    end
                end
                
                if !isempty(leaves_new)
                    if r.lb < l.ub 
                        bounded_constr = Dict(:forestindex=>forestindex_new, :leaves=>leaves_new, :lb=>r.lb, :ub=>l.ub)
                        append!(df, bounded_constr, cols=:union)
                    else
                        println(l)
                        println(r)
                        @assert false "there is empty box in the tree"
                    end
                end
            end

            ##### Remove implied constraints
            ### - There is no same RHS
            ### - If LHS is same but RHS is larger -> remove
            ### - for every constr, collect all smaller RHS constraints and check leaves are same or not
            df.index = 1:nrow(df)

            implied = zeros(Int, nrow(df))
            for i in 1:nrow(df)
                df_within = filter( row -> ( row.index != i && row.lb >= df.lb[i] && row.ub <= df.ub[i] ), df )

                for r in eachrow(df_within)
                    if issubset(df.leaves[i], r.leaves)
                        implied[i] = 1
                        break
                    end
                end
            end

            df = df[implied .== 0, :]
            df.interval_length = df.ub .- df.lb
            sort!(df, [:interval_length], rev=true)
            select!(df, [:forestindex, :leaves, :lb, :ub])            

            df.treeid = treeid * ones(Int, nrow(df))
            df.splitvar1 = splitvar1 * ones(Int, nrow(df))

            append!(df_liftall, df, cols=:union)
        end

        df_splitconstr = df_liftall
    end

    df_splitconstr.index = 1:nrow(df_splitconstr)

    df_forest = update_splitconstr_info_at_df_forest(df_forest, df_splitconstr)

    elapsed_time = time() - t_begin
    @info "get_df_splitconstr($method):: $(elapsed_time) seconds"

    df_splitconstr, df_forest
end

function covert_df_splitconstr_liftall(df_splitconstr, cols_groupby=nothing)
    @assert all([col in names(df_splitconstr) for col in ["treeid", "splitvar1", "leaves", "lb", "ub"]]) "missing columns"
    t_begin = time()

    if isnothing(cols_groupby)
        cols_groupby = ["treeid", "splitvar1"]
    else
        @assert all([col in names(df_splitconstr) for col in cols_groupby]) "missing columns"
    end
    

    df_splitconstr = filter(row -> length(row.leaves) > 0, df_splitconstr)

    df_splitconstr_new = DataFrame()
    for df in groupby(df_splitconstr, cols_groupby)
        # df = deepcopy(df)
        df_liftall = convert_to_liftall(df)
        for col in cols_groupby
            df_liftall[:, col] = df[1, col] .* ones(Int, nrow(df_liftall))
        end
        append!(df_splitconstr_new, df_liftall, cols=:union)
    end
    df_splitconstr_new.index = 1:nrow(df_splitconstr_new)

    elapsed_time = time() - t_begin
    @info "covert_df_splitconstr_liftall():: $(elapsed_time) seconds"

    df_splitconstr_new
end

function add_splitconstr(model, df_splitconstr, lazylevel)
    @info "lazylevel = $(lazylevel)"
    variables = JuMP.all_variables(model)

    cnt = 0
    for r in eachrow(df_splitconstr)
        if length(r.leaves) == 0
            nothing
        elseif r.lb < r.ub
            lb = r.lb > 0 ? variables[r.lb] : 0
            ub = r.ub < typemax(Int) ? variables[r.ub] : 1
            c = @constraint(model, sum(variables[j] for j in r.leaves) <= ub - lb, base_name="constr_branching_$cnt")
            cnt += 1
            if lazylevel > 0
                MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), c, lazylevel)
            end
        else
            for j in r.leaves
                JuMP.fix(variables[j], 0, force=true)
            end
        end
    end
    
    model
end

function get_model(data::Dict, optimizer_attributes::Dict, formulation::String, lazylevel::Int, is_LP::Bool)
    @assert formulation in ["TEOM", "TEOR", "TEOC"]
    @assert lazylevel in [0, 2, 3]

    model = initialize_model(optimizer_attributes)
    model = add_basic_formulation(model, data)

    formulation_to_internal_name = Dict("TEOM" => "misic", "TEOR" => "liftx", "TEOC" => "liftall")
    df_splitconstr, df_forest = get_df_splitconstr(formulation_to_internal_name[formulation], data["df_forest"], data["df_vars"])
    model = add_splitconstr(model, df_splitconstr, lazylevel)
    data["df_splitconstr"] = deepcopy(df_splitconstr)
    data["df_forest"] = deepcopy(df_forest)

    if is_LP
        model = relax_integrality(model)
        model = covert_LP_to_MIP(model)
    end

    model, data
end
