# Julia v1.4.2
# JuMP v.0.21.3
# Gurobi v0.8.1 (optimizer v9.0.3)

using Formatting, JuMP, MathOptInterface, DataFrames, DataStructures
# using Gurobi # CPLEX does not have lazy constraint option in Julia
using SCIP, GLPK, Gurobi
const MOI = MathOptInterface
const GRB_ENV = Gurobi.Env()

# include("utils.jl") # collection of utility functions such as read_forest()
include("check_feasibility.jl")
include("extended.jl")
include("constraint.jl")
include("extended_optbox.jl")
include("partition_leaves.jl")
include("utils_JuMP.jl")
include("splitconstr.jl")

# const GRB_ENV = Gurobi.Env()


##### Rewrite the code 

function initialize_model(optimizer_attributes)
    ### Initialize model with Gurobi
    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV))
    for (name, value) in optimizer_attributes
        JuMP.set_optimizer_attribute(model, name, value)
    end

    ### Initialize model with SCIP
    # model = JuMP.Model(SCIP.Optimizer)
    # JuMP.set_optimizer_attribute(model, "limits/gap", optimizer_attributes["MIPGap"])
    # JuMP.set_optimizer_attribute(model, "limits/time", optimizer_attributes["TimeLimit"])

    # model = Model(GLPK.Optimizer)
    # set_attribute(model, "tm_lim", time_limit_sec * 1_000)
    # set_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)

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

function get_df_splitconstr_ver2(method, df_forest, df_vars)
    function get_empty_row()
        Dict(:forestindex=>Set(), :treeid=>0, :splitvar1=>0, :leaves=>SortedSet(), :lb=>0, :ub=>typemax(Int))
    end

    function drop_variables_with_pred_zero(df_splitconstr, df_vars)
        nr = nrow(df_splitconstr)

        yvars_positive = SortedSet(filter(row -> row.vartype == "y" && row.pred > 0, df_vars)[:, "varindex"])
        df_splitconstr.leaves = [intersect(row.leaves, yvars_positive) for row in eachrow(df_splitconstr)]
        filter!(row -> length(row.leaves) > 0, df_splitconstr)

        @info format("num constraints from {} -> {}", nr, nrow(df_splitconstr))

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
                constr[:leaves] = SortedSet(r.leafbegin:r.leafend)
                if df_forest.lchild[p] == r.nodeid
                    constr[:ub] = df_forest.varindex[p]
                else
                    constr[:lb] = df_forest.varindex[p]
                end
                append!(df_splitconstr, constr)
            end
        end
    elseif occursin("liftx", method)
        df_splitconstr, df_forest = get_df_splitconstr_ver2("misic", df_forest, df_vars)

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
            row.leaves = SortedSet([leaf for leaf in row.leaves if df_vars[leaf, "pred"] > 0.0])

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
        df_splitconstr, df_forest = get_df_splitconstr_ver2("misic", df_forest, df_vars)
        df_splitconstr = drop_variables_with_pred_zero(df_splitconstr, df_vars)
        df_splitconstr = covert_df_splitconstr_liftall(df_splitconstr)

        elapsed_time = time() - t_begin
        @info "get_df_splitconstr_ver2($method):: $(elapsed_time) seconds"
        return df_splitconstr, df_forest

        @assert true "not using this part"
        df_liftall = DataFrame()
        for _df1 in groupby(df_splitconstr, [:treeid, :splitvar1])
            treeid = _df1.treeid[1] 
            splitvar1 = _df1.splitvar1[1]
            df = DataFrame()

            ##### Combine splitconstriant with same RHS
            for _df2 in groupby(_df1, [:lb, :ub])
                forestindex_new = SortedSet()
                leaves_new = SortedSet()
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

                forestindex_new = SortedSet()
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
    @info "get_df_splitconstr_ver2($method):: $(elapsed_time) seconds"

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


# function get_df_splitconstr_liftall_tree(df_tree)
#     splitvars = df_tree.splitvar1
# end

# function get_df_splitconstr_liftall(df_forest)
#     arr_df_tree = df_forest_to_arr_df_tree(df_forest)
#     arr_df_splitconstr = [get_df_splitconstr_liftall_tree(df_tree) for df_tree in arr_df_tree]

# end

function get_df_splitconstr_with_liftall_where_rhs_is_liftx(df_splitconstr_liftall, df_splitconstr_liftx)
    # function filter_logic(row)
    #     for r_x in eachrow(df_splitconstr_liftx)
    #         if r_x.treeid == row.treeid && r_x.splitvar1 == row.splitvar1 && r_x.lb == row.lb && r_x.ub == row.ub
    #             return true
    #         end
    #     end
    #     return false
    # end

    ConstrIndices = Dict()
    for row in eachrow(df_splitconstr_liftx)
        if !haskey(ConstrIndices, row.treeid)
            ConstrIndices[row.treeid] = Dict()
        end
        if !haskey(ConstrIndices[row.treeid], row.splitvar1)
            ConstrIndices[row.treeid][row.splitvar1] = Dict()
        end
        if !haskey(ConstrIndices[row.treeid][row.splitvar1], row.lb)
            ConstrIndices[row.treeid][row.splitvar1][row.lb] = Set()
        end
        push!(ConstrIndices[row.treeid][row.splitvar1][row.lb], row.ub)
    end

    function filter_logic(row)
        return haskey(ConstrIndices, row.treeid) &&
            haskey(ConstrIndices[row.treeid], row.splitvar1) &&
            haskey(ConstrIndices[row.treeid][row.splitvar1], row.lb) &&
            row.ub in ConstrIndices[row.treeid][row.splitvar1][row.lb]
    end

    df_splitconstr = filter(row -> filter_logic(row), df_splitconstr_liftall)

    return df_splitconstr
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



function get_model(data, optimizer_attributes, formulation, lazylevel)
    @assert formulation in ["plain", "misic", "liftx", "liftall", "minimal"]
    @assert lazylevel in [0, 2, 3]

    model = initialize_model(optimizer_attributes)
    model = add_basic_formulation(model, data)

    if formulation in ["misic", "liftx", "liftall"] 
        df_splitconstr, df_forest = get_df_splitconstr_ver2(formulation, data["df_forest"], data["df_vars"])
        model = add_splitconstr(model, df_splitconstr, lazylevel)
        data["df_splitconstr"] = deepcopy(df_splitconstr)
        data["df_forest"] = deepcopy(df_forest)
    elseif formulation == "minimal"
        @info "Formulation is minimal"

        df_splitconstr_liftx, df_forest_liftx = get_df_splitconstr_ver2("liftx", data["df_forest"], data["df_vars"])
        df_splitconstr_liftall, df_forest_liftall = get_df_splitconstr_ver2("liftall", data["df_forest"], data["df_vars"])
        df_splitconstr = get_df_splitconstr_with_liftall_where_rhs_is_liftx(df_splitconstr_liftall, df_splitconstr_liftx)
        @info "liftx = $(nrow(df_splitconstr_liftx)) liftall = $(nrow(df_splitconstr_liftall)) minimal = $(nrow(df_splitconstr))"

        model = add_splitconstr(model, df_splitconstr, lazylevel)
        data["df_splitconstr"] = deepcopy(df_splitconstr)
        data["df_forest"] = deepcopy(df_forest_liftx)
    end

    model, data
end

function get_disjunctive_model(data, optimizer_attributes, formulation, lazylevel, boxes_splitval)
    """
    Disjuctive programming with boxes
    box 1, box2, ..., box B, complements
    """    
    function get_empty_row()
        Dict("treeid"=>0, "nodeid"=>0, "splitvar1"=>0, "splitval1"=>0, "vartype"=>"unknown")
    end

    function add_x_disj_variables(model, df_vars, boxes_splitval1)
        num_splitvars = maximum(df_vars.splitvar1)
        for (b, box) in enumerate(boxes_splitval1)
            for v = 1:num_splitvars
                lb = 0
                ub = typemax(Int)
                if v in keys(box)
                    lb = box[v][1]
                    ub = box[v][2]
                end
                df = filter(row -> row.vartype == "x" && row.splitvar1 == v && row.splitval1 >= lb && row.splitval1 <= ub, df_vars)
                df.vartype .= "x-disj"
                df.boxindex = [b for i=1:nrow(df)]
                append!(df_vars, df, cols=:union)
            end
        end
        df_vars.varindex = 1:nrow(df_vars)

        df_vars_boxes = filter(row -> row.vartype == "x-disj", df_vars)
        num_new_variables = nrow(df_vars_boxes)
        @variable(model, [i=1:num_new_variables], lower_bound=0, upper_bound=1)

        model, df_vars
    end

    function add_z_variables(model, df_vars, num_disj)
        @variable(model, z[i=1:num_disj], lower_bound=0, upper_bound=1)
        # @variable(model, z[i=1:num_disj], Bin)
        for b = 1:num_disj
            row = get_empty_row()
            row = merge(row, Dict("vartype"=>"z", "boxindex"=>b%num_disj))
            append!(df_vars, row, cols=:union)
        end
        df_vars.varindex = 1:nrow(df_vars)
        model, df_vars, z
    end    

    function add_logical_constraints_complement(model, df_vars, z)
        df_vars_x = filter(row -> row.vartype == "x", df_vars)
        variables = JuMP.all_variables(model)
        for df in groupby(df_vars_x, ["splitvar1"])
            @constraint(model, variables[df.varindex[end]] <= z[end])
        end
        model
    end    

    function add_logical_constraints_boxes(model, df_vars, z)
        df_vars_boxes = filter(row -> row.vartype == "x-disj", df_vars)
        variables = JuMP.all_variables(model)
    
        for df in groupby(df_vars_boxes, ["boxindex", "splitvar1"])
            @constraint(model, [i=1:nrow(df)-1], variables[df.varindex[i]] <= variables[df.varindex[i]+1])
            @constraint(model, variables[df.varindex[end]] <= z[df.boxindex[1]])
        end

        model
    end

    df_vars = data["df_vars"]
    df_vars.boxindex = zeros(Int, nrow(df_vars))
    
    boxes_splitval1 = [convert_box_format(box, "splitval1", df_vars) for box in boxes_splitval]
    boxes_varindex = [convert_box_format(box, "varindex", df_vars) for box in boxes_splitval]
    num_disj = length(boxes_splitval) + 1

    @info "boxes_splitval" boxes_splitval
    @info "boxes_splitval1" boxes_splitval1
    @info "boxes_varindex" boxes_varindex

    model, data = get_model(data, optimizer_attributes, "plain", 0)
    model = relax_integrality(model)
   
    ### Add variables
    model, df_vars = add_x_disj_variables(model, df_vars, boxes_splitval1) 
    model, df_vars, z = add_z_variables(model, df_vars, num_disj)

    ### Add sum(z) = 1 and logical constraints
    @constraint(model, sum(z[i] for i=1:num_disj) == 1)
    model = add_logical_constraints_complement(model, df_vars, z)
    model = add_logical_constraints_boxes(model, df_vars, z)

    ### Add split constraints
    t_begin = time()
    # (1) compute boxes for each leaf
    # println("(1) compute boxes for each leaf")
    df_boxes = deepcopy(data["df_forest"])
    df_boxes.box_varindex = get_boxes_df_forest(df_boxes, "splitvar1", "varindex", 0, typemax(Int))
    df_boxes = select(df_boxes, ["treeid", "nodeid", "box_varindex"])
    df_vars = leftjoin(df_vars, df_boxes, on=["treeid", "nodeid"])
    for row in eachrow(df_vars)
        if row.vartype == "y"
            for (b, box) in enumerate(boxes_varindex)
                if issubset_box(row.box_varindex, box)
                    row.boxindex = b
                    break
                end

                if !isempty_box(intersect_boxes([row.box_varindex, box]))
                    println(row)
                    println(b)
                    @assert true "error in partition leaves"
                end
            end
        end
    end

    # CSV.write("debug/df_vars.csv", df_vars)

    # (2) Create split constraitns for each box
    # (2-1) Generate split constrs 
    # println("(2-1) Generate split constrs ")
    df_splitconstr, df_forest = get_df_splitconstr_ver2(formulation, data["df_forest"], data["df_vars"])
    
    # (2-2) Get leaves by box
    # println("(2-2) Get leaves by box")
    df_leaves = filter(row -> row.vartype == "y", df_vars)
    # CSV.write("df_leaves.csv", df_leaves)
    leaves_by_treeid_box = Dict((df.treeid[1], df.boxindex[1]) => SortedSet(df.varindex) for df in groupby(df_leaves, ["treeid", "boxindex"]))
    keys_by_box = Dict(b => Set(collect(keys(box))) for (b,box) in enumerate(boxes_splitval))
    keys_by_box[0] = Set()
    for b = 1:num_disj-1
        keys_by_box[0] = union(keys_by_box[0], keys_by_box[b])
    end



    # (2-3) Split leaves by box for each constraint
    # println("(2-3) Split leaves by box for each constraint")
    # CSV.write("df_splitconstr.csv", df_splitconstr)

    ### Create dictionary for varindex
    x_disj_varindex = Dict((row.splitvar1, row.splitval1, row.boxindex) => row.varindex for row in eachrow(df_vars) if row.vartype in ["x", "x-disj"])
    df_z = filter(row -> row.vartype == "z", df_vars)
    z_varindex = df_z.varindex
    function find_varindex_in_b(df_vars, varindex, b)
        if varindex == 0
            return 0
        elseif varindex == typemax(Int)                     
            return z_varindex[b == 0 ? num_disj : b]

        #     df = filter(row -> row.vartype == "z", df_vars)
        #     if b == 0
        #         @assert x_disj_varindex[0, typemax(Int), b] == df.varindex[end]
        #         return df.varindex[end]
        #     else
        #         @assert x_disj_varindex[0, typemax(Int), b] == df.varindex[b]
        #         return df.varindex[b]
        #     end
        elseif b == 0
            return varindex
        else
            return x_disj_varindex[df_vars.splitvar1[varindex], df_vars.splitval1[varindex], b]
            # df = filter(row -> row.splitvar1 == df_vars.splitvar1[varindex]
            #                    && row.splitval1 == df_vars.splitval1[varindex] 
            #                    && row.boxindex == b, 
            #                    df_vars)
            # @assert x_disj_varindex[df_vars.splitvar1[varindex], df_vars.splitval1[varindex], b] == df.varindex[1]
            # return df.varindex[1]
        end
    end

    df_splitconstr.boxindex = zeros(Int, nrow(df_splitconstr))
    df_splitconstr_boxes = DataFrame(treeid=Int[], splitvar1=Int[], lb=Int[], ub=Int[], boxindex=Int[], leaves=SortedSet[])
    for row in eachrow(df_splitconstr), b = 0:num_disj-1
        if (row.treeid, b) in keys(leaves_by_treeid_box)
            row_new = Dict("lb"=>0, "ub"=>typemax(Int), 
                            "leaves"=>intersect(row.leaves, leaves_by_treeid_box[row.treeid, b]), 
                            "treeid"=>row.treeid, "splitvar1"=>row.splitvar1, "boxindex"=>b)
            if length(row_new["leaves"]) > 0
                if row["splitvar1"] in keys_by_box[b]
                    boxes = [df_vars.box_varindex[leaf] for leaf in row_new["leaves"]]
                    lb, ub = minimal_interval_containing_boxes(boxes, row.splitvar1, 0, typemax(Int))
                    @assert row.lb <= lb && ub <= row.ub "error"
                else
                    lb = row.lb
                    ub = row.ub 
                end
                
                row_new["lb"] = find_varindex_in_b(df_vars, lb, b)
                row_new["ub"] = find_varindex_in_b(df_vars, ub, b)

                append!(df_splitconstr_boxes, row_new)
            end
        end
    end
    
    model = add_splitconstr(model, df_splitconstr_boxes, lazylevel)
    elapsed_time = time() - t_begin
    @info "time to add split constraints:: $elapsed_time seconds"

    ### Convert model to MIP with dummy binary variables (to use lazy constraint feature when solving LP)
    model = covert_LP_to_MIP(model)

    header = ["index", "treeid", "forestindex", "splitvar1", "boxindex", "lb", "ub", "leaves"]
    select!(df_splitconstr_boxes, intersect(header, names(df_splitconstr_boxes)))
    data["df_vars"] = df_vars
    data["df_splitconstr"] = df_splitconstr
    data["df_splitconstr_boxes"] = df_splitconstr_boxes

    model, data
end

function fix_box(data, model, box)
    df_vars = data["df_vars"]

    for (var, lbub) in box
        lb = lbub[1]
        ub = lbub[2]
        variables = JuMP.all_variables(model)

        ## Fix lb part
        df = filter(row -> row.vartype == "x" && row.splitvar==var && row.splitval <= lb, df_vars)
        for i in df.varindex
            JuMP.fix(variables[i], 0.0, force=true)
        end

        ## Fix ub part
        df = filter(row -> row.vartype == "x" && row.splitvar==var && row.splitval >= ub, df_vars)
        for i in df.varindex
            JuMP.fix(variables[i], 1.0, force=true)
        end
    end

    model
end

function unfix_varialbes(data, model)
    df_vars = data["df_vars"]
    df = filter(row -> row.vartype in ["x", "y", "x-disj", "z"], df_vars)
    variables = JuMP.all_variables(model)
    for row in eachrow(df)
        if JuMP.is_fixed(variables[row.varindex])
            JuMP.unfix(variables[row.varindex])
            JuMP.set_lower_bound(variables[row.varindex], 0.0)
            JuMP.set_upper_bound(variables[row.varindex], 1.0)
        end
    end
    model
end
