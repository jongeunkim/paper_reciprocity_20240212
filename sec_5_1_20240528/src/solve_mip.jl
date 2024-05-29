# Julia v1.4.2
# JuMP v.0.21.3
# Gurobi v0.8.1 (optimizer v9.0.3)

using Formatting, JuMP, MathOptInterface, DataFrames, DataStructures
using Gurobi # CPLEX does not have lazy constraint option in Julia
const MOI = MathOptInterface

# include("utils.jl") # collection of utility functions such as read_forest()
include("check_feasibility.jl")
include("extended.jl")
include("constraint.jl")
include("extended_optbox.jl")
include("partition_leaves.jl")
include("utils_JuMP.jl")
include("splitconstr.jl")

const GRB_ENV = Gurobi.Env()

function add_constraint_box(model, box, df_forest, vartype)
    # @info box

    @assert vartype in ['B', 'C']
    if vartype == 'C'
        @variable(model, [1:1], lower_bound = 0, upper_bound = 1) 
    else
        @variable(model, [1:1], Bin) 
    end
    z = JuMP.all_variables(model)[end]
    
    df_forest.box = get_boxes_df_forest(df_forest, "splitvar", "splitval", -Inf, Inf)
    arr_df_tree = df_forest_to_arr_df_tree(df_forest)

    variables = JuMP.all_variables(model)

    for df_tree in arr_df_tree
        boxes = [row.box for row in eachrow(df_tree) if row.lchild == 0 && issubset_box(row.box, box)]
        # @info "selected boxes" boxes
        minimal_box = minimal_box_containing_boxes(boxes)
        @assert box == minimal_box format("{}", minimal_box)

        leaves = [row.varindex for row in eachrow(df_tree) if row.lchild == 0 && issubset_box(row.box, box)]

        @constraint(model, z == sum(variables[j] for j in leaves))
    end

    model
end

function add_logarithmic_model(model, df_splits)
    function get_logarithm_constraints(num_elements)
        function get_recursive_half_splits(num_elements)
            function split_in_half(arr)
                if length(arr) == 1
                    return [arr, [0]]
                else
                    split = div(length(arr) + 1, 2)
                    [arr[1:split], arr[(split + 1):end]]
                end
            end
        
            @assert num_elements > 0
        
            collection = [1:num_elements]
        
            check = true
            while check
                collection_new = []
                for arr in collection
                    collection_new = vcat(collection_new, split_in_half(arr))
                end
                collection = deepcopy(collection_new)
        
                check = any([length(arr) > 1 for arr in collection])
            end
        
        
            collection = [arr[1] for arr in collection]
        
            return collection
        end
        collection = get_recursive_half_splits(num_elements)
    
        num_constrs = Int(log2(length(collection)))
        
        output = []
        for i = 1:num_constrs
            block_size = 2^(i - 1)
            
            elements = []
            j = 1
            while j <= length(collection)
                elements = vcat(elements, collection[j:j + block_size - 1])
                j += 2 * block_size
            end
            push!(output, elements)    
        end
    
        output = [setdiff(Set(arr), Set(0)) for arr in output]
        output
    end
    

    ### Relax integrality (binary -> [0,1]) 
    @debug "relax model"
    # model = JuMP.relax_integrality(model)
    variables = JuMP.all_variables(model)
    binary_variables = [v for v in variables if JuMP.is_binary(v)]
    for v in binary_variables
        JuMP.unset_binary(v)
        JuMP.set_upper_bound(v, 1)
        JuMP.set_lower_bound(v, 0)
    end


    ### 
    @debug "add var and constraints"
    for df in groupby(df_splits, ["splitvar1"])
        @debug "add var and constraints" df.splitvar1[1]
        num_intervals = nrow(df) + 1
        lhs_intervals = get_logarithm_constraints(num_intervals)

        @variable(model, [1:length(lhs_intervals)], Bin)
        variables = JuMP.all_variables(model)

        for (i, intervals) in enumerate(lhs_intervals)

            lhs = 0
            for j in intervals
                if j == 1
                    lhs += variables[df.varindex[j]]
                elseif j == num_intervals
                    lhs += 1 - variables[df.varindex[j - 1]]
                else
                    lhs += variables[df.varindex[j]] - variables[df.varindex[j - 1]]
                end
            end
            @constraint(model, lhs <= variables[end - i + 1])

            lhs = 0
            for j in setdiff(Set(1:num_intervals), intervals)
                if j == 1
                    lhs += variables[df.varindex[j]]
                elseif j == num_intervals
                    lhs += 1 - variables[df.varindex[j - 1]]
                else
                    lhs += variables[df.varindex[j]] - variables[df.varindex[j - 1]]
                end
            end
            @constraint(model, lhs <= 1 - variables[end - i + 1])
        end
    end

    model
end


function solve_mip(dir::String, data, method::String, timelimit::Number; 
                    MIPGap=1e-04, df_sol_fixed=DataFrame(), optbox=Dict(), 
                    boxes=[], z_vartype='C') 
    """
    Function milp_lazycb() solves a tree ensemble optimization problem.
    
    INPUT:
    dir = directory to save the result
    data = "forest" and "ptr_forest"
    """
    
    t_begin = time()
    output = Dict()

    ptr_forest = data["ptr_forest"]
    df_forest = data["df_forest"]
    df_leaves = data["df_leaves"]
    df_splits = data["df_splits"]

    # Initialize MIP model
    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV))
    JuMP.set_optimizer_attributes(model, "LogFile" => dir * "gurobi.log", "LogToConsole" => 0, "TimeLimit" => timelimit, "MIPGap" => MIPGap)
    
    # Add variables
    num_trees = length(ptr_forest) - 1
    num_leaves = nrow(df_leaves)
    num_splits = nrow(df_splits)
    num_vars = num_leaves + num_splits
    @variable(model, variables[1:num_vars], lower_bound = 0, upper_bound = 1)
    if !occursin("-LP", method)
        JuMP.set_binary.(variables[num_leaves + 1:end])
    end

    if occursin("-log", method)
        @info "add_logarithmic_model"
        model = add_logarithmic_model(model, df_splits)
        @info "add_logarithmic_model - done"
    end

    if !isempty(df_sol_fixed)
        for row in eachrow(df_sol_fixed)
            if row.treeid == 0
                JuMP.fix(variables[row.varindex], row.sol, force=true)
            end
        end
    end

    # Set objective
    @objective(model, Max, data["intercept"] + sum(data["df_vars"].pred[i] * variables[i] for i = 1:num_leaves))

    # Set constraints: one node per tree
    leaves_in_tree = [df_leaves.varindex[df_leaves.treeid .== t] for t = 1:num_trees]
    @constraint(model, constr_one_node[t=1:num_trees], sum(variables[i] for i in leaves_in_tree[t]) <= 1)

    # Set constraints: between two split variables
    ind = df_splits.varindex[df_splits.splitval1 .> 1, :]
    @constraint(model, constr_splitvar[i=1:length(ind)], variables[ind[i] - 1] <= variables[ind[i]])
    
    # Set constraints: if there are leaves identical to a box, then add equality constraints
    if length(boxes) > 0
        @info "add_constraint_box()"
        @info boxes
        for box in boxes
            model = add_constraint_box(model, box, df_forest, z_vartype)
        end

        # for subboxes in powerset(boxes)
        #     if length(subboxes) > 0
        #         box = intersect_boxes(subboxes)
        #         if !isempty_box(box)
        #             model = add_constraint_box(model, box, df_forest, z_vartype)
        #         end
        #     end
        # end
    end

    ###### Callback / Split Constr
    @debug "solve_mip():: callback"
    num_call = 0
    num_lazy = 0
    num_user = 0
    time_cb = 0
    if method == "plain"
        nothing
    elseif method in ["-extvars", "-extvars-LP"]
        df_extvars = get_df_extvars(ptr_forest, df_forest, df_leaves, df_splits)
        CSV.write(dir * "df_extvars.csv", df_extvars)
        @variable(model, extvars[1:nrow(df_extvars)], lower_bound = 0, upper_bound = 1)
        model = add_constraints_for_extvars(model, extvars, variables, df_extvars)
    elseif occursin("-LP", method) || occursin("-fixed", method)
        df_splitconstr, df_forest, arr_zero_leaves = get_df_splitconstr(method, df_forest)
        @info format("length(arr_zero_leaves) = {}", length(arr_zero_leaves))
        for i in arr_zero_leaves
            JuMP.set_upper_bound(variables[i], 0)
        end

        select!(df_splitconstr, [:index, :treeid, :splitvar1, :forestindex, :leaves, :lb, :ub])
        CSV.write(dir * "df_splitconstr.csv", df_splitconstr)
        CSV.write(dir * "df_forest.csv", df_forest)

        model = add_constraints(model, variables, df_splitconstr)
    elseif method in ["liftall", "liftall-nocut", "liftall-optbox"]
        df_splitconstr, df_forest, arr_zero_leaves = get_df_splitconstr(method, df_forest)
        select!(df_splitconstr, [:index, :treeid, :splitvar1, :forestindex, :leaves, :lb, :ub])
        CSV.write(dir * "df_splitconstr.csv", df_splitconstr)
        CSV.write(dir * "df_forest.csv", df_forest)

        @info format("length(arr_zero_leaves) = {}", length(arr_zero_leaves))
        for i in arr_zero_leaves
            JuMP.set_upper_bound(variables[i], 0)
        end

        df_vars = data["df_vars"]
        if method == "liftall-optbox"
            @info method
            @info optbox
            
            # println("debug 1")
            df_extvars_optbox = get_extvars_optbox(optbox, df_leaves, df_forest)
            # @info df_extvars_optbox
            CSV.write(dir * "df_extvars_optbox.csv", df_extvars_optbox)

            # println("debug 2")
            df_vars = merge_df_vars_df_extvars_optbox(df_vars, df_extvars_optbox)
            df_vars.varindex = 1:nrow(df_vars)
            CSV.write(dir * "df_vars.csv", df_vars)

            model = add_variables_for_df_extvars_optbox(model, df_vars)
            variables = JuMP.all_variables(model)

            model = add_y_constraint_for_df_extvars_optbox(model, df_vars)
            df_splitconstr = lift_df_extvars_optbox_to_df_splitconstr(df_splitconstr, df_vars)
            CSV.write(dir * "df_splitconstr2.csv", df_splitconstr)

            # println("debug end")
        end

        df_sol = df_vars
        df_sol.sol = zeros(Int, size(df_sol)[1])
        df_sol.fsol = zeros(Float64, size(df_sol)[1])
        min_violation = 0.1

        function my_callback_checkall(cb_data, cb_where::Cint)
            t_begin = time()

            ##### You can select where the callback is run
            if cb_where != Gurobi.GRB_CB_MIPSOL && cb_where != Gurobi.GRB_CB_MIPNODE
                return
            end

            if method == "liftall-nocut" && cb_where == Gurobi.GRB_CB_MIPNODE
                return
            end

            ##### You can query a callback attribute using GRBcbget
            if cb_where == GRB_CB_MIPNODE
                resultP = Ref{Cint}()
                GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
                if resultP[] != GRB_OPTIMAL
                    return  # Solution is something other than optimal.
                end

                resultP = Ref{Cdouble}()
                GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_NODCNT, resultP)
                if resultP[] < 1
                    @info format("resultP[] = {}", resultP[])
                    min_violation = 0.01
                else
                    return
                end

            end

            num_call += 1

            if cb_where == Gurobi.GRB_CB_MIPSOL
                ##### Before querying `callback_value`, you must call:
                Gurobi.load_callback_variable_primal(cb_data, cb_where)
                for i = 1:size(df_sol)[1]
                    df_sol[i, :sol] = convert(Int, round(callback_value(cb_data, variables[i])))
                end        

                ##### Collect active leaves
                df_active_leaves = df_sol[1:num_leaves, :]
                df_active_leaves = df_active_leaves[df_active_leaves.sol .== 1, :]

                ##### Collect violated constraints
                max_depth = typemax(Int)
                max_violated_per_tree = 1
                if method == "maxclq"
                    max_violated_per_tree = typemax(Int)
                else
                    nothing
                end

                arr_conindex = []
                for leaf in eachrow(df_active_leaves)
                    t = leaf.treeid
                    leafindex = leaf.varindex
                    n = ptr_forest[t]
                    if leafindex < df_forest[n, "leafbegin"] || leafindex > df_forest[n, "leafend"]
                        @info "error in df_active_leaves"
                    end
                                        
                    while df_forest.lchild[n] > 0
                        lch = n + df_forest.lchild[n] - df_forest.nodeid[n]
                        if leafindex <= df_forest[lch, "leafend"]
                            ### leaf is on the left
                            if df_sol.sol[df_forest.varindex[n]] == 0
                                push!(arr_conindex, df_forest.conindex[lch])
                                break
                            end

                            ### move to the right
                            n = lch
                        else
                            ### leaf is on the right
                            rch = n + df_forest.rchild[n] - df_forest.nodeid[n]
                            if df_sol.sol[df_forest.varindex[n]] == 1
                                push!(arr_conindex, df_forest.conindex[rch])
                                break
                            end

                            ### move to the right
                            n = rch
                        end
                    end
                end

                ##### Add lazy constraints
                for i in arr_conindex
                    lb = df_splitconstr.lb[i] > 0 ? variables[df_splitconstr.lb[i]] : 0
                    ub = df_splitconstr.ub[i] < typemax(Int) ? variables[df_splitconstr.ub[i]] : 1
                    con = @build_constraint(sum(variables[i] for i in df_splitconstr.leaves[i]) <= ub - lb)
                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)
                end

                num_lazy += length(arr_conindex)
            
            else
                ### Fractional Solution

                ##### Before querying `callback_value`, you must call:
                Gurobi.load_callback_variable_primal(cb_data, cb_where)
                for i = 1:size(df_sol)[1]
                    df_sol[i, :fsol] = callback_value(cb_data, variables[i])
                end
                
                arr_conindex = []
                for df in groupby(df_splitconstr, [:treeid])
                    maxvio = 0
                    conindex = 0

                    for r in eachrow(df)
                        lhs = sum(df_sol.fsol[i] for i in r.leaves)
                        lb = r.lb > 0 ? df_sol.fsol[r.lb] : 0
                        ub = r.ub < typemax(Int) ? df_sol.fsol[r.ub] : 1
                        vio = lhs - (ub - lb)
                        
                        if vio > maxvio
                            maxvio = vio
                            conindex = r.index
                        end
                    end

                    if maxvio > min_violation
                        push!(arr_conindex, conindex)
                    end
                end

                ##### Add user cuts
                for i in arr_conindex
                    lb = df_splitconstr.lb[i] > 0 ? variables[df_splitconstr.lb[i]] : 0
                    ub = df_splitconstr.ub[i] < typemax(Int) ? variables[df_splitconstr.ub[i]] : 1
                    con = @build_constraint(sum(variables[i] for i in df_splitconstr.leaves[i]) <= ub - lb)
                    MOI.submit(model, MOI.UserCut(cb_data), con)
                end

                num_user += length(arr_conindex)
            end

            time_cb += time() - t_begin
        end
            
        JuMP.set_optimizer_attribute(model, "LazyConstraints", 1)
        MOI.set(model, Gurobi.CallbackFunction(), my_callback_checkall)
    else
        @debug "callback-else"
        df_sol = data["df_vars"]
        df_sol.sol = zeros(Int, size(df_sol)[1])

        active_leaves_history = [Set() for t = 1:num_trees]
        if method == "activeleaf-init"
            for t = 1:num_trees
                df = df_sol[df_sol.treeid .== t, :]
                df = df[df.pred .>= 0.9 * (maximum(df.pred) - minimum(df.pred)) + minimum(df.pred), :]
                active_leaves_history[t] = Set(df.varindex)
            end
            @info [length(i) for i in active_leaves_history]
        elseif method == "liftall"
            for t = 1:num_trees
                df = df_sol[df_sol.treeid .== t, :]
                active_leaves_history[t] = Set(df.varindex)
            end
            @info [length(i) for i in active_leaves_history]
        end

        function my_callback(cb_data, cb_where::Cint)
            t_begin = time()

            ##### You can select where the callback is run
            if cb_where != Gurobi.GRB_CB_MIPSOL
                return
            end

            num_call += 1

            ##### Before querying `callback_value`, you must call:
            Gurobi.load_callback_variable_primal(cb_data, cb_where)
            for i = 1:size(df_sol)[1]
                df_sol[i, :sol] = convert(Int, round(callback_value(cb_data, variables[i])))
            end        

            ##### Collect active leaves
            df_active_leaves = df_sol[1:num_leaves, :]
            df_active_leaves = df_active_leaves[df_active_leaves.sol .== 1, :]

            ##### Collect violated constraints
            max_depth = typemax(Int)
            max_violated_per_tree = 1
            if method == "maxclq"
                max_violated_per_tree = typemax(Int)
            else
                nothing
            end

            df_lazy = DataFrame()
            if occursin("-unsplit", method)
                for leaf in eachrow(df_active_leaves)
                    n = ptr_forest[leaf.treeid] + leaf.nodeid - 1
                    rec, con_inds = get_rectangle(df_forest, n; get_con_ind=true)
                    for v in keys(rec)
                        ## 
                        nothing
                    end


                end
            else
                for leaf in eachrow(df_active_leaves)
                    t = leaf.treeid
                    leafindex = leaf.varindex
                    n = ptr_forest[t]
                    if leafindex < df_forest[n, "leafbegin"] || leafindex > df_forest[n, "leafend"]
                        @info "error in df_active_leaves"
                    end
                    
                    current_depth = 0
                    num_violated_per_tree = 0
                    while df_forest.lchild[n] > 0 && current_depth <= max_depth && num_violated_per_tree < max_violated_per_tree
                        lch = n + df_forest.lchild[n] - df_forest.nodeid[n]
                        if leafindex <= df_forest[lch, "leafend"]
                            ### leaf is on the left
                            if df_sol.sol[df_forest.varindex[n]] == 0
                                ### violation found: sum(y) <= x
                                num_violated_per_tree += 1                        
                                push!(df_lazy, Dict(:leaves => df_forest.leafbegin[lch]:df_forest.leafend[lch], :lb => 0, :ub => df_forest.varindex[n],
                                                    :con_ind => lch, :treeid => t, :splitvar1 => df_forest.splitvar1[n], :contype => "<=x", :conlevel => num_violated_per_tree), cols=:union)
                            end

                            ### move to the right
                            n = lch
                        else
                            ### leaf is on the right
                            rch = n + df_forest.rchild[n] - df_forest.nodeid[n]
                            if df_sol.sol[df_forest.varindex[n]] == 1
                                ### violation found: sum(y) <= 1 - x
                                num_violated_per_tree += 1                        
                                push!(df_lazy, Dict(:leaves => df_forest.leafbegin[rch]:df_forest.leafend[rch], :lb => df_forest.varindex[n], :ub => typemax(Int),
                                                    :con_ind => rch, :treeid => t, :splitvar1 => df_forest.splitvar1[n], :contype => "<=1-x", :conlevel => num_violated_per_tree), cols=:union)
                            end

                            ### move to the right
                            n = rch
                        end

                        current_depth += 1
                    end
                end
            end

            ##### Refine df_lazy depending on the method
            if nrow(df_lazy) > 0
                ##### Lift x
                if method == "misic"
                    nothing
                else
                    df_lazy = df_lazy_misic2bounded(df_forest, df_lazy)
                end

                ##### Lift y
                if occursin("liftx", method)
                    nothing
                elseif method == "maxclq"
                    df_lazy = df_lazy_misic2bounded(df_forest, df_lazy)

                    df_maxclq = DataFrame()
                    for row in eachrow(df_lazy)
                        if row.conlevel == 1
                            df = DataFrame()
                            push!(df, Dict(:con_ind => row.con_ind, :leaves => row.leaves, :lb => row.lb, :ub => row.ub), cols=:union)

                            for r in eachrow(df_lazy)
                                if r.treeid != row.treeid
                                    lb, ub = get_interval(df_forest, r.con_ind, row.splitvar1)
                                    if lb >= row.lb && ub <= row.ub
                                        push!(df, Dict(:con_ind => r.con_ind, :leaves => r.leaves, :lb => lb, :ub => ub), cols=:union)
                                    end
                                end
                            end

                            df = get_maxclq_cuts(df_forest, df, must_include=[1], maxcuts=1)
                            append!(df_maxclq, df)
                        end
                    end
                elseif method in ["activeleaf", "activeleaf-init", "activeleaf-only", "liftall", "activeleaf-unsplit", "activeleaf-init-unsplit", "activeleaf-only-unsplit", "liftall-unsplit"]
                    df_lazy = df_lazy_misic2bounded(df_forest, df_lazy)
                    
                    df_lazy_leaves = []
                    for r in eachrow(df_lazy)
                        t = r.treeid
                        con_leaves = Array(r.leaves)

                        for leaf in active_leaves_history[t]
                            if !(leaf in con_leaves)
                                lb, ub = get_interval(df_forest, ptr_forest[t] + df_sol.nodeid[leaf] - 1, r.splitvar1)
                                if lb >= r.lb && ub <= r.ub
                                    push!(con_leaves, leaf)
                                    @debug("num_call = $num_call, lifted! $(r.con_ind)")
                                end
                            end
                        end
                        push!(df_lazy_leaves, con_leaves)
                    end

                    df_lazy.leaves = df_lazy_leaves
                else
                    @assert true "should not be here (wrong method name)"
                end

                ##### Log active leaves
                for r in eachrow(df_active_leaves)
                    union!(active_leaves_history[r.treeid], r.varindex)
                end

                ##### Add lazy constraints
                for r in eachrow(df_lazy)
                    leaves = r.leaves
                    if method == "activeleaf-only"
                        leaves = r.leaves[:]
                        intersect!(leaves, active_leaves_history[r.treeid])
                    end

                    if r.lb >= r.ub
                        @info "lb >= ub found"
                        con = @build_constraint(sum(variables[i] for i in leaves) <= 0)
                        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
                    else
                        con = @build_constraint(sum(variables[i] for i in leaves) <= 
                                                (r.ub < typemax(Int) ? variables[r.ub] : 1) - (r.lb > 0 ? variables[r.lb] : 0))
                        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
                    end
                end

                num_lazy += nrow(df_lazy)

                ##### Heuristic solution        
                # status = MOI.submit( model, MOI.HeuristicSolution(cb_data), [x], [floor(Int, x_val)] )
                # @debug "num_call = $num_call, I submitted a heuristic solution, and the status was: $status"
            end


            time_cb += time() - t_begin

        end
            
        JuMP.set_optimizer_attribute(model, "LazyConstraints", 1)
        MOI.set(model, Gurobi.CallbackFunction(), my_callback)
    end

    @debug "save formulate time"
    output[:time_formulate] = time() - t_begin

    ## optimize!(model)
    @debug("Begin optimize!()")
    optimize!(model)

    ##### Get output
    output[:num_constrs] = MOI.get(model, Gurobi.ModelAttribute("NumConstrs"))
    output[:num_vars] = MOI.get(model, Gurobi.ModelAttribute("NumVars"))
    output[:time_solve] = JuMP.solve_time(model)
    output[:num_call] = num_call
    output[:num_lazy] = num_lazy
    output[:num_user] = num_user
    output[:time_cb] = time_cb

    if MOI.get(model, Gurobi.ModelAttribute("isMIP")) == 1
        output[:MIPGap] = MOI.get(model, Gurobi.ModelAttribute("MIPGap"))
        output[:ObjBound] = MOI.get(model, Gurobi.ModelAttribute("ObjBound"))
        output[:MIPNodes] = MOI.get(model, Gurobi.ModelAttribute("NodeCount"))
    else
        output[:MIPGap] = missing
        output[:ObjBound] = missing
        output[:MIPNodes] = missing
    end

    if JuMP.has_values(model)
        output[:objval] = JuMP.objective_value(model)

        df_sol = data["df_vars"]
        if nrow(df_sol) < length(JuMP.all_variables(model))
            num_missing = length(JuMP.all_variables(model)) - nrow(df_sol)
            for i = 1:num_missing 
                append!(df_sol, Dict("varindex" => nrow(df_sol) + 1, "vartype" => "unknown"), cols=:union)
            end
        end

        if occursin("-LP", method)
            df_sol.sol = JuMP.value.(JuMP.all_variables(model))
        else
            df_sol.sol = convert(Array{Int}, round.(JuMP.value.(JuMP.all_variables(model))))
        end
        CSV.write(dir * "df_sol.csv", df_sol)
        output[:check_feasibility] = check_feasibility(data["df_forest"], df_sol[1:length(variables), :])
    else
        output[:objval] = missing
        output[:check_feasibility] = missing
    end

    output 
end


##### Rewrite the code 

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

function add_splitconstr(model, df_splitconstr, lazylevel)
    variables = JuMP.all_variables(model)

    for r in eachrow(df_splitconstr)
        if length(r.leaves) == 0
            nothing
        elseif r.lb < r.ub
            lb = r.lb > 0 ? variables[r.lb] : 0
            ub = r.ub < typemax(Int) ? variables[r.ub] : 1
            c = @constraint(model, sum(variables[j] for j in r.leaves) <= ub - lb)
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
    @assert formulation in ["plain", "misic", "liftx", "liftall"]
    @assert lazylevel in [0, 2, 3]

    model = initialize_model(optimizer_attributes)
    model = add_basic_formulation(model, data)

    if formulation in ["misic", "liftx", "liftall"] 
        df_splitconstr, df_forest = get_df_splitconstr_ver2(formulation, data["df_forest"], data["df_vars"])
        model = add_splitconstr(model, df_splitconstr, lazylevel)
        data["df_splitconstr"] = deepcopy(df_splitconstr)
        data["df_forest"] = deepcopy(df_forest)
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

# function get_model_output(model, data)
#     @assert "df_forest" in keys(data)
#     @assert "df_vars" in keys(data)

#     ### Initialize
#     output = Dict()

#     ### Get output from gurobi report
#     output[:num_constrs] = MOI.get(model, Gurobi.ModelAttribute("NumConstrs"))
#     output[:num_vars] = MOI.get(model, Gurobi.ModelAttribute("NumVars"))
#     output[:num_binary] = get_num_binary_variables(model)
#     output[:time_solve] = JuMP.solve_time(model)
#     # output[:num_call] = num_call
#     # output[:num_lazy] = num_lazy
#     # output[:num_user] = num_user
#     # output[:time_cb] = time_cb

#     if MOI.get(model, Gurobi.ModelAttribute("isMIP")) == 1
#         output[:MIPGap] = MOI.get(model, Gurobi.ModelAttribute("MIPGap"))
#         output[:ObjBound] = MOI.get(model, Gurobi.ModelAttribute("ObjBound"))
#         output[:MIPNodes] = MOI.get(model, Gurobi.ModelAttribute("NodeCount"))
#     else
#         output[:MIPGap] = missing
#         output[:ObjBound] = missing
#         output[:MIPNodes] = missing
#     end

#     ### Get output from solution
#     df_sol = data["df_vars"]
#     if JuMP.has_values(model)
#         output[:objval] = JuMP.objective_value(model)

#         if nrow(df_sol) < length(JuMP.all_variables(model))
#             num_missing = length(JuMP.all_variables(model)) - nrow(df_sol)
#             for i = 1:num_missing 
#                 append!(df_sol, Dict("varindex" => nrow(df_sol) + 1, "vartype" => "unknown"), cols=:union)
#             end
#         end

#         df_sol.sol = JuMP.value.(JuMP.all_variables(model))
#         # df_sol.sol = convert(Array{Int}, round.(JuMP.value.(JuMP.all_variables(model))))
#         output[:check_feasibility] = check_feasibility(data["df_forest"], df_sol)
#     else
#         output[:objval] = missing
#         output[:check_feasibility] = missing
#     end
#     data["df_sol"] = df_sol

#     output, data
# end

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
