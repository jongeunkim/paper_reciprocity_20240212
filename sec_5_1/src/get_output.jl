function check_feasibility(df_forest, df_sol)
    df_sol = filter(row -> row.vartype in ["x", "y"], df_sol)
    num_trees = maximum(df_forest.treeid)

    ### Binary check
    EPS = 0.01
    df_nonbinary = filter(row -> row.sol > EPS && row.sol < 1-EPS, df_sol)
    if !isempty(df_nonbinary)
        return false
    end

    ### Compute sol by y 
    df_sol_y = filter(row -> row.vartype == "y", df_sol)
    obj_y = sum([row.pred * row.sol for row in eachrow(df_sol_y)])

    ### Compute sol by x 
    obj_x = 0
    for t = 1:num_trees
        df_tree = df_forest[df_forest.treeid .== t, :]

        n = 1
        while n > 0
            if df_tree.lchild[n] == 0
                break
            end

            if df_sol[df_tree.varindex[n], "sol"] > 1-EPS
                n = df_tree.lchild[n]
            else
                n = df_tree.rchild[n]
            end
        end

        obj_x += df_tree.pred[n]
    end

    if abs(obj_y - obj_x) / abs(obj_y) > 1e-04
        return false
    else
        return true
    end
end

function get_output(model, data)
    @assert "df_forest" in keys(data)
    @assert "df_vars" in keys(data)

    ### Initialize
    output = Dict()

    num_constr = 0
    for constr_type in JuMP.list_of_constraint_types(model)
        num_constr += JuMP.num_constraints(model, constr_type[1], constr_type[2])
    end
    output[:num_constrs] = num_constr
    output[:num_vars] = JuMP.num_variables(model)
    output[:num_binary] = length([v for v in JuMP.all_variables(model) if JuMP.is_binary(v)])
    output[:time_solve] = JuMP.solve_time(model)
    output[:mipgap] = JuMP.relative_gap(model)
    output[:objbound] = JuMP.objective_bound(model)        
    output[:nodecnt] = JuMP.node_count(model)

    ### Get output from solution
    df_sol = data["df_vars"]
    if JuMP.has_values(model)
        output[:objval] = JuMP.objective_value(model)

        if nrow(df_sol) < length(JuMP.all_variables(model))
            num_missing = length(JuMP.all_variables(model)) - nrow(df_sol)
            for i = 1:num_missing 
                append!(df_sol, Dict("varindex" => nrow(df_sol) + 1, "vartype" => "unknown"), cols=:union)
            end
        end

        df_sol.sol = JuMP.value.(JuMP.all_variables(model))
        output[:check_feasibility] = check_feasibility(data["df_forest"], df_sol)
    else
        output[:objval] = missing
        output[:check_feasibility] = missing
    end
    data["df_sol"] = df_sol

    output, data
end
