using JuMP

function get_model_output(model, data)
    @assert "df_forest" in keys(data)
    @assert "df_vars" in keys(data)

    ### Initialize
    output = Dict()

    ### Get output from gurobi report
    # output[:num_constrs] = MOI.get(model, Gurobi.ModelAttribute("NumConstrs"))
    # output[:num_vars] = MOI.get(model, Gurobi.ModelAttribute("NumVars"))
    # output[:num_binary] = get_num_binary_variables(model)
    # output[:time_solve] = JuMP.solve_time(model)
    # output[:num_call] = num_call
    # output[:num_lazy] = num_lazy
    # output[:num_user] = num_user
    # output[:time_cb] = time_cb

    num_constr = 0
    for constr_type in JuMP.list_of_constraint_types(model)
        num_constr += JuMP.num_constraints(model, constr_type[1], constr_type[2])
    end
    output[:num_constrs] = num_constr
    output[:num_vars] = JuMP.num_variables(model)
    output[:num_binary] = length([v for v in JuMP.all_variables(model) if JuMP.is_binary(v)])
    output[:time_solve] = JuMP.solve_time(model)
    output[:MIPGap] = JuMP.relative_gap(model)
    output[:ObjBound] = JuMP.objective_bound(model)        
    output[:MIPNodes] = JuMP.node_count(model)

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
        # df_sol.sol = convert(Array{Int}, round.(JuMP.value.(JuMP.all_variables(model))))
        output[:check_feasibility] = check_feasibility(data["df_forest"], df_sol)
    else
        output[:objval] = missing
        output[:check_feasibility] = missing
    end
    data["df_sol"] = df_sol

    output, data
end