using JuMP

function get_binary_variables(model)
    [v for v in JuMP.all_variables(model) if JuMP.is_binary(v)]
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

function get_num_binary_variables(model)
    length(get_binary_variables(model))
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

function set_branch_priority(model, df_vars)
    """
    Change Gurobi's Variable Attribute `BranchPriority`.

    The description of BranchPriority extracted from Gurobi Reference Manual:
        Variable branching priority. The value of this attribute is used as the primary criterion for selecting a fractional variable 
        for branching during the MIP search. Variables with larger values always take priority over those with smaller values. 
        Ties are broken using the standard branch variable selection criteria. The default variable branch priority value is zero.

        Note that deleting variables from your model will cause several attributes to be discarded (variable hints and branch priorities). 
        If you'd like them to persist, your program will need to repopulate them after deleting the variables and making a subsequent model update call.

        Extracted from https://www.gurobi.com/documentation/9.1/refman/branchpriority.html.
    """
    
    @assert "varindex" in names(df_vars)
    @assert "branchpriority" in names(df_vars)

    variables = all_variables(model)
    for row in eachrow(df_vars)
        MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), variables[row.varindex], row.branchpriority)
    end

    model
end

function fix_xvar_by_box(model, data, box)
    # @info "fix_xvar_by_box(model, data, box)"

    df_vars = data["df_vars"]
    variables = JuMP.all_variables(model)

    for row in eachrow(df_vars)
        # println(row.vartype)
        if row.vartype == "x"
            var = row.splitvar
            val = row.splitval
            xvar = variables[row.varindex]

            if val <= get(box, var, [-Inf, Inf])[1]
                ### Fix to 0
                # println("fix 0 $var $val")
                JuMP.fix(xvar, 0.0, force=true)
            elseif val >= get(box, var, [-Inf, Inf])[2]
                ### Fix to 1
                # println("fix 1 $var $val")
                JuMP.fix(xvar, 1.0, force=true)
            else
                ### Unfix to [0,1]
                # println("unfix $var $val")
                if JuMP.is_fixed(xvar)
                    JuMP.unfix(xvar)
                    JuMP.set_upper_bound(xvar, 1)
                    JuMP.set_lower_bound(xvar, 0)
                end
            end

        end
    end

    model
end