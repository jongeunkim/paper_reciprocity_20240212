using JuMP, Glob

include("utils.jl")

function get_extvars_optbox(optimal_box, df_leaves, df_forest)
    """
    optimal_box: Dict(indvar => [lb, ub])
    """

    df_extvars_optbox = DataFrame()
    for row in eachrow(df_leaves)
        rec = get_rectangle_by_treeid_node_id(df_forest, row.treeid, row.nodeid)

        if !isConflict_tworec(optimal_box, rec)
            for (k,v) in optimal_box
                lb = v[1]
                ub = v[2]

                ### Left
                if k in keys(rec) && rec[k][1] < lb
                    extvar = Dict("treeid"=>row.treeid, "nodeid"=>row.nodeid, "leafvarindex"=>row.varindex, "splitvar1"=>k, "splitval1_lb"=>rec[k][1], "splitval1_ub"=>lb, "exttype"=>'L')
                    append!(df_extvars_optbox, extvar)
                end

                ### Right
                if k in keys(rec) && rec[k][2] > ub
                    extvar = Dict("treeid"=>row.treeid, "nodeid"=>row.nodeid, "leafvarindex"=>row.varindex, "splitvar1"=>k, "splitval1_lb"=>ub, "splitval1_ub"=>rec[k][2], "exttype"=>'R')
                    append!(df_extvars_optbox, extvar)
                end

                ### Middle
                if k in keys(rec) 
                    extvar = Dict("treeid"=>row.treeid, "nodeid"=>row.nodeid, "leafvarindex"=>row.varindex, "splitvar1"=>k, "splitval1_lb"=>max(lb, rec[k][1]), "splitval1_ub"=>min(ub, rec[k][2]), "exttype"=>'M')
                    append!(df_extvars_optbox, extvar)
                end
            end
        end
    end

    df_extvars_optbox.vartype = ["ext" for i=1:nrow(df_extvars_optbox)]
    header = ["treeid", "nodeid", "leafvarindex", "splitvar1", "splitval1_lb", "splitval1_ub", "exttype", "vartype"]
    select!(df_extvars_optbox, header)
    # CSV.write("text.csv", df_extvars_optbox)

    df_extvars_optbox
end

function merge_df_vars_df_extvars_optbox(df_vars, df_extvars_optbox)
    append!(df_vars, df_extvars_optbox, cols=:union)
    df_vars
end

function add_variables_for_df_extvars_optbox(model, df_vars)
    df_extvars = filter(row -> row.vartype == "ext", df_vars)
    @variable(model, extvars[1:nrow(df_extvars)], lower_bound=0, upper_bound=1)
    model
end

function add_y_constraint_for_df_extvars_optbox(model, df_vars)
    variables = JuMP.all_variables(model)

    df_extvars = filter(row -> row.vartype == "ext", df_vars)
    for df in groupby(df_extvars, ["treeid", "nodeid"])
        @constraint(model, variables[df.leafvarindex[1]] == sum(variables[i] for i in df.varindex))
    end

    model
end

function lift_df_extvars_optbox_to_df_splitconstr(df_splitconstr, df_vars)

    df_extvars_optbox = filter(row -> row.vartype=="ext", df_vars)

    for df_constr in groupby(df_splitconstr, ["treeid", "splitvar1"])
        treeid = df_constr.treeid[1]
        splitvar1 = df_constr.splitvar1[1]

        df_var = filter(row -> row.treeid == treeid && row.splitvar1 == splitvar1, df_extvars_optbox)
        
        for constr in eachrow(df_constr)
            df_var_selected = filter(row -> row.splitval1_lb >= constr.lb && row.splitval1_ub <= constr.ub && !(row.leafvarindex in constr.leaves), df_var)

            # if nrow(df_var_selected) > 0
            #     println(treeid, ", ", splitvar1, ", ", constr.index) 

            # end

            union!(constr.leaves, Set(df_var_selected.varindex))
        end
    end

    df_splitconstr
end

# dir = "debug/"
# foreach(rm, Glob.glob(dir*"*"))

# dataset = "Redwine_T5000_GBT"
# num_trees = 50
# method = "liftall"
# dir_LP = "result_1201/$(dataset)_$(num_trees)_0_$(method)-LP_3600_1/"
# dir_IP = "result_1201/$(dataset)_$(num_trees)_0_$(method)_3600_1/"

# df_splits = DataFrame(CSV.File(dir_IP * "df_splits.csv"))
# num_splitvars = maximum(df_splits.splitvar1)

# df_sol = DataFrame(CSV.File(dir_IP * "df_sol.csv"))
# df_forest = DataFrame(CSV.File(dir_IP * "df_forest.csv"))
# df_splits = DataFrame(CSV.File(dir_IP * "df_splits.csv"))
# df_leaves = DataFrame(CSV.File(dir_IP * "df_leaves.csv"))
# optbox = get_optimal_box(df_sol, df_forest, df_splits, splitval_type="default")

# # get_extvars_optbox(optimal_box, df_leaves, df_forest)



# ensembleopt(dir, "dataset/$(dataset)/", Array(1:num_trees), "liftall-optbox", 1200; loglevel="Info", preprocessing_level=1, optbox=optbox)