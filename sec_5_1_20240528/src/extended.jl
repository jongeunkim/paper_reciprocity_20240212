include("read_data.jl")
include("utils.jl")

using JuMP, Gurobi, Random

# const GRB_ENV = Gurobi.Env()

function get_df_extvars(ptr_forest, df_forest, df_leaves, df_splits)
    num_trees = maximum(df_forest.treeid)
    num_splitvars = maximum(df_splits.splitvar1)

    ## number of intervals for each split variable
    varindex_first = [df.varindex[1] for df in groupby(df_splits, [:splitvar1])]
    varindex_last = [df.varindex[end] for df in groupby(df_splits, [:splitvar1])]


    df_extvars = DataFrame()

    # rectangles = []
    # num_extendedvars = []
    for r in eachrow(df_leaves)
        forest_index = get_forest_index(ptr_forest, r.treeid, r.nodeid)
        rec = get_rectangle(df_forest, forest_index)

        num_evars = 0
        df = DataFrame()
        newvar = Dict(:treeid=>r.treeid, :nodeid=>r.nodeid, :leafindex=>r.varindex, :rec=>rec)
        for v = 1:num_splitvars
            newvar[:splitvar1] = v

            if v in keys(rec)
                if rec[v][1] >= rec[v][2]
                    df = DataFrame()
                    @info "error"
                    break
                end

                if rec[v][1] == 0
                    newvar[:interval_index_lb] = 0
                    newvar[:interval_index_ub] = varindex_first[v]
                    append!(df, Dict(newvar), cols=:union)
                    for i = varindex_first[v]:(rec[v][2]-1)
                        newvar[:interval_index_lb] = i
                        newvar[:interval_index_ub] = i+1
                        append!(df, Dict(newvar), cols=:union)
                    end
                elseif rec[v][2] == typemax(Int)
                    for i = rec[v][1]:(varindex_last[v]-1)
                        newvar[:interval_index_lb] = i
                        newvar[:interval_index_ub] = i+1
                        append!(df, Dict(newvar), cols=:union)
                    end
                    newvar[:interval_index_lb] = varindex_last[v]
                    newvar[:interval_index_ub] = typemax(Int)
                    append!(df, Dict(newvar), cols=:union)
                else
                    for i = rec[v][1]:(rec[v][2]-1)
                        newvar[:interval_index_lb] = i
                        newvar[:interval_index_ub] = i+1
                        append!(df, Dict(newvar), cols=:union)
                    end
                end
            else
                # nothing

                newvar[:interval_index_lb] = 0
                newvar[:interval_index_ub] = varindex_first[v]
                append!(df, Dict(newvar), cols=:union)
                for i = varindex_first[v]:(varindex_last[v]-1)
                    newvar[:interval_index_lb] = i
                    newvar[:interval_index_ub] = i+1
                    append!(df, Dict(newvar), cols=:union)
                end
                newvar[:interval_index_lb] = varindex_last[v]
                newvar[:interval_index_ub] = typemax(Int)
                append!(df, Dict(newvar), cols=:union)
            end
        end

        append!(df_extvars, df, cols=:union)
    end

    select!(df_extvars, [:treeid, :nodeid, :splitvar1, :leafindex, :interval_index_lb, :interval_index_ub, :rec])

    df_extvars.varindex = 1:nrow(df_extvars)

    df_extvars
end

function add_constraints_for_extvars(model, extvars, variables, df_extvars)
    @info "add_constraints_for_extvars()"

    ##### Constraints between a leaf and extvars 
    for df in groupby(df_extvars, [:leafindex, :splitvar1])
        @constraint(model, variables[df.leafindex[1]] == sum(extvars[i] for i in df.varindex))
    end

    ##### Constraints between an adjacent interval and extvars
    for df in groupby(df_extvars, [:treeid, :interval_index_lb, :interval_index_ub])
        lb = df.interval_index_lb[1] > 0 ? variables[df.interval_index_lb[1]] : 0
        ub = df.interval_index_ub[1] < typemax(Int) ? variables[df.interval_index_ub[1]] : 1
        @constraint(model, sum(extvars[k] for k in df.varindex) == ub - lb)
    end

    return model
end

function main()
    df_info = DataFrame()
    # for treedir in ["dataset/Redwine_T5000_GBT/", "dataset/Concrete_T5000_GBT/", "dataset/Redwine_T5000/", "dataset/Concrete_T5000/"], num_trees = 100:100:1000
    # for treedir in ["dataset/Redwine_T5000_GBT/", "dataset/Concrete_T5000_GBT/", "dataset/Redwine_T5000/", "dataset/Concrete_T5000/"], num_trees = 100:100:1000
    for treedir in ["dataset/Redwine_T5000_GBT/"], num_trees = 10:10:10
    
        treeids = Array(1:num_trees)
        data = read_data([treedir], [treeids])

        df_extvars = get_df_extvars(data["ptr_forest"], data["df_forest"], data["df_leaves"], data["df_splits"])


        result = Dict()
        result["dataset"] = treedir
        result["num_extendedvars"] = nrow(df_extvars)
        result["num_positioningvars"] = nrow(data["df_splits"])
        result["num_leafvars"] = nrow(data["df_leaves"])
        result["num_trees"] = num_trees

        append!(df_info, result, cols=:union)

        df_extvars = df_extvars[df_extvars.treeid .== 1, :]
        CSV.write("debug/df_extvars.csv", df_extvars)

    end

    df_info = select(df_info, ["dataset", "num_trees", "num_positioningvars", "num_leafvars", "num_extendedvars"])
    CSV.write("debug/df_info.csv", df_info)

end

function test_gurobi()
    num_vars = Int(1e+9)

    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV))
    JuMP.set_optimizer_attributes(model, "LogFile" => "debug/gurobi.log")
    
    @variable(model, variables[1:num_vars], lower_bound=0, upper_bound=1)
    # JuMP.set_binary.(variables[num_leaves+1:end])

    # Set objective
    objcoeff = Random.rand(num_vars)
    @objective(model, Max, sum(objcoeff[i] * variables[i] for i=1:num_vars))

    for i = 1:100
        subset = Random.randperm(num_vars)[1:100]
        coeff = Random.rand(100)
        @constraint(model, sum(coeff[i] * variables[subset[i]] for i=1:100) <= 1)
    end

    optimize!(model)
end

# test_gurobi()


# main()


