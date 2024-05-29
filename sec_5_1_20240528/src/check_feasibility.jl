using CSV, DataFrames

function check_feasibility_old(df_forest, df_sol)
    EPS = 0.01
    num_trees = maximum(df_forest.treeid)
    # old
    # df_sol = df_sol[1:maximum(df_forest.varindex), :]
    df_sol = filter(row -> row.vartype in ["x", "y"], df_sol)
    df_nonbinary = filter(row -> row.sol > EPS && row.sol < 1-EPS, df_sol)
    if !isempty(df_nonbinary)
        return false
    end

    # old
    # df_ysol = filter(row -> row.treeid > 0 && row.sol > 0.5 && row.splitvar1 == 0, df_sol)
    
    # if size(df_ysol)[1] != num_trees
    #     return false
    # elseif length(unique(df_ysol.treeid)) != num_trees
    #     return false
    # end

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

        df = filter(row -> row.vartype in ["y"] && row.treeid == t && (row.sol > 1-EPS || row.pred <= 0.0), df_sol)
        if !(n in df.nodeid)
            @info "feasibility failed at tree $t node $n"
            return false
        end    

        # leafids = df_ysol[t, "nodeid"]
        # if n != leafid
        #     @info "feasibility failed at tree $t node $n"
        #     return false
        # end    
    end

    return true
end

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

        # df = filter(row -> row.vartype in ["y"] && row.treeid == t && (row.sol > 1-EPS || row.pred <= 0.0), df_sol)
        # if !(n in df.nodeid)
        #     @info "feasibility failed at tree $t node $n"
        #     return false
        # end    

        # leafids = df_ysol[t, "nodeid"]
        # if n != leafid
        #     @info "feasibility failed at tree $t node $n"
        #     return false
        # end    
    end

    if abs(obj_y - obj_x) / abs(obj_y) > 1e-04
        return false
    else
        return true
    end
end

# df_forest = DataFrame(CSV.File("result_MIP_1111/Redwine_T5000_200_0_liftall_3600/df_forest.csv"))
# df_sol = DataFrame(CSV.File("result_MIP_1111/Redwine_T5000_200_0_liftall_3600/df_sol.csv"))
# check_feasibility(df_forest, df_sol)