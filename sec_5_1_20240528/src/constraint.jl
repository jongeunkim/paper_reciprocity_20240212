using DataFrames, DataStructures, JuMP

function get_df_splitconstr(method, df_forest)
    @info "get_df_splitconstr(method=$method)"
    df_splitconstr = DataFrame()
    arr_zero_leaves = Set()

    if occursin("misic", method)
        for r in eachrow(df_forest)
            if r.nodeid > 1
                p = r.index - (r.nodeid - r.parent)

                constr = Dict()
                if df_forest.lchild[p] == r.nodeid
                    constr = Dict(:forestindex=>Set([r.index]), :treeid=>r.treeid, :splitvar1=>df_forest.splitvar1[p], :leaves=>SortedSet(r.leafbegin:r.leafend), :lb=>0, :ub=>df_forest.varindex[p])
                else
                    constr = Dict(:forestindex=>Set([r.index]), :treeid=>r.treeid, :splitvar1=>df_forest.splitvar1[p], :leaves=>SortedSet(r.leafbegin:r.leafend), :lb=>df_forest.varindex[p], :ub=>typemax(Int))
                end

                append!(df_splitconstr, constr)
            end
        end
    elseif occursin("liftx", method)
        df_splitconstr, df_forest_new, arr_zero_leaves = get_df_splitconstr("misic", df_forest)
        @info "liftx (1)"
        for row in eachrow(df_splitconstr)
            row.lb, row.ub = get_interval(df_forest, collect(row.forestindex)[1], row.splitvar1)
        end
        @info "liftx (2)"

    elseif occursin("liftall", method)
        df_splitconstr, df_forest_new, arr_zero_leaves = get_df_splitconstr("misic", df_forest)

        df_liftall = DataFrame()
        for _df1 in groupby(df_splitconstr, [:treeid, :splitvar1])
            # isdebug = false
            # if _df1.treeid[1] == 6 && _df1.splitvar1[1] == 6
            #     isdebug = true
            # end

            treeid = _df1.treeid[1] 
            splitvar1 = _df1.splitvar1[1]
            df = DataFrame()

            # if isdebug
            #     CSV.write("debug_LP_1116_2/debug1.csv", _df1)
            # end

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

            # if isdebug
            #     CSV.write("debug_LP_1116_2/debug2.csv", df)
            # end

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

            # if isdebug
            #     CSV.write("debug_LP_1116_2/debug3.csv", df)
            # end

            ##### Create constraints with bounded intervals
            df_left = df[df.lb .== 0, :]
            df_right = df[df.ub .== typemax(Int), :]
            for l in eachrow(df_left), r in eachrow(df_right)
                leaves_new = intersect(l.leaves, r.leaves)

                forestindex_new = SortedSet()
                # if l.leaves == leaves_new
                #     union!(forestindex_new, l.forestindex)
                # end
                # if r.leaves == leaves_new
                #     union!(forestindex_new, r.forestindex)
                # end
                for row in eachrow(_df1)
                    if issubset(row.leaves, leaves_new)
                        union!(forestindex_new, row.forestindex)
                    end
                end
                
                if !isempty(leaves_new)
                    if r.lb < l.ub && !isempty(leaves_new)
                        bounded_constr = Dict(:forestindex=>forestindex_new, :leaves=>leaves_new, :lb=>r.lb, :ub=>l.ub)
                        append!(df, bounded_constr, cols=:union)
                    else
                        ##### all leaves to be zero
                        union!(arr_zero_leaves, leaves_new)
                    end
                end
            end

            # if isdebug
            #     CSV.write("debug_LP_1116_2/debug4.csv", df)
            # end

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

            # if isdebug
            #     CSV.write("debug_LP_1116_2/debug5.csv", df)
            # end

            df.treeid = treeid * ones(Int, nrow(df))
            df.splitvar1 = splitvar1 * ones(Int, nrow(df))

            append!(df_liftall, df, cols=:union)
        end

        df_splitconstr = df_liftall
    end

    df_splitconstr.index = 1:nrow(df_splitconstr)
    df_forest_new = df_forest[:,:]
    df_forest_new.conindex = zeros(Int, nrow(df_forest_new))
    for r in eachrow(df_splitconstr)
        for i in r.forestindex
            df_forest_new.conindex[i] = r.index
        end
    end

    df_splitconstr, df_forest_new, arr_zero_leaves
end

function add_constraints(model, variables, df_constr)
    @info "add_constraints()"

    for r in eachrow(df_constr)
        lb = r.lb > 0 ? variables[r.lb] : 0
        ub = r.ub < typemax(Int) ? variables[r.ub] : 1
        @constraint(model, sum(variables[j] for j in r.leaves) <= ub - lb)
    end
    
    return model
end