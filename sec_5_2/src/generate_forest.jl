function compute_cost(flow, capa)
    denom = 0.25 * sum(capa)
    p1 = max(0.0, 0.25 * sum(capa) - sum(flow)) / denom
    p2 = 5 * max(0.0, sum(flow) - 0.75 * sum(capa)) / denom
    return (1 + p1 + p2) * LinearAlgebra.norm(flow, 2)
end

function generate_cost_forest(capa, max_depth, n_trees, seed)
    """
    Input:
        C: # of commodities
        capa: capacity cost vector 
        b: unit size vector 
        B: shipping size 
    Output:
        df_forest
    """

    features = collect(Iterators.product([0:c for c in capa]...))[:]
    features = vcat(reduce.(hcat, features)...) # array of tuple to 2D-array
    features = convert(Array{Float64}, features)
    # function compute_cost(flow)
    #     discount = sum(flow .> 0) <= 1 ? 0.0 : 0.25
    #     return (1 - discount) * sqrt(B) * ceil(sum(b .* flow) / B)
    # end

    labels = [compute_cost(flow, capa) for flow in eachrow(features)]    

    n_subfeatures=size(features, 2); partial_sampling=0.7
    min_samples_leaf=5; min_samples_split=2; min_purity_increase=0.0; seed=seed

    model = DT.build_forest(labels, features,
                        n_subfeatures,
                        n_trees,
                        partial_sampling,
                        max_depth,
                        min_samples_leaf,
                        min_samples_split,
                        min_purity_increase;
                        rng = seed)

    function tree2dataframe(tree)
        """
        df = nodeid | lchild | rchild | featid | featval | predval
        """
    
        df = DF.DataFrame(nodeid = Int[], lchild = Int[], rchild = Int[], parent = Int[], featid = Int[], featval = Float64[], predval = Float64[])
    
        function update_info(tree, df, parent)
            nodeid = size(df, 1) + 1
            row = Dict(:nodeid => nodeid, :lchild => 0, :rchild => 0, :parent => parent, :featid => 0, :featval => -9999.0, :predval => -9999.0)
    
            if DT.is_leaf(tree)
                row[:predval] = tree.majority
                push!(df, row)
                return df, nodeid
            end
    
            row[:featid] = tree.featid
            row[:featval] = tree.featval
            push!(df, row)
    
            df, lchild = update_info(tree.left, df, nodeid)
            df[nodeid, :lchild] = lchild
    
            df, rchild = update_info(tree.right, df, nodeid)
            df[nodeid, :rchild] = rchild
    
            return df, nodeid
        end
    
        df, _ = update_info(tree, df, 0)
        df
    end

    arr_df_tree = []
    for (t,tree) in enumerate(model.trees)
        df_tree = tree2dataframe(tree)
        df_tree[:,:treeid] .= t #* ones(nrow(df_tree))
        push!(arr_df_tree, df_tree)
    end
    # df_forest = vcat(arr_df_tree...)
    # df_forest = DF.select(df_forest, [:treeid, :nodeid, :lchild, :rchild, :featid, :featval, :predval])

    return arr_df_tree, model
end