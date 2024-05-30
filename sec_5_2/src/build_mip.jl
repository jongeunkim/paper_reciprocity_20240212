function treedf2boxes(df)
    """
    Collect all boxes using DFS
    """

    function collect_box(df, nodeid, current_box, boxes, predval)
        if df[nodeid, :lchild] == 0
            push!(boxes, current_box)
            push!(predval, df[nodeid, :predval])
            return boxes, predval           
        end

        # Collect boxes from the left child
        box = deepcopy(current_box)
        lbub = get(box, df[nodeid,:featid], [-Inf,Inf])
        lbub[2] = min(lbub[2], df[nodeid,:featval])
        box[df[nodeid,:featid]] = lbub
        boxes, predval = collect_box(df, df[nodeid,:lchild], box, boxes, predval)

        # Collect boxes from the right child
        box = deepcopy(current_box)
        lbub = get(box, df[nodeid,:featid], [-Inf,Inf])
        lbub[1] = max(lbub[1], df[nodeid,:featval])
        box[df[nodeid,:featid]] = lbub
        boxes, predval = collect_box(df, df[nodeid,:rchild], box, boxes, predval)

        # Return boxes
        return boxes, predval
    end

    boxes, predval = collect_box(df, 1, Dict(), [], [])
    boxes, predval
end

function add_cost_for_each_tree(formulation, model, z, wt, C, capa, df_tree, eps)
    @assert formulation in ["MC", "DCC", "DLog"] "wrong formulation"

    boxes, predval = treedf2boxes(df_tree)
    n_boxes = length(boxes)

    y = nothing
    if formulation in ["MC", "DCC"]
        # Variable y to indicate which box is active
        y = @variable(model, [b=1:n_boxes], Bin)
        @constraint(model, sum(y) == 1)
    else
        # DLog
        y = @variable(model, [b=1:n_boxes])
        @constraint(model, [b=1:n_boxes], y[b] >= 0)
        @constraint(model, sum(y) == 1)

        L = Int(ceil(log2(n_boxes)))
        ylog = @variable(model, [i=1:L], Bin)

        for i = 1:L
            @constraint(model, sum(y[b] for b=1:n_boxes if b % 2^i <= 2^(i-1)-1) <= ylog[i])
            @constraint(model, sum(y[b] for b=1:n_boxes if b % 2^i >= 2^(i-1)) <= 1 - ylog[i])
        end
    end
    
    # Define cost for this tree
    @constraint(model, wt == sum(predval[b] * y[b] for b=1:n_boxes))

    if formulation == "MC"
        # Define variable zy for flows for each box  
        zy = @variable(model, [b=1:n_boxes, c=1:C])
        @constraint(model, [c=1:C], z[c] == sum(zy[:,c]))

        for (b, box) in enumerate(boxes), c=1:C
            lbub = get(box, c, [-Inf, Inf])
            lbub[1] = lbub[1] > 0.0 ? lbub[1] + eps : 0.0
            lbub[2] = lbub[2] < capa[c] ? lbub[2] - eps : capa[c]

            @constraint(model, zy[b,c] >= y[b] * lbub[1])
            @constraint(model, zy[b,c] <= y[b] * lbub[2])
        end
    elseif formulation in ["DCC", "DLog"]
        lam = @variable(model, [b=1:n_boxes, k=1:2^C])
        @constraint(model, [b=1:n_boxes, k=1:2^C], lam[b,k] >= 0)
        @constraint(model, sum(lam) == 1)
        @constraint(model, [b=1:n_boxes], sum(lam[b,:]) <= y[b])

        corners = []
        for (b, box) in enumerate(boxes)
            lbubs = [get(box, c, [-Inf, Inf]) for c = 1:C]
            for c = 1:C
                lbubs[c][1] = lbubs[c][1] > 0.0 ? lbubs[c][1] + eps : 0.0
                lbubs[c][2] = lbubs[c][2] < capa[c] ? lbubs[c][2] - eps : capa[c]
            end
            push!(corners, collect(Iterators.product(lbubs...))[:])
        end
        @constraint(model, [c=1:C], z[c] == sum(corners[b][k][c] * lam[b,k] for b=1:n_boxes, k=1:2^C))
    end

    model
end

function add_cost_for_ensemble(formulation, model, z, wt, capa, arr_df_tree, eps)
    # Variable x
    df_vars = vcat(arr_df_tree...)
    df_vars = DF.filter!(row -> row.lchild > 0, df_vars)
    DF.select!(df_vars, [:featid, :featval])
    df_vars = unique(df_vars)
    DF.sort!(df_vars, [:featid, :featval])
    num_xvars = DF.nrow(df_vars)
    df_vars[!,:xidx] = 1:num_xvars

    x = @variable(model, [1:num_xvars], Bin)
    for i = 1:num_xvars-1
        if df_vars[i,:featid] == df_vars[i+1,:featid]
            @constraint(model, x[i] <= x[i+1])
        end
    end

    # # x and z constraint
    # for i = 1:num_xvars
    #     featid = df_vars[i,:featid]
    #     featval = df_vars[i,:featval]

    #     lb = featval + eps
    #     @constraint(model, z[featid] >= lb * (1 - x[i]))

    #     ub = featval - eps
    #     @constraint(model, z[featid] <= ub * x[i] + capa[featid] * (1 - x[i]))
    # end

    # x and z constraint (new)
    for df in DF.groupby(df_vars, [:featid])
        featid = df.featid[1]
        xvars = vcat([0.0], [x[row.xidx] for row in eachrow(df)], [1.0])
        vals_lb = vcat([0.0], df.featval .+ eps, [capa[featid]])
        vals_ub = vcat([0.0], df.featval .- eps, [capa[featid]])
        n = length(xvars)

        # println(df)
        # println("featid = $featid")
        # println("xvars = $xvars")
        # println("vals_lb = $vals_lb")
        # println("vals_ub = $vals_ub")

        @constraint(model, z[featid] >= sum(vals_lb[i]   * (xvars[i+1] - xvars[i]) for i = 1:n-1))
        @constraint(model, z[featid] <= sum(vals_ub[i+1] * (xvars[i+1] - xvars[i]) for i = 1:n-1))
    end

    # For each tree
    xinfo2xidx = Dict((df_vars[i,:featid], df_vars[i,:featval]) => i for i = 1:num_xvars)
    for (t, df_tree) in enumerate(arr_df_tree)
        num_nodes = DF.nrow(df_tree)
        
        # Declare variable y for each node
        y = @variable(model, [1:num_nodes])
        @constraint(model, [i=1:num_nodes], y[i] >= 0)
        @constraint(model, y[1] == 1)
        for row in eachrow(df_tree)
            if row.lchild > 0
                @constraint(model, y[row.nodeid] == y[row.lchild] + y[row.rchild])
            end
        end

        # y and cost
        @constraint(model, wt[t] == sum(row.predval * y[row.nodeid] for row in eachrow(df_tree) if row.lchild == 0))

        # y <= x constraints
        if formulation in ["TreeO", "TreeM"]
            for row in eachrow(df_tree)
                if row.lchild > 0
                    xvarlb, xvarub = nothing, nothing
                    
                    if formulation == "TreeO"
                        n = row.nodeid
                        p = row.parent
                        while p > 0
                            if df_tree[p, :featid] == row.featid
                                # update xvarlb, xvarub
                                if df_tree[p, :lchild] == n
                                    if isnothing(xvarub)
                                        xvarub = x[xinfo2xidx[df_tree[p,:featid], df_tree[p,:featval]]]
                                    end
                                else
                                    if isnothing(xvarlb)
                                        xvarlb = x[xinfo2xidx[df_tree[p,:featid], df_tree[p,:featval]]]
                                    end
                                end

                                if !isnothing(xvarlb) && !isnothing(xvarub)
                                    break 
                                end
                            end

                            # Go up
                            n = p
                            p = df_tree[n, :parent]
                        end
                    end

                    if isnothing(xvarlb)
                        xvarlb = 0
                    end

                    if isnothing(xvarub)
                        xvarub = 1
                    end

                    xvar = x[xinfo2xidx[row.featid, row.featval]]
                    @constraint(model, y[row.lchild] <= xvar - xvarlb)
                    @constraint(model, y[row.rchild] <= xvarub - xvar)
                end
            end
        else
            # TreeCvx
            boxes, predval = treedf2boxes(df_tree)
            bidx2yidx = [row.nodeid for row in eachrow(df_tree) if row.lchild == 0]

            for featid in unique(df_vars.featid)
                intervals = [get(box, featid, [-Inf, Inf]) for box in boxes]
                intervals = vcat(transpose.(intervals)...)

                if length(unique(intervals)) <= 2
                    continue
                end
                # println(intervals)

                list_lb_ub_indices = generate_bndconstr(intervals)

                for (lb, ub, indices) in list_lb_ub_indices
                    if lb <= -Inf && ub >= Inf
                        continue
                    end
                    # println("$lb $ub $indices")
                    
                    xvarlb = lb <= -Inf ? 0 : x[xinfo2xidx[featid, lb]]
                    xvarub = ub >= Inf ? 1 : x[xinfo2xidx[featid, ub]]

                    @constraint(model, sum(y[bidx2yidx[i]] for i in indices) <= xvarub - xvarlb)
                end
            end
        end
    end

    model
end

function build_mip(formulation, M, N, C, supply, demand, capacity, cost_forest)
    @assert formulation in ["MC", "DCC", "DLog", "TreeM", "TreeO", "TreeCvx"]

    ########## Basic Model with z and rho ##########
    model = direct_model(Gurobi.Optimizer(GRB_ENV))

    # Flow Variable z for each arc and for each commodity
    @variable(model, 0 <= z[i=1:M, j=1:N, c=1:C] <= capacity[i,j,c])

    # Cost Variable w for each arc
    @variable(model, 0 <= w[i=1:M, j=1:N])

    # Constr: Supply and Demand
    @constraint(model, [i=1:M, c=1:C], sum(z[i,j,c] for j=1:N) <= supply[i,c]) 
    @constraint(model, [j=1:N, c=1:C], sum(z[i,j,c] for i=1:M) == demand[j,c]) 

    # Obj:
    @objective(model, Min, sum(w[i,j] for i=1:M, j=1:N))
    ########## Basic Model with z and rho (end) ##########

    @variable(model, 0 <= wt[i=1:M, j=1:N, t=1:length(cost_forest[i,j])])
    @constraint(model, [i=1:M, j=1:N], w[i,j] == sum(wt[i,j,t] for t=1:length(cost_forest[i,j])) / length(cost_forest[i,j]))
    eps = 0.00 #1e-02
    if formulation in ["MC", "DCC", "DLog"]
        for i=1:M, j=1:N, (t,df_tree) in enumerate(cost_forest[i,j])
            model = add_cost_for_each_tree(formulation, model, z[i,j,:], wt[i,j,t], C, capacity[i,j,:], df_tree, eps)
        end
    elseif formulation in ["TreeM", "TreeO", "TreeCvx"]
        for i=1:M, j=1:N
            model = add_cost_for_ensemble(formulation, model, z[i,j,:], [wt[i,j,t] for t=1:length(cost_forest[i,j])], 
                capacity[i,j,:], cost_forest[i,j], eps)
        end
    else
        @assert false "wrong formulation, $formulation"
    end


    model, z, w, wt
end