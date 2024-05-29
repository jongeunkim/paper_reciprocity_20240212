# Julia v1.4.2
# JuMP v.0.21.3
# Gurobi v0.8.1 (optimizer v9.0.3)

using Formatting, JuMP, MathOptInterface
using Gurobi # CPLEX does not have lazy constraint option in Julia
const MOI = MathOptInterface

include("utils.jl") # collection of utility functions such as read_forest()

# column name-number pairs of the forest table
const (NODEID, LCHILD, RCHILD, SPLITVAR, SPLITVAL, NODETYPE, YVAL, TREEID, PARENT, DEPTH, MAXYVAL) = (1,2,3,4,5,6,7,8,9,10,11)

const GRB_ENV = Gurobi.Env()

function solve_MIP_xy(directory::String, num_trees::Int; seed=0, formulation="unbounded", lazylevel=2, relax=false,
                      TimeLimit=Inf, MIPGap=1e-04, logfile="optimizer.log", return_solution=false, verbose=1,
                      initsol=false, fixsol=false, solution_X=nothing) 
    """
    Function solve_MIP_xy() solves a tree ensemble optimization problem using xy-formulation.
    
    INPUT:
    directory = the path of the directory that contains random forest data files
    num_tree = the number of trees
    seed: if seed == 0, then read trees in order (1,2,3,...). if seed > 0, then randomly choose trees to read
    formulation = unbounded, bounded, unbounded_min, bounded_min, unbounded_max, bounded_max
    lazylevel = 0, 1, 2, 3 *please refer http://www.gurobi.com/documentation/8.1/refman/lazy.html
    relax = true (LP relaxation), false (MIP)
    TimeLimit
    MIPGap
    logfile = the path where to save the log of the optimizer 
    return_solution = false (return output), true (return output, sol_X, sol_Y)
    verbose = 0 (nothing), 1 (general info), 2 (debugging purpose)

    * option 1) provide an initial solution: initsol = true, solution_X = solution_X
    * option 2) fix (a subset of) X variables: fixsol = true, solution_X = solution_X
    """
    
    t_begin_formulate = time()
    t_begin = time()
    output = Dict()
    
    dataname = split(directory, "/", keepempty=false)[end]
    lazylevel = relax ? 2 : lazylevel

    print(verbose > 0 ? "\n" : "")
    print(verbose > 0 ? "solve_MIP_xy()\n" : "")
    print(verbose > 0 ? format("directory:   {}\n", directory) : "")
    print(verbose > 0 ? format("dataname:    {}\n", dataname) : "")
    print(verbose > 0 ? format("num_trees:   {}\n", num_trees) : "")
    print(verbose > 0 ? format("seed:        {}\n", seed) : "")
    print(verbose > 0 ? format("formulation: {}\n", formulation) : "")
    print(verbose > 0 ? format("lazylevel:   {}\n", lazylevel) : "")
    print(verbose > 0 ? format("relax:       {}\n", relax) : "")
    print(verbose > 0 ? format("TimeLimit:   {}\n", TimeLimit) : "")
    print(verbose > 0 ? format("MIPGap:      {}\n", MIPGap) : "")
    print(verbose > 0 ? format("logfile:     {}\n", logfile) : "")
    print(verbose > 0 ? format("verbose:     {}\n", verbose) : "")
    
    # Read forest from txt files
    forest = read_forest(directory, num_trees, seed=seed, verbose=verbose)
    print(verbose > 0 ? format("Read data (duration: {:.2f} sec)\n", time() - t_begin) : "")
    t_begin = time()
    
    # Convert forest to a single array [num_nodes, 8] (8th column is treeId)
    forest, ptr_forest = convert_forest_to_single_array(forest)
    
    # Add parent and depth columns [num_nodes, 10] (9th, 10th and 11th are parent, depth and maxyval)
    forest = add_columns_to_forest(forest, ptr_forest; add_parent=true, add_depth=true, add_maxyVal=true)

    # Create leaf_nodes and split_nodes from forest
    split_nodes = forest[forest[:, LCHILD] .> 0, :]
    leaf_nodes = forest[forest[:, LCHILD] .== 0, :]
    num_split_nodes = size(split_nodes, 1)
    num_leaf_nodes = size(leaf_nodes, 1)
    ptr_leaf_nodes = get_pointers(leaf_nodes, TREEID)
    print(verbose > 1 ? format("ptr_leaf_nodes: {}\n", ptr_leaf_nodes) : "")
    
    # Get splitVars, splits, ptr_splits
    # splitVars: list of unique splitVars
    # splits: list of unique pair of splitVar and splitVal
    # ptr_splits: pointer for splits
    splitVars, splits, ptr_splits = get_splits(split_nodes, verbose=verbose)
    print(verbose > 0 ? format("Preprocessed data (duration: {:.2f} sec)\n", time() - t_begin) : "")
    t_begin = time()
    
    # Initialize MIP model
    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV))
    JuMP.set_optimizer_attributes(model, "LogFile" => logfile, "LogToConsole" => 0, "Presolve" => 0, "TimeLimit" => TimeLimit, "MIPGap" => MIPGap, "LazyConstraints" => min(lazylevel, 1))


    print(verbose > 0 ? format("Initialized a MIP model (duration: {:.2f} sec)\n", time() - t_begin) : "")
    t_begin = time()
    
    # Add variables (x:split, y:leaf)
    @variable(model, y[i=1:num_trees, j=1:(ptr_leaf_nodes[i+1]-ptr_leaf_nodes[i])], lower_bound=0, upper_bound=1)
    @variable(model, x[i=1:length(splitVars), j=1:(ptr_splits[i+1]-ptr_splits[i])], Bin)
    variables = all_variables(model)
    variables_y = variables[1:num_leaf_nodes]
    print(verbose > 0 ? format("Added variables (duration: {:.2f} sec)\n", time() - t_begin) : "")
    t_begin = time()

    # Set objective
    yVal = leaf_nodes[:,7]
    @objective(model, Max, sum(yVal[i] / num_trees * variables_y[i] for i=1:num_leaf_nodes))
    print(verbose > 0 ? format("Set the objective function (duration: {:.2f} sec)\n", time() - t_begin) : "")
    t_begin = time()

    # Set constraints: one node per tree
    @constraint(model, constr_one_node[t = 1:num_trees], sum(variables_y[ptr_leaf_nodes[t]:ptr_leaf_nodes[t+1]-1]) <= 1)
    print(verbose > 0 ? format("Added one_node_per_tree constraints, sum(y[]) <= 1 (duration: {:.2f} sec)\n", time() - t_begin) : "")
    t_begin = time()

    # Set constraints: between two split variables
    @constraint(model, constr_splitVar[i=1:length(splitVars), j=1:(ptr_splits[i+1]-ptr_splits[i]-1)], x[i,j] <= x[i,j+1])
    print(verbose > 0 ? format("Added split variable constraints, x[i,j] <= x[i,j+1] (duration: {:.2f} sec)\n", time() - t_begin) : "")
    t_begin = time()

    # Define all_varindex (will be used when add split constraints)
    print(verbose > 1 ? "Define all_varindex\n" : "")
    all_varindex = zeros(Int, num_leaf_nodes + num_split_nodes)
    print(verbose > 1 ? "Define all_varindex 1\n" : "")
    all_varindex[forest[:, 2] .== 0] = Array(1:num_leaf_nodes)
    print(verbose > 1 ? "Define all_varindex 2\n" : "")
    for i in 1:size(forest, 1)
        if forest[i, 2] > 0
            all_varindex[i] = num_leaf_nodes + get_split_index(splitVars, splits, ptr_splits, forest[i, 4], forest[i, 5])
        end
    end
    
    # Count the number of bounded split constarints
    num_split_constr = 0
    num_bounded_split_constr = 0    
    
    # Define a container for the first split constraint
    first_split_constr = 0
    first_split_constr_updated = false
    split_constr_lazylevel = Array{Int}(undef, 0)
    # split_constr = []
    
    # max_depth_split_constr
    max_depth_split_constr = Inf

    # split constraint info (treeId, nodeId, splitVar, splitVal, leaves)
    split_constr_info = []
        
    
    # Set constraint: split nodes (left, right)
    for t in 1:num_trees
        print(verbose > 1 && t % 10 == 1 ? format("Add split constraints for tree {}\n", t) : "")
        
        # list of tuple (splitVar, lb, ub, leaves)
        split_constr_in_tree = []

        # Save tree info for DFS search
        tree = get_tree_from_forest(forest, ptr_forest, t)
        lchild = convert(Array{Int}, tree[:, 2])
        rchild = convert(Array{Int}, tree[:, 3])
        splitVar = convert(Array{Int}, tree[:, 4])
        splitVal = tree[:, SPLITVAL]
        varindex = all_varindex[ptr_forest[t]:ptr_forest[t+1]-1]
        depth = convert(Array{Int}, tree[:, 10])
        maxyVal = tree[:, MAXYVAL]
        parent = convert(Array{Int}, tree[:, PARENT])

        # Initialize two containers used in DFS search
        # splits_above is the sequence of split info tuple from the root to the parent of the current node
        # A tuple in splits_above contains two elements ("L" or "R" + variable number, variable index of the split)
        # leaves_above is the sequnce of a set of leaves 
        # If the current node is a leaf, it contains the set of leaves in the left branch from the root to the parent of the current node
        # If the current node is a non-leaf, it is same as the case of a leaf except the last two elements
        # The second from the last element is the set of leaves in the left branch
        # The last element is the set of leaves in the right branch
        splits_above = []
        leaves_above = []
        
        function DFS_search(node)
            # (Step 1) Move to the left child to complete left branch
            if lchild[node] > 0
                # Add the current split info
                push!(splits_above, (splitVar[node], 'L', splitVal[node]))
                DFS_search(lchild[node])
            end

            # (Step 2) Move to the right child to complete right branch
            if rchild[node] > 0
                # Add the current split info
                push!(splits_above, (splitVar[node], 'R', splitVal[node]))
                DFS_search(rchild[node])
            end
            
            # (Step 3) Add constraint for the current node
            if lchild[node] == 0
                # If the current node is a leaf, update leaves_above
                push!(leaves_above, [varindex[node]])
            else
                # Compute lower and upper splitVal
                (lower, upper) = (-Inf, Inf)
                for i in reverse(1:length(splits_above))
                    current_split = splits_above[i]
                    if current_split[1] == splitVar[node]
                        if current_split[2] == 'R'
                            lower = max(lower, current_split[3])
                        else
                            upper = min(upper, current_split[3])
                        end
                    end
                        
                    # if lower > -Inf && upper < Inf
                    #     break
                    # end
                end

                if depth[node] <= max_depth_split_constr
                    # Update split_constr_in_tree (left)
                    if lower < splitVal[node]
                        lower = formulation in ["bounded", "bounded2"] ? lower : -Inf
                        push!(split_constr_in_tree, (splitVar[node], lower, splitVal[node], vec(leaves_above[end-1])))
                    else
                        JuMP.set_upper_bound.(variables[vec(leaves_above[end-1])], 0)
                    end

                    # Update split_constr_in_tree (right)
                    if splitVal[node] < upper
                        upper = formulation in ["bounded", "bounded2"] ? upper : Inf
                        push!(split_constr_in_tree, (splitVar[node], splitVal[node], upper, vec(leaves_above[end])))
                    else
                        JuMP.set_upper_bound.(variables[vec(leaves_above[end])], 0)
                    end
                end
                
                # Update leaves_above: combine last two leaves into one
                current_leaves = vcat([leaves_above[end-1] leaves_above[end]])
                leaves_above = leaves_above[1:end-2]
                push!(leaves_above, current_leaves)               
            end

            # (Step 4) Remove the last split info since it moves to parent 
            splits_above = splits_above[1:end-1]
        end
        
        # Start DFS_search from the root node
        DFS_search(1)


        # If formulation is _min or _max, then adjust constraints
        function merge_identical_split_constraints(constrs; keep_order=false)
            num_constrs = length(constrs)
            
            # Sort constrs and keep argsort for the case of keep_order=true
            p = sortperm(constrs)
            constrs = constrs[p]

            # (1) Combine identical splits into one split
            ptr_to_remove = []
            for i = 1:num_constrs-1
                if constrs[i][1:3] == constrs[i+1][1:3]
                    leaves = union(constrs[i][4], constrs[i+1][4])
                    constrs[i+1] = (constrs[i+1][1], constrs[i+1][2], constrs[i+1][3], leaves)
                    push!(ptr_to_remove, i)
                end
            end
            constrs = constrs[setdiff(1:num_constrs, ptr_to_remove)]
            
            if keep_order
                p = p[setdiff(1:num_constrs, ptr_to_remove)]
                p = sortperm(p)
                constrs = constrs[p]
            end
            
            num_constrs = length(constrs)
            constrs
        end

        function unbnded_to_unbndedmax(constrs)
            num_constrs = length(constrs)

            # (1) Combine identical splits into one split
            constrs = merge_identical_split_constraints(constrs)
            num_constrs = length(constrs)

            # (2-1) Add all leaves that satisfiy the split - Forward / less than intervals
            for i = 1:num_constrs-1
                if constrs[i][1:2] == constrs[i+1][1:2]
                    leaves = union(constrs[i][4], constrs[i+1][4])
                    constrs[i+1] = (constrs[i+1][1], constrs[i+1][2], constrs[i+1][3], leaves)
                end
            end

            # (2-2) Add all leaves that satisfiy the split - BackwardForward / greater than intervals
            for i = num_constrs:-1:2
                if constrs[i][[1,3]] == constrs[i-1][[1,3]]
                    leaves = union(constrs[i][4], constrs[i-1][4])
                    constrs[i-1] = (constrs[i-1][1], constrs[i-1][2], constrs[i-1][3], leaves)
                end
            end

            # printfmt("num_constrs = {}\n", num_constrs)
            # for i = 1:num_constrs
            #     if constrs[i][1] == 10
            #         println(constrs[i])
            #     end
            # end

            constrs
        end

        function unbndedmax_to_unbndedmin(constrs)
            ## Assume that constrs are sorted
            num_constrs = length(constrs)

            # (1-1) Backward / less than intervals
            for i = num_constrs:-1:2
                if constrs[i][1:2] == constrs[i-1][1:2]
                    leaves = setdiff(constrs[i][4], constrs[i-1][4])
                    constrs[i] = (constrs[i][1], constrs[i][2], constrs[i][3], leaves)
                end
            end

            # (1-2) Forward / greater than intervals
            for i = 1:num_constrs-1
                if constrs[i][[1,3]] == constrs[i+1][[1,3]]
                    leaves = setdiff(constrs[i][4], constrs[i+1][4])
                    constrs[i] = (constrs[i][1], constrs[i][2], constrs[i][3], leaves)
                end
            end

            # printfmt("num_constrs = {}\n", num_constrs)
            # for i = 1:num_constrs
            #     if constrs[i][1] == 10
            #         println(constrs[i])
            #     end
            # end

            constrs
        end

        function unbndedmin_to_bndedmin(constrs)
            ## Assume that constrs are sorted
            num_constrs = length(constrs)

            constrs_bnded = []
            ub_begin = 1
            ub_end = 0
            lb_begin = 0
            lb_end = 0
            while ub_begin < num_constrs
                # Current splitVar
                splitVar = constrs[ub_begin][1]

                # Compute ub_end, lb_begin, lb_end
                for i = ub_begin:num_constrs
                    if constrs[i][1:2] == (splitVar, -Inf)
                        ub_end = i
                    else
                        break
                    end
                end
                lb_begin = ub_end + 1
                for i = lb_begin:num_constrs
                    if constrs[i][[1,3]] == (splitVar, Inf)
                        lb_end = i
                    else
                        break
                    end
                end
                
                ## Get bounded constr
                for i = lb_begin:lb_end
                    for j = ub_begin:ub_end
                        leaves = intersect(constrs[i][4], constrs[j][4])
                        if length(leaves) > 0
                            # If bounded-interval exists, then add it to constrs_bnded
                            push!(constrs_bnded, (splitVar, constrs[i][2], constrs[j][3], leaves))
                            
                            # If half-interval only includes bounded leaves, then remove it
                            for k in [i, j]
                                constrs[k] = (constrs[k][1], constrs[k][2], constrs[k][3], setdiff(constrs[k][4], leaves))
                            end
                        end
                    end
                end

                # Prepare for the next
                ub_begin = lb_end + 1
            end

            # Remove empty constraints
            ptr_to_remove = []
            for i = 1:num_constrs
                if length(constrs[i][4]) == 0
                    push!(ptr_to_remove, i)
                end
            end
            constrs = constrs[setdiff(1:num_constrs, ptr_to_remove)]

            # Combine half-intervals and bounded-intervals
            constrs = vcat(constrs, constrs_bnded)
            num_constrs = length(constrs)

            # printfmt("num_constrs = {}\n", num_constrs)
            # for i = 1:num_constrs
            #     if constrs[i][1] == 10
            #         println(constrs[i])
            #     end
            # end

            constrs
        end 

        function unbndedmax_to_bndedmax(constrs)
            ## Assume that constrs are sorted
            num_constrs = length(constrs)

            ptr_to_remove = []
            constrs_bnded = []
            ub_begin = 1
            ub_end = 0
            lb_begin = 0
            lb_end = 0
            while ub_begin < num_constrs
                # Current splitVar
                splitVar = constrs[ub_begin][1]

                # Compute ub_end, lb_begin, lb_end
                for i = ub_begin:num_constrs
                    if constrs[i][1:2] == (splitVar, -Inf)
                        ub_end = i
                    else
                        break
                    end
                end
                lb_begin = ub_end + 1
                for i = lb_begin:num_constrs
                    if constrs[i][[1,3]] == (splitVar, Inf)
                        lb_end = i
                    else
                        break
                    end
                end
                
                ## Get bounded constr
                for i = lb_begin:lb_end
                    num_leaves_prev = 0
                    for j = ub_begin:ub_end
                        leaves = intersect(constrs[i][4], constrs[j][4])
                        num_leaves = length(leaves)
                        if num_leaves > num_leaves_prev
                            # If bounded-interval exists, then add it to constrs_bnded
                            push!(constrs_bnded, (splitVar, constrs[i][2], constrs[j][3], leaves))
                            
                            # If half-interval only includes bounded leaves, then add it to the list to remove
                            for k in [i, j]
                                if length(constrs[k][4]) <= num_leaves
                                    push!(ptr_to_remove, k)
                                end
                            end

                            num_leaves_prev = length(leaves)
                        end
                    end
                end

                # Prepare for the next
                ub_begin = lb_end + 1
            end

            # Remove half-intervals if some bounded-interval dominates
            constrs = constrs[setdiff(1:num_constrs, ptr_to_remove)]
            num_constrs = length(constrs)

            # Filter bounded-intervals
            num_constrs_bnded = length(constrs_bnded)
            ptr_to_remove = []
            var_begin = 1
            while var_begin < num_constrs_bnded
                # Current splitVar
                splitVar = constrs_bnded[var_begin][1]
                
                # Compute var_end
                var_end = var_begin
                for i = var_begin:num_constrs_bnded
                    if constrs_bnded[i][1] == splitVar
                        var_end = i
                    else
                        break
                    end
                end

                # Find a pair of constraints that have identical leaves and add the larger interval to the list
                for (i, j) in Iterators.product(var_begin:var_end, var_begin:var_end)
                    if i < j
                        constr_i = constrs_bnded[i]
                        constr_j = constrs_bnded[j]
                        if constr_i[4] == constr_j[4]
                            if constr_i[2] <= constr_j[2] && constr_i[3] >= constr_j[3]
                                push!(ptr_to_remove, i)
                            elseif constr_i[2] >= constr_j[2] && constr_i[3] <= constr_j[3]
                                push!(ptr_to_remove, j)
                            else
                                println("WARNING:: one interval should include the other but NOT HERE")
                                printfmt("constr_i   {}\n", constr_i)
                                printfmt("constr_j   {}\n", constr_j)
                            end
                        end
                    end
                end

                # Prepare for the next
                var_begin = var_end + 1
            end
            constrs_bnded = constrs_bnded[setdiff(1:num_constrs_bnded, ptr_to_remove)]
            num_constrs_bnded = length(constrs_bnded)

            # printfmt("num_constrs_bnded = {}\n", num_constrs_bnded)
            # for i = 1:num_constrs_bnded
            #     if constrs_bnded[i][1] == 10
            #         println(constrs_bnded[i])
            #     end
            # end

            # Combine half-intervals and bounded-intervals
            constrs = vcat(constrs, constrs_bnded)
            num_constrs = length(constrs)

            constrs
        end

        if formulation == "unbounded"
            nothing
        elseif formulation == "bounded"
            nothing
        elseif formulation == "bounded2"
            split_constr_in_tree = merge_identical_split_constraints(split_constr_in_tree, keep_order=false)
        elseif formulation == "unbounded_max"
            ## unbnded -> unbounded_max
            split_constr_in_tree = unbnded_to_unbndedmax(split_constr_in_tree)
        elseif formulation == "unbounded_min"
            ## unbnded -> unbounded_max -> unbounded_min
            split_constr_in_tree = unbnded_to_unbndedmax(split_constr_in_tree)
            split_constr_in_tree = unbndedmax_to_unbndedmin(split_constr_in_tree)
        elseif formulation == "bounded_max"
            ## unbnded -> unbounded_max -> bnded_max
            split_constr_in_tree = unbnded_to_unbndedmax(split_constr_in_tree)
            split_constr_in_tree = unbndedmax_to_bndedmax(split_constr_in_tree)
        elseif formulation == "bounded_min"
            ## unbnded -> unbounded_max -> unbounded_min -> bounded_min
            split_constr_in_tree = unbnded_to_unbndedmax(split_constr_in_tree)
            split_constr_in_tree = unbndedmax_to_unbndedmin(split_constr_in_tree)
            split_constr_in_tree = unbndedmin_to_bndedmin(split_constr_in_tree)
        end
        
        ## Debugging 
        if verbose >= 2 && t == 1
            for (i, j) in Iterators.product(1:length(split_constr_in_tree), 1:length(split_constr_in_tree))
                if i < j                        
                    if formulation == "unbounded_min"
                        if split_constr_in_tree[i][[1,2,4]] == split_constr_in_tree[j][[1,2,4]]
                            println("WARNING:: in split constraints")
                        elseif split_constr_in_tree[i][[1,3,4]] == split_constr_in_tree[j][[1,3,4]]
                            println("WARNING:: in split constraints")
                        end
                    else
                        if split_constr_in_tree[i][[1,4]] == split_constr_in_tree[j][[1,4]]
                            println("WARNING:: in split constraints")
                        end
                    end
                end
            end
        end

        ## Debugging 
        if verbose >= 2 && t == 1
            printfmt("Show the sample of split constraints in {} formulation for indepvar = 1\n", formulation)
            for i = 1:length(split_constr_in_tree)
                if split_constr_in_tree[i][1] == 1
                    println(split_constr_in_tree[i])
                end
            end
        end

        ## Add the split constraints to model
        num_split_constr += length(split_constr_in_tree)
        for (splitVar, lower, upper, leaves) in split_constr_in_tree
            # if maximum(yVal[leaves]) < M || minimum(yVal[leaves]) >= M
            if lower <= -Inf
                # no lowerbound info: leaves <= x
                upper_index = num_leaf_nodes + get_split_index(splitVars, splits, ptr_splits, splitVar, upper)
                constr = @constraint(model, sum(variables[i] for i in leaves) <= variables[upper_index])
            elseif upper >= Inf
                # no upperbound info: leaves <= 1 - x
                lower_index = num_leaf_nodes + get_split_index(splitVars, splits, ptr_splits, splitVar, lower)
                constr = @constraint(model, sum(variables[i] for i in leaves) <= 1 - variables[lower_index])
            else
                # bounded: leaves <= x_j - x_i
                upper_index = num_leaf_nodes + get_split_index(splitVars, splits, ptr_splits, splitVar, upper)
                lower_index = num_leaf_nodes + get_split_index(splitVars, splits, ptr_splits, splitVar, lower)
                constr = @constraint(model, sum(variables[i] for i in leaves) <= variables[upper_index] - variables[lower_index])
                num_bounded_split_constr += 1
            end

            
            MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), constr, lazylevel)

            if ~first_split_constr_updated
                first_split_constr = constr
                first_split_constr_updated = true
            end
        end
    end
    print(verbose > 0 ? format("Added split constraints, left <= x, right <= 1-x (duration: {:.2f} sec)\n", time() - t_begin) : "")
    
    if relax
        t_begin = time()
        JuMP.unset_binary.(x)
        JuMP.set_lower_bound.(x, 0)
        JuMP.set_upper_bound.(x, 1)
        if lazylevel > 0
            # Add dummy binary variables and a constraint to use lazy constraint feature
            @variable(model, dummy[i=1:2], Bin)
            @constraint(model, constr_dummy, sum(dummy[1:2]) <= 1)
        end
        print(verbose > 0 ? format("Relax the problem to LP (duration: {:.2f} sec)\n", time() - t_begin) : "")
    end

    
    if (initsol || fixsol) && ~isequal(solution_X, nothing)
        t_begin = time()
        grb = JuMP.backend(model)
        variable_index = num_leaf_nodes
        for i=1:length(splitVars)
            var = splitVars[i]
            if var in keys(solution_X)
                (val_lb, val_ub) = solution_X[var]
                for j=1:(ptr_splits[i+1]-ptr_splits[i])
#                     variable_index = grb[JuMP.optimizer_index(x[i,j])]
                    variable_index += 1
                    if val_ub <= splits[ptr_splits[i] + j - 1, 2]
                        if initsol
                            Gurobi.set_dblattrelement!(grb.inner, "Start", variable_index, 1)
                        else # fixsol
                            JuMP.unset_binary(x[i,j])
                            JuMP.fix(x[i,j], 1)
                            # @constraint(model, x[i,j] >= 1)
#                                 printfmt("x[{},{}] = 1\n", i, j)
                        end
                    elseif splits[ptr_splits[i] + j - 1, 2] <= val_lb
                        if initsol
                            Gurobi.set_dblattrelement!(grb.inner, "Start", variable_index, 0)
                        else # fixsol
                            JuMP.unset_binary(x[i,j])
                            JuMP.fix(x[i,j], 0)
                            # @constraint(model, x[i,j] <= 0)
#                                 printfmt("x[{},{}] = 0\n", i, j)
                        end
                    else
                        nothing
                    end
                end
            else
                variable_index += ptr_splits[i+1]-ptr_splits[i]
            end
        end
        print(verbose > 0 ? format("Set initial solution or set-up a box (duration: {:.2f} sec)\n", time() - t_begin) : "")
    end
    
    output[:num_splitVars] = length(splitVars)
    output[:num_splitConstrs] = num_split_constr
    output[:bounded_ratio] = num_bounded_split_constr / num_split_constr
    output[:relax] = relax
    output[:runtime_formulate] = time() - t_begin_formulate
    print(verbose > 0 ? format("Total runtime from the beginning to before optimize! (duration: {:.2f} sec)\n", output[:runtime_formulate]) : "")
    
    ## optimize!(model)
    print(verbose > 0 ? "Begin optimize!()\n\n" : "")
    optimize!(model)
    output[:objval] = JuMP.objective_value(model)
    
    (solution_X, solution_Y) = (nothing, nothing)
    if has_values(model)
        solution_X = get_solution_X(JuMP.value.(x), splitVars, splits, ptr_splits)
        solution_Y = get_solution_Y(JuMP.value.(y), num_trees, leaf_nodes, ptr_leaf_nodes)
    end

    grb = JuMP.backend(model)
    output[:num_constrs] = Gurobi.get_intattr(grb.inner, "NumConstrs")
    output[:num_vars] = Gurobi.get_intattr(grb.inner, "NumVars")
    output[:runtime] = Gurobi.get_dblattr(grb.inner, "Runtime")

    if Gurobi.get_intattr(grb.inner, "isMIP") == 1
        output[:MIPGap] = Gurobi.get_dblattr(grb.inner, "MIPGap")
        output[:MIPNodes] = Int(Gurobi.get_dblattr(grb.inner, "NodeCount"))
    else
        output[:MIPGap] = missing
        output[:MIPNodes] = missing
    end
    
    print(verbose > 0 ? format("# of bounded / # of unbounded split constraints = {:.4f}\n", output[:bounded_ratio]) : "")
    print(verbose > 0 ? format("runtime = {:.2f}\n", output[:runtime]) : "")
    print(verbose > 0 ? format("objval = {:.4f}\n", output[:objval]) : "")
    print(verbose > 0 ? "\n" : "")

    if return_solution
        return output, solution_X, solution_Y
    else
        return output
    end
end








directory = "sample_dataset/"
# directory = "../datasets/T5000/redwine_T5000/"
num_trees = 10
verbose = 0

output = solve_MIP_xy(directory, num_trees, formulation="unbounded", lazylevel=2, verbose=verbose)

for (index, value) in output
    println("$index: $value")
end