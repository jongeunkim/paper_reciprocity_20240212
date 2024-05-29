using DataFrames, CSV, Formatting, DataStructures, Statistics

include("utils.jl")
include("callbacks.jl")

function get_num_interval(df_splits, d)
    nrow(df_splits[df_splits.splitvar1 .== d, :]) + 1
end

function get_grid_x_1D(d, df_sol, df_forest, df_splits)
    ptr_forest = get_pointers(df_forest, "treeid")

    ##### Create a grid
    num_intervals = get_num_interval(df_splits, d)

    ##### 
    grid = zeros(num_intervals)
    sol = df_sol[df_sol.splitvar1 .== d, "sol"]

    for i in 1:num_intervals
        frac1 = 0
        if i == 1
            grid[i] = sol[1]
        elseif i == num_intervals
            grid[i] = 1 - sol[end]
        else
            grid[i] = sol[i] - sol[i-1]
        end
    end

    grid
end

function get_grid_x_2D(d1, d2, df_sol, df_forest, df_splits, celltype)
    grid1 = get_grid_x_1D(d1, df_sol, df_forest, df_splits)
    grid2 = get_grid_x_1D(d2, df_sol, df_forest, df_splits)

    n1 = length(grid1)
    n2 = length(grid2)

    grid = zeros((n1,n2))

    for i=1:n1, j=1:n2
        if celltype == "min"
            grid[i,j] = min(grid1[i], grid2[j])
        elseif celltype == "product"
            grid[i,j] = grid1[i] * grid2[j]
        end
    end

    grid
end


# function get_optimal_box(df_sol, df_forest, df_splits; splitval_type="splitval1")
#     df_optleaves = filter(row -> row.treeid > 0 && row.sol > 0.5, df_sol)
#     ptr_forest = get_pointers(df_forest, "treeid")

#     optimal_box = SortedDict()
#     for row in eachrow(df_optleaves)
#         forest_index = get_forest_index(ptr_forest, row.treeid, row.nodeid)
#         rec = get_rectangle(df_forest, forest_index)
#         optimal_box = intersect_rectangles(optimal_box, rec)
#     end

#     if splitval_type == "splitval1"
#         for (k,v) in optimal_box
#             num_intervals = get_num_interval(df_splits, k)
#             new_v = [-1,-1]
#             for i in 1:2
#                 if v[i] == 0
#                     new_v[i] = 0
#                 elseif v[i] == typemax(Int)
#                     new_v[i] = num_intervals
#                 else
#                     new_v[i] = df_splits[df_splits.varindex .== v[i], :].splitval1[1]
#                 end
#             end
#             optimal_box[k] = new_v
#         end
#     end

#     optimal_box
# end

function optimal_box_to_grid(d1, d2, optimal_box, df_splits)
    num_intervals = [get_num_interval(df_splits, d1), get_num_interval(df_splits, d2)]
    grid = zeros((num_intervals[1], num_intervals[2]))

    dims = [d1, d2]
    lbs = [-1, -1]
    ubs = [-1, -1]

    for k in 1:2
        d = dims[k]
        if d in keys(optimal_box)
            (lbs[k], ubs[k]) = (optimal_box[d][1]+1, optimal_box[d][2])
        else
            (lbs[k], ubs[k]) = (1, num_intervals[k])
        end
    end

    for i in lbs[1]:ubs[1]
        for j in lbs[2]:ubs[2]
            grid[i,j] = 1
        end
    end

    grid
end

function print_grid(grid)
    (M, N) = size(grid)

    for i in 1:M
        for j in 1:N
            print(format("{}\t", grid[i,j]))
        end
        println()
    end

end

# dir_LP = "result_1201/Concrete_T5000_GBT_50_0_liftall-LP_3600_1/"
# dir_IP = "result_1201/Concrete_T5000_GBT_50_0_liftall_3600_1/"
# num_splitvar = 8
# arr_perc = []
# for d1 = 1:num_splitvar, d2=(d1+1):num_splitvar
#     # println(d1, ", ", d2)
#     df_sol = DataFrame(CSV.File(dir_LP * "df_sol.csv"))
#     df_forest = DataFrame(CSV.File(dir_LP * "df_forest.csv"))
#     df_splits = DataFrame(CSV.File(dir_LP * "df_splits.csv"))
    

#     grid_x_LP = get_grid_x_2D(d1, d2, df_sol, df_forest, df_splits, "product")
#     # print_grid(grid_x_LP)

#     df_sol = DataFrame(CSV.File(dir_IP * "df_sol.csv"))
#     df_forest = DataFrame(CSV.File(dir_IP * "df_forest.csv"))
#     df_splits = DataFrame(CSV.File(dir_IP * "df_splits.csv"))
#     optimal_box = get_optimal_box(df_sol, df_forest, df_splits)

#     grid_y = optimal_box_to_grid(d1, d2, optimal_box, df_splits)

#     (m,n) = size(grid_y)
#     perc = 0
#     for i=1:m, j=1:n
#         if grid_y[i,j] > 0.5
#             perc += grid_x_LP[i,j]
#         end
#     end

#     # println(perc)
#     push!(arr_perc, perc)
# end

# println(mean(arr_perc))
# println(arr_perc)

function main_LP_in_optimalbox()
    ##### Maximize y variables in the optimal box
    df_result = DataFrame()
    for dataset in ["concrete_T1000_Depth10", "Redwine_T5000_GBT", "Whitewine_T5000_GBT", "Concrete_T5000_GBT", "redwinepair", "assortment_T5000_depth10"], num_trees in [40,50,60,70,80,90,100], method in ["misic", "liftx", "liftall"]
        dir_LP = "result_1201/$(dataset)_$(num_trees)_0_$(method)-LP_3600_1/"
        dir_IP = "result_1201/$(dataset)_$(num_trees)_0_$(method)_3600_1/"

        if !isdir(dir_LP) || !isdir(dir_IP)
            continue
        end

        df_splits = DataFrame(CSV.File(dir_IP * "df_splits.csv"))
        num_splitvars = maximum(df_splits.splitvar1)

        df_sol = DataFrame(CSV.File(dir_IP * "df_sol.csv"))
        df_forest = DataFrame(CSV.File(dir_IP * "df_forest.csv"))
        df_splits = DataFrame(CSV.File(dir_IP * "df_splits.csv"))
        optimal_box = get_optimal_box(df_sol, df_forest, df_splits)

        df_sol = DataFrame(CSV.File(dir_LP * "df_sol.csv"))
        df_forest = DataFrame(CSV.File(dir_LP * "df_forest.csv"))
        df_splits = DataFrame(CSV.File(dir_LP * "df_splits.csv"))
        grids = []
        for d = 1:num_splitvars
            push!(grids, get_grid_x_1D(d, df_sol, df_forest, df_splits))
        end

        value_1Davg = 0
        value_FD = 1
        for d = 1:num_splitvars
            if d in keys(optimal_box)
                value_1D = sum([grids[d][i] for i = optimal_box[d][1]+1:optimal_box[d][2]])
                value_FD = min(value_1D, value_FD)
                value_1Davg += value_1D / num_splitvars
            end
        end

        append!(df_result, Dict("dataset"=>dataset, "num_trees"=>num_trees, "method"=> method, "value_1Davg"=>value_1Davg, "value_FD"=>value_FD), cols=:union)
    end

    header = ["dataset","num_trees","method","value_1Davg", "value_FD"]
    select!(df_result, header)
    sort!(df_result, header)
    CSV.write("LP_in_optimalbox.csv", df_result)
end
main_LP_in_optimalbox()

# function num_unitbox(rec)
#     num = 1
#     for (k,v) in rec
#         num *= (rec[k][2] - rec[k][1])
#     end
#     num
# end

# for dataset in ["Redwine_T5000_GBT", "Whitewine_T5000_GBT", "Concrete_T5000_GBT", "redwinepair"], method in ["misic", "liftx", "liftall"], num_trees in [50]
#     dir_LP = "result_1201/$(dataset)_$(num_trees)_0_$(method)-LP_3600_1/"
#     dir_IP = "result_1201/$(dataset)_$(num_trees)_0_$(method)_3600_1/"

#     df_splits = DataFrame(CSV.File(dir_IP * "df_splits.csv"))
#     num_splitvars = maximum(df_splits.splitvar1)

#     df_sol = DataFrame(CSV.File(dir_IP * "df_sol.csv"))
#     df_forest = DataFrame(CSV.File(dir_IP * "df_forest.csv"))
#     df_splits = DataFrame(CSV.File(dir_IP * "df_splits.csv"))
#     optimal_box = get_optimal_box(df_sol, df_forest, df_splits)
#     println(num_unitbox(optimal_box))
# end