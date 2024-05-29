using CSV, DataFrames, Glob, Formatting

include("read_data.jl")
include("utils.jl")

function compute_nonlinear_term(box, df_sol)
    value = 1
    for (k,v) in box
        lb = v[1] > 0 ? df_sol.sol[v[1]] : 0.0
        ub = v[2] < typemax(Int) ? df_sol.sol[v[2]] : 1.0

        value *= ub - lb
    end

    value
end

function compute_nonlinear_objective(result_folder)

    """
    (1) compute all rectangles of nodes where bounds are represented by varindex
        - DFS(for i = 1:nrow), intersect(rec_parent, new_split)
    (2) for each leaf, compute the nonlinear term
    (3) product sum of (pred) and (nonlinear term)
    """

    ### (1) Read data from the folder
    df_forest = DataFrame(CSV.File(result_folder * "df_forest.csv"))
    df_sol = DataFrame(CSV.File(result_folder * "df_sol.csv"))
    intercept = readdlm(result_folder * "intercept.txt")[1,1]

    ### (2) Compute boxes for all nodes (we need boxes of leaves)
    boxes = get_boxes_df_forest(df_forest, "splitvar1", "varindex", 0, typemax(Int))
    
    ### (3) Compute the nonlinear objective value = intercept + sum(pred * nonlinear_term for leaf in all_leaves)


    ### compute nonlinear term and add those as a column of df_sol
    df_forest.ysol = zeros(Float64, nrow(df_forest))
    df_forest.nlsol = zeros(Float64, nrow(df_forest))
    df_forest.ynlgap = zeros(Float64, nrow(df_forest))
    df_sol.nlsol = zeros(Float64, nrow(df_sol))
    df_sol.box = [Dict() for i = 1:nrow(df_sol)]
    for i in 1:nrow(df_forest)
        ysol = sum(df_sol.sol[df_forest.leafbegin[i]:df_forest.leafend[i]])
        nlsol = compute_nonlinear_term(boxes[i], df_sol)
        df_forest[i, "ysol"] = ysol
        df_forest[i, "nlsol"] = nlsol
        df_forest.ynlgap[i] = abs(ysol - nlsol)

        if df_forest[i, "lchild"] == 0
            df_sol.nlsol[df_forest[i, "varindex"]] = nlsol
            df_sol.box[df_forest[i, "varindex"]] = boxes[i]
        end
    end

    df_forest.box = boxes

    CSV.write(result_folder * "df_forest.csv", df_forest)
    CSV.write(result_folder * "df_sol.csv", df_sol)

    df_leaves = filter(row -> row.vartype in ["x or y", "x", "y"], df_sol)
    df_leaves = filter(row -> row.treeid > 0 && row.splitvar1 == 0, df_leaves)

    obj = intercept + sum(row.pred * row.nlsol for row in eachrow(df_leaves))

    obj
end

function get_objective(result_folder)
    DataFrame(CSV.File(result_folder * "result.csv"))[end, "objval"]
end


function get_x_interval_solution(df_sol)
    """
    df = | splitvar1 | interval | sol |
    interval 1's sol    = x_{v,1}
    interval j's sol    = x_{v,j} - x_{v,j-1}
    interval end's sol  = 1 - x_{v,end}
    """
    
    ### Preprocess df_sol (select rows and columns)
    colnames = ["splitvar1", "splitval1", "sol"]
    @assert issubset(colnames, names(df_sol)) "there is a missing column in df_sol"
    select!(df_sol, colnames)
    filter!(row -> row.splitvar1 > 0 && row.splitval1 > 0, df_sol)
    sort!(df_sol, colnames)
    # @info df_sol[1:5,:]

    ### Create a column to save interval solution
    df_sol.intervalsol = zeros(Float64, nrow(df_sol))

    for i in 1:nrow(df_sol)
        if df_sol.splitval1[i] == 1
            ### First, x_{v,1}
            df_sol.intervalsol[i] = df_sol.sol[i]
        else
            ### Middle, x_{v,i} - x_{v,i-1}
            @assert df_sol.splitvar1[i] == df_sol.splitvar1[i - 1]
            @assert df_sol.splitval1[i] == df_sol.splitval1[i - 1] + 1
            df_sol.intervalsol[i] = df_sol.sol[i] - df_sol.sol[i - 1]
        end
        
        ### Last, 1 - x{v,end}
        ### Since the number of intervals is one more than the number of splitval1,
        ### we append the row corresponding to the last interval
        if i == nrow(df_sol) || df_sol.splitvar1[i] != df_sol.splitvar1[i + 1]
            row = Dict("splitvar1" => df_sol.splitvar1[i], "splitval1" => df_sol.splitval1[i] + 1, 
                        "sol" => -1, "intervalsol" => 1 - df_sol.sol[i])
            append!(df_sol, row)
        end
    end

    select!(df_sol, Not("sol"))
    rename!(df_sol, Dict("splitval1" => "interval", "intervalsol" => "sol"))

    # append!(df_sol, df_sol_lastinterval)
    sort!(df_sol, ["splitvar1", "interval", "sol"])

    df_sol
end

function get_x_interval_solution_OPT(df_sol, df_forest)
    df_splits = df_forest_to_df_splits(df_forest)

    optbox = get_optbox(df_sol, df_forest, df_splits; splitval_type="splitval1")
    # println(optbox)

    df = select(filter(row -> row.treeid == 0, df_sol), ["splitvar1", "splitval1", "sol"])
    rename!(df, Dict("splitval1"=>"interval"))
    df.sol .= 0

    for i in 1:nrow(df)
        var = df[i, "splitvar1"]
        val = df[i, "interval"]

        if var in keys(optbox)
            if val - 1 >= optbox[var][1] && val <= optbox[var][2]
                df[i,"sol"] = 1
            end
        else
            df[i,"sol"] = 1
        end

        ### Last, 1 - x{v,end}
        ### Since the number of intervals is one more than the number of splitval1,
        ### we append the row corresponding to the last interval
        if i == nrow(df) || df.splitvar1[i] != df.splitvar1[i + 1]
            row = Dict("splitvar1" => var, "interval" => val + 1, "sol" => 0)
            if var in keys(optbox)
                if val < optbox[var][2]
                    row["sol"] = 1
                end
            else
                row["sol"] = 1
            end
            append!(df, row)
        end
    end

    sort!(df, ["splitvar1", "interval"])

    df
end

function get_x_interval_solution_with_LP_OPT(LP_result_folder, OPT_result_folder)

    df_sol_LP = DataFrame(CSV.File(LP_result_folder * "df_sol.csv"))
    df_LP = get_x_interval_solution(df_sol_LP)

    df_sol = DataFrame(CSV.File(OPT_result_folder * "df_sol.csv"))
    df_forest = DataFrame(CSV.File(OPT_result_folder * "df_forest.csv"))
    df_OPT = get_x_interval_solution_OPT(df_sol, df_forest)
    
    rename!(df_LP, Dict("sol" => "sol_LP"))
    rename!(df_OPT, Dict("sol" => "sol_OPT"))
    df = leftjoin(df_LP, df_OPT, on=["splitvar1", "interval"])    

    df
end

function update_LBnl_result_csv(result_folder)
    if !isfile(result_folder * "result.csv")
        return
    end

    df_result = DataFrame(CSV.File(result_folder * "result.csv"))
    if "LBnl" in names(df_result)
        return
    end
    
    LB_nonlinear = compute_nonlinear_objective(result_folder)
    df_result.LBnl = LB_nonlinear .* ones(nrow(df_result))
    CSV.write(result_folder * "result.csv", df_result)
end

function update_LBnl(dir)
    result_folders = Glob.glob(dir * "**/")
    for result_folder in result_folders
        println(result_folder)
        update_LBnl_result_csv(result_folder)
    end
end

# compute_nonlinear_objective("result_1211/redwine_T1000_Depth10_20_0_liftall-LP_3600_1/")
# result_folders = Glob.glob("result_1211/" * "**/")
# println(format("num folders = {}", length(result_folders)))
# for (i, result_folder) in enumerate(result_folders)
#     # if i <= 318
#     #     continue
#     # end

#     println(i, ", ", result_folder)
#     if isfile(result_folder * "result.csv")
#         update_LBnl_result_csv(result_folder)
#     end
# end
