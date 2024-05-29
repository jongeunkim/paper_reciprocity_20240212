using DataFrames, DataStructures

function remove_dominated(df_splitconstr)
    df_splitconstr.dominated = [false for i = 1:nrow(df_splitconstr)]
    for df in groupby(df_splitconstr, "leaves")
        lb_max = maximum(df.lb)
        ub_min = minimum(df.ub)

        if lb_max == 0 && ub_min == typemax(Int)
            df.dominated .= true
        else
            for row in eachrow(df)
                if row.lb < lb_max || row.ub > ub_min
                    row.dominated = true
                end
            end
        end
    end
    filter!(row -> !row.dominated, df_splitconstr)
    select!(df_splitconstr, Not("dominated"))
    df_splitconstr
end


function convert_to_half_interval(df)
    ### lb
    arr_lbs = sort(unique(df.lb))
    num_lbs = length(arr_lbs)
    df_lb = DataFrame(lb=arr_lbs, ub=maximum(df.ub) * ones(Int, num_lbs), leaves=[SortedSet() for i = 1:num_lbs])
    for row in eachrow(df), row_lb in eachrow(df_lb)
        if row.lb >= row_lb.lb
            union!(row_lb.leaves, row.leaves)
        end
    end
    df_lb = remove_dominated(df_lb)

    ### ub
    arr_ubs = sort(unique(df.ub), rev=true)
    num_ubs = length(arr_ubs)
    df_ub = DataFrame(lb=minimum(df.lb) * ones(Int, num_ubs), ub=arr_ubs, leaves=[SortedSet() for i = 1:num_ubs])
    for row in eachrow(df), row_ub in eachrow(df_ub)
        if row.ub <= row_ub.ub
            union!(row_ub.leaves, row.leaves)
        end
    end
    df_ub = remove_dominated(df_ub)

    df_lb, df_ub
end


function convert_to_bounded_interval(df_lb, df_ub)
    df_lb = select(df_lb, ["lb", "leaves"])
    rename!(df_lb, "leaves" => "leaves_lb")
    df_lb.dummy = ones(Int, nrow(df_lb))
    df_ub = select(df_ub, ["ub", "leaves"])
    rename!(df_ub, "leaves" => "leaves_ub")
    df_ub.dummy = ones(Int, nrow(df_ub))
    
    df_bounded = outerjoin(df_lb, df_ub, on="dummy")
    df_bounded.leaves = [intersect(row.leaves_lb, row.leaves_ub) for row in eachrow(df_bounded)]
    select!(df_bounded, ["lb", "ub", "leaves"])
    filter!(row -> length(row.leaves) > 0, df_bounded)

    ### remove dominated ones
    df_bounded = remove_dominated(df_bounded)

    df_bounded
end


function convert_to_liftall(df)
    df_lb, df_ub = convert_to_half_interval(df)
    df_bounded = convert_to_bounded_interval(df_lb, df_ub)

    df_all = DataFrame(lb=Int[], ub=Int[], leaves=SortedSet[])
    append!(df_all, df_lb)
    append!(df_all, df_ub)
    append!(df_all, df_bounded)
    df_all = remove_dominated(df_all)

    df_all.index = 1:nrow(df_all)

    df_all
end

