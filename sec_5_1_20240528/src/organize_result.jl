using CSV, DataFrames, Glob



function main()
    infile = "result_2020-10-30/result.csv"
    outfile = "result.csv"
    df = DataFrame(CSV.File(infile))

    ##### Fill missing to 0
    # println(names(df))
    # rename!(df, "seed"=>:seed)
    df[ismissing.(df.seed), :seed] .= 0

    # filter!(row -> (row.dataset != "sample"), df)
    filter!(row -> (row.dataset != "sample") && (row.method in ["misic", "bounded", "activeleaf", "liftall", "maxclq"]) && row.num_trees >= 100, df)

    selected = [(df.dataset[1], df.num_trees[1], df.seed[1]) for df in groupby(df, [:dataset, :num_trees, :seed]) if maximum(df.time_solve) >= 0]
    filter!(row -> (row.dataset, row.num_trees, row.seed) in selected, df)

    select!(df, [:cnt, :dataset, :num_trees, :seed, :method, :num_vars, :num_constrs, :check_feasibility, :objval, :MIPGap, :time_formulate, :time_solve, :time_cb, :MIPNodes, :num_lazy, :num_call])

    df.GBT = occursin.("_GBT", df.dataset)
    sort!(df, [:GBT, :dataset, :num_trees, :seed, :method])
    df.cnt = 1:nrow(df)
    
    CSV.write(outfile, df)

    ##### method comparison: time, bnbnodes, num_lazy, num_call
    rename!(df, :time_solve => :time)
    df.time2 = df.time - df.time_cb

    df_all = DataFrame(df)

    select!(df, [:dataset, :num_trees, :seed])
    unique!(df)

    for col in [:time, :time2, :MIPNodes, :num_lazy, :num_call], method in unique(df_all.method)
        _df = filter(row -> row.method == method, df_all)
        select!(_df, [:dataset, :num_trees, :seed, col])
        rename!(_df, col => "$(col)_$(method)")
        df = leftjoin(df, _df, on=[:dataset, :num_trees, :seed])
    end


    CSV.write("result_compare.csv", df)
end

main()