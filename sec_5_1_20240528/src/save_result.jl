using CSV, DataFrames

function save_result(output, file)
    df = DataFrame()
    cnt = 1

    if isfile(file)
        df = DataFrame(CSV.File(file))
        cnt = df.cnt[end] + 1
    end

    output["cnt"] = cnt

    output = Dict( Symbol(k)=>v for (k,v) in output )
    push!(df, output, cols=:union)
    CSV.write(file, df)
end