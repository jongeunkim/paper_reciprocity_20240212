using CSV, DataFrames, Dates, DelimitedFiles, Formatting, Glob, Logging, Missings, Printf, Random

include("ensembleopt2.jl")
include("save_result.jl")

function main(args)
    ##### input
    dir = args[1]
    datasetdir = args[2]
    num_trees = parse(Int, args[3])
    seed = parse(Int, args[4])
    method = args[5]
    timelimit = parse(Float64, args[6])
    preprocessing_level = parse(Int, args[7])

    try
        mkdir(dir)
    catch
        nothing
    end


    treeids = Array(1:num_trees)
    if seed > 0
        T = length(Glob.glob(datasetdir *  "*tree*.txt"))
        Random.seed!(seed)
        treeids = Random.randperm(T)[1:num_trees]
    end

    ## dummy (compile)
    if num_trees >= 100
        output = ensembleopt(dir, datasetdir, Array(1:10), method, timelimit, loglevel="Info")
    end
    
    ## main
    output = ensembleopt(dir, datasetdir, treeids, method, timelimit, loglevel="Debug", preprocessing_level=preprocessing_level)
    merge!(output, Dict(:dataset=>datasetdir, :num_trees=>num_trees, :seed=>seed, :method=>method, :timelimit=>timelimit, :preprocessing_level=>preprocessing_level))
    save_result(output, dir * "result.csv")
end

main(ARGS)
