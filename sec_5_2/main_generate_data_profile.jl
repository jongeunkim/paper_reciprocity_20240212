using DataFrames, Plots, Plots.PlotMeasures
import CSV
ENV["GKSwstype"] = "nul"

function generate_data_profile(   X, fig_filename;
                                title="", label=nothing, xlabel=nothing, ylabel=nothing, 
                                legend=:bottomright, left_margin=5mm, bottom_margin=5mm,
                                xaxis=nothing, yaxis=nothing)
    (m,n) = size(X)
    Y = 0:1/m:1
    X_sorted = vcat(zeros(1,n), sort(X, dims=1))

    fig = plot( X_sorted, Y, linetype=:steppre, title=title, xlabel=xlabel, ylabel=ylabel, label=reshape(label, 1, :), 
                legend=legend, left_margin=left_margin, bottom_margin=bottom_margin)
    
    isnothing(xaxis) ? xaxis!([0,maximum(X)]) : xaxis!(xaxis)
    isnothing(yaxis) ? yaxis!([0,1]) : yaxis!(yaxis)
    
    savefig(fig_filename)
    
    fig
end

function main()
    # Parameters for easy instances
    result_filename = "result_data/easy.csv"
    fig_filename = "result_data/easy.png" # png or eps
    formulations = ["TEOC", "TEOR", "TEOM", "MC", "DCC", "DLog"]
    timelimit_sec = 3600.

    # Parameters for hard instances
    # result_filename = "result_data/hard.csv"
    # fig_filename = "result_data/hard.eps" # png or eps
    # formulations = ["TEOC", "TEOR", "TEOM"]
    # timelimit_sec = 10800.

    # Construct solution time matrix X whose rows are instances and columns are formulation 
    df = DataFrame(CSV.File(result_filename))
    filter!(row -> row.formulation in formulations, df)
    df_f = [filter(row -> row.formulation == f, df) for f in formulations]
    @assert length(unique(nrow.(df_f))) == 1 "Not the same number of instances solved, $formulations = $(nrow.(df_f))"
    X = hcat([df1.runtime for df1 in df_f]...)

    title = ""
    xlabel = "Solution time (sec)"
    ylabel = "Fraction of solved instances"
    generate_data_profile(X, fig_filename; 
                        title=title, label=formulations, xlabel=xlabel, ylabel=ylabel, 
                        legend=:bottomright, left_margin=5mm, bottom_margin=5mm,
                        xaxis=[0., timelimit_sec], yaxis=[0., 1.])
end

main()