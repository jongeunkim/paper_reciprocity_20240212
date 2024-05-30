using JuMP, Statistics
import Gurobi, Random, CSV, Glob, MathOptInterface, DecisionTree, DataFrames, LinearAlgebra
const DT = DecisionTree
const DF = DataFrames
const MOI = MathOptInterface
const GRB_ENV = Gurobi.Env()

include("generate_instance.jl")
include("generate_forest.jl")
include("build_mip.jl")
include("generate_bndconstr.jl")

function output_model(model, isMIP)
    objval = -1
    bnbcnt = -1
    mipgap = -1
    objbnd = -1

    if isMIP
        bnbcnt = JuMP.node_count(model)
        objbnd = JuMP.objective_bound(model)
    end

    if JuMP.has_values(model)
        objval = JuMP.objective_value(model)
        if isMIP
            mipgap = JuMP.relative_gap(model)
        end
    else
        """
        infeasible or not found a solution
        """
    end

    return objval, objbnd, mipgap, bnbcnt 
end

function run(formulation_name_in_paper, seed, M, N, C, max_depth, n_trees; timelimit=3600, verbose=0)
    formulation_name_in_paper_to_formulation = Dict("TEOM" => "TreeM", "TEOR" => "TreeO", "TEOC" => "TreeCvx", "DCC" => "DCC", "DLog" => "DLog", "MC" => "MC")
    @assert formulation_name_in_paper in keys(formulation_name_in_paper_to_formulation)
    formulation  = formulation_name_in_paper_to_formulation[formulation_name_in_paper]

    s, d, c = generate_instance(M, N, C, seed; verbose=verbose)

    cost_forest = Dict()
    cost_forestmodel = Dict()
    for i=1:M, j=1:N
        cost_forest[i,j], cost_forestmodel[i,j] = generate_cost_forest(c[i,j,:], max_depth, n_trees, seed + 10^2 * i + 10^4 * j)
    end

    # Input for an optimization problem
    # M, N, C
    # s, d, c
    # cost_forest

    model, z, w, wt = build_mip(formulation, M, N, C, s, d, c, cost_forest)
    
    JuMP.set_optimizer_attribute(model, "TimeLimit", timelimit)
    JuMP.set_optimizer_attribute(model, "LogToConsole", 0)
    t_begin = time()
    optimize!(model)
    runtime = time() - t_begin

    objval, objbnd, mipgap, bnbcnt = output_model(model, true)

    output = Dict()
    output["runtime"] = runtime
    output["objval"] = objval
    output["objbnd"] = objbnd
    output["mipgap"] = mipgap 
    output["bnbcnt"] = bnbcnt

    ## Obtain LP value
    var = JuMP.all_variables(model)
    for v in var
        if JuMP.is_binary(v)
            JuMP.unset_binary(v)
            JuMP.set_lower_bound(v, 0.0)
            JuMP.set_upper_bound(v, 1.0)
        end
    end
    optimize!(model)
    output["objval_LP"] = JuMP.objective_value(model)

    return output
end