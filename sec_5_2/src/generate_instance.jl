
function generate_instance1(M, N, C, K, seed)
    @assert C in [2,3] "wrong dimension"

    """
    Input:
        M: # of supply nodes
        N: # of demand nodes
        C: # of commodities
        K: grid density (K x K or K x K x K)
    Output:
        s: supply
        d: demand
        c: capacity
        cost for each arc
            for each arc, generate cost data -> train tree ensemble
    """

    Random.seed!(seed)

    s = rand(K:(2*K-1), M, C)
    d = rand(K:(2*K-1), N, C)
    c = reshape([Int(ceil((s[i,k] + d[j,k])/4.0)) for i=1:M, j=1:N, k=1:C], (M, N, C))
    B = rand([3*K, 4*K, 5*K], M, N)
    b = C == 2 ? [3, 8] : [3, 5, 8] 

    s, d, c, B, b
end

function generate_instance(M, N, C, seed; verbose=0)
    Random.seed!(seed)
    
    # Total supply and demand
    # total_demand = 100 .* ((1:C).^2)
    total_demand = 100 .* [3^(i-1) for i in 1:C]
    total_supply = 1.2 .* total_demand

    # Generate supply and demand distributions
    ps = 1 .+ rand(M,C)
    pd = 1 .+ rand(N,C)

    # Compute supply, demand, capacity
    s = Int.(ceil.([total_supply[k] * ps[i,k] / sum(ps[:,k]) for i in 1:M, k in 1:C]))
    d = Int.(floor.([total_demand[k] * pd[i,k] / sum(pd[:,k]) for i in 1:N, k in 1:C]))
    c = Int.(ceil.([0.25 * (s[i,k] + d[j,k]) for i=1:M, j=1:N, k=1:C]))

    # Check
    @assert size(s) == (M,C)
    @assert size(d) == (N,C)
    @assert size(c) == (M,N,C)

    if verbose >= 1
        for k in 1:C
            println("s[:,$k] = $(s[:,k])")
            println("d[:,$k] = $(d[:,k])")
        end
        for i=1:M, j=1:N
            println("c[$i,$j,:] = $(c[i,j,:])")
        end
    end

    s, d, c
end

