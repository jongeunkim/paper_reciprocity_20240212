function generate_bndconstr(intervals)
    """
    intervals : 2-dim matrix
        row - intervals
        col = [lb, ub]
    """

    num_intervals = size(intervals, 1)
    vals = sort(unique(intervals))
    vals_lb = sort(unique(intervals[:,1]))
    vals_ub = sort(unique(intervals[:,2]))

    
    list_bndconstr1 = [(lb, ub, [i for i in 1:num_intervals if lb <= intervals[i,1] && intervals[i,2] <= ub])
        for lb in vals_lb, ub in vals_ub if lb < ub]
    filter!(r -> length(r[3]) > 0, list_bndconstr1)

    # Remove redundant: if there is val that separate intervals
    list_bndconstr2 = []
    for (lb, ub, indices) in list_bndconstr1
        is_redundant = false
        vals_btw = filter(r -> lb < r && r < ub, vals)
        for v in vals_btw
            indices_left = [i for i in indices if intervals[i,2] <= v]
            indices_right = [i for i in indices if intervals[i,1] >= v]

            if length(indices_left) + length(indices_right) == length(indices)
                # this constraint is redundant
                is_redundant = true
                break
            end
        end

        if !is_redundant
            push!(list_bndconstr2, (lb,ub,indices))
        end
    end

    return list_bndconstr2
end
