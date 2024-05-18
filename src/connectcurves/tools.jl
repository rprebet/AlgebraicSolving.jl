function order_permut2d(L)
    # Create a list of tuples with elements and their corresponding indices
    LL = [(L[i][j], (i, j)) for i in eachindex(L) for j in eachindex(L[i])]
    # Sort the enumerated list based on the values
    sorted_LL = sort(LL, by = x -> x[1])
    # Extract the sorted values and their corresponding indices
    #sorted_L = [pair[1] for pair in sorted_LL]
    sorted_ind = [pair[2] for pair in sorted_LL]
    
    return sorted_ind
end

function diff(p, v, n)
    dp = p
    for j in 1:n
        dp = derivative(dp, v)
    end
    return dp
end

function diff_list(p, v, n)
    Ldp = [p]
    for j in 1:n
        push!(Ldb,derivative(Ldp[end], v))
    end
    return Ldp
end