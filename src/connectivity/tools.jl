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
        push!(Ldp, derivative(Ldp[end], v))
    end
    return Ldp
end

function trimat_rand(A, n; up=true, range=-100:100)
    if up
        return [ i==j ? one(A) : (i<j ? A(rand(range)) : zero(A)) for i in 1:n, j in 1:n ]
    else
        return [ i==j ? one(A) : (i>j ? A(rand(range)) : zero(A)) for i in 1:3, j in 1:3 ]
    end
end