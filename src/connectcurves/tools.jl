function order_permut2d(L)
    # Create a list of tuples with elements and their corresponding indices
    LL = [(L[i][j], (i, j)) for i in eachindex(L) for j in eachindex(L[i])]
    # Sort the enumerated list based on the values
    sorted_LL = sort(LL, by = x -> x[1])
    # Extract the sorted values and their corresponding indices
    sorted_L = [pair[1] for pair in sorted_LL]
    sorted_ind = [pair[2] for pair in sorted_LL]
    
    return sorted_L, sorted_ind
end

function Arb_to_rat(x)
	r = radius(x)
	return map(simplest_rational_inside, [x-2*r, x+2*r])
end

function rat_to_Arb(x, prec)
    x1,x2 = x
    xm, xd = ArbField(prec)((x1+x2)/2), ArbField(prec)(x2-x1)
    return ball(xm,xd)
end

function evaluate_Arb(f, x, prec)
	cf = coefficients_of_univariate(f)
	return evalpoly(ArbField(prec)(x), cf) 
end