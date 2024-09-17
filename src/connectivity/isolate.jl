include("src/usolve/usolve.jl")

function isolate(f; prec = 32, software="usolve")
    # Univariate isolation using usolve
    if total_degree(f) == 0
        return []
    end
    @assert is_univariate(f) "Not univariate polynomial"
	if software == "usolve"
		return usolve(f, precision = prec, uspath="AlgebraicSolving.jl/src/connectivity/src/usolve/usolve", output="inter")
    else
        #println(f)
        sols = real_solutions(Ideal([change_ringvar(f)]), precision=prec, interval=true)
        return getindex.(sols, 1)
    end
end

function isolate_eval(f, ivar, val; prec=64, software="usolve")
    #print(val)
    # univariate isolation of roots of a bivariate polynomial f
    # whose ivar-th variable is evaluated at val
    # uses msolve
    fev = change_ringvar(evaluate(f, [ivar], [val]))
    CD = lcm(map(denominator, collect(coefficients(fev))))
    fev *= CD

    return isolate(fev, prec=prec, software=software)

    
end