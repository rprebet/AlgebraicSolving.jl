include("src/usolve/usolve.jl")

function isolate(f; prec = 32, software="usolve")
    # Univariate isolation using usolve
    @assert is_univariate(f) "Not univariate polynomial"
	if software == "usolve"
		return usolve(f, precision = prec, uspath="src/usolve/usolve")
    end
end

function isolate_eval(f, ivar, val; prec=64)
    # univariate isolation of roots of a bivariate polynomial f
    # whose ivar-th variable is evaluated at val
    # uses msolve
    sols = inter_solutions(Ideal([f, gens(R)[ivar] - val]), precision=prec)
    return [ s[ivar%2+1] for s in sols ]
end