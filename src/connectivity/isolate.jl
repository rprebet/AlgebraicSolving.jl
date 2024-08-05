include("src/usolve/usolve.jl")

function isolate(f; prec = 32, software="msolve")
    # Univariate isolation using usolve
    @assert is_univariate(f) "Not univariate polynomial"
	if software == "usolve"
		return getindex.(usolve(f, precision = prec, uspath="src/usolve/usolve"),2)
    else
        fvar = vars(f)[1]
        sols = real_solutions(Ideal(vcat([f], [ t for t in gens(parent(f)) if t!=fvar ])), precision=prec, interval=true)
        return getindex.(sols, var_index(fvar))
    end
end

function isolate_eval(f, ivar, val; prec=64)
    # univariate isolation of roots of a bivariate polynomial f
    # whose ivar-th variable is evaluated at val
    # uses msolve
    sols = real_solutions(Ideal([f, gens(parent(f))[ivar] - val]), precision=prec, interval=true)
    return [ s[ivar%2+1] for s in sols ]
end