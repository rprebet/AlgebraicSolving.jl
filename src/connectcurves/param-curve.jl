
function repl_ind(L,i,x)
    # replace by element x at index i
    return [ L[1:i-1]; x; L[i+1:end] ]
end

function rem_var(f, i, A)
    CE = [collect(coefficients(f)), collect(exponent_vectors(f))]
    CE[2] = [ rem_ind(e, i) for e in CE[2]]
    return array_to_poly(CE, A)
end

function deg_Alg(F, dim)
    # Compute the degree of the the alg set of dimension dim
    # defined by F
    if dim <= 0
        r = rational_parametrization(AlgebraicSolving.Ideal(F))
    else
        vars = gens(parent(first(F)))
        planes = [ sum([ rand(-100:100)*v for v in vcat(vars,1) ])  for _ in 1:dim]
        r = rational_parametrization(AlgebraicSolving.Ideal(vcat(F,planes)))
    end
    return degree(r.elim)
end


function compute_param(F)
    varias = gens(parent(first(F)))
    N = length(varias)

    if N == 2
        return F
    end

    DEG = deg_Alg(F,1)
    if deg_Alg(vcat(F,varias[N-1]-rand(-100:100)),0) != DEG
        error("The curve is not in generic position")
    end

    PARAM = []
    rr = 0
    _values, i = [], 1
    A, (s,t) = polynomial_ring(QQ,[:s,:t])
    while length(_values) < DEG+2
        try
            Feval = [ evaluate(f, repl_ind(varias, N-1, QQ(i))) for f in F ]
            Feval = [ rem_var(f, N-1, A) for f in Feval ]
            r = rational_parametrization(AlgebraicSolving.Ideal(Feval))
            #ivev = findfirst(t->t==R.S[end-1],r.vars)
            lc = leading_coefficient(r.elim)
            rr = [ p/lc for p in vcat(r.elim,r.param) ]
            push!(PARAM, rr)
            push!(_values, QQ(i))
            i+=1
        catch
            error("bad specialization choice")
        end
    end

    CURVEPARAM = []

    T, (x1,x2) = polynomial_ring(QQ, [:x1,:x2])

    for count in 1:N-1
        COEFFS = []
        for deg in 0:DEG
            _evals = [ coeff(PARAM[i][count],deg) for i in 1:length(PARAM)]
            S, u = polynomial_ring(QQ,"u")
            push!(COEFFS, interpolate(S, map(QQ,_values[1:length(PARAM)]), _evals))
        end

        C = [ collect(coefficients(c)) for c in COEFFS ]
        POL = sum([ C[i][j]*x1^(j-1)*x2^(i-1) for i in 1:length(C) for j in 1:length(C[i])])
        
        push!(CURVEPARAM, POL)
    end

    return CURVEPARAM
end
