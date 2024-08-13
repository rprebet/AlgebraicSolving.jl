#println("Loading libraries and data..")
#using Pkg
#using Revise
#Pkg.activate("AlgebraicSolving.jl-official")
#using AlgebraicSolving
#import Nemo: coefficients, exponent_vectors, coeff, interpolate, MPolyBuildCtx, push_term!, finish, gens

export compute_param

function rem_var(f, i)
    # Remove the occurence of the ith variable
    R = parent(f)
    A, Avars = polynomial_ring(base_ring(R), rem_ind(R.S, i))
    CE = [collect(coefficients(f)), collect(exponent_vectors(f))]
    CE[2] = [ rem_ind(e, i) for e in CE[2]]

    # Reconstruct the polynomial in the ring A
    C = MPolyBuildCtx(A)
    R = base_ring(A)
    for i in eachindex(CE[1])
        push_term!(C, R(CE[1][i]), CE[2][i]);
    end
    return finish(C)
end

function rem_ind(L,i::Int)
    # remove element at index i
    return [ L[1:i-1]; L[i+1:end] ]
end

function deg_Alg(F, dim)
    if dim <= 0
        r = rational_parametrization(Ideal(F))
    else
        varias = gens(parent(first(F)))
        planes = [ sum([ rand(-100:100)*v for v in vcat(varias,1) ])  for _ in 1:dim]
        r = rational_parametrization(Ideal(vcat(F,planes)))
    end
    return degree(r.elim)
end


function compute_param(F; lfs = false)
    R = parent(first(F))
    varias, N = gens(R), nvars(R)
    DEG = deg_Alg(F,1)

    if !lfs
        # Identification of two generic variables for the parametrization")
        ivarias_gen = Vector{Int}(undef,0)
        ind = N
        while true
            NEW_DEG = deg_Alg(vcat(F, varias[ind]-rand(-100:100)), 0)
            if NEW_DEG == DEG
                push!(ivarias_gen, ind)
            end
            (length(ivarias_gen) < 2 && ind > 0) || break
            ind -= 1
        end
        reverse!(ivarias_gen)
        nb_lf = 2-length(ivarias_gen)

        if ivarias_gen != [N-1,N]
            # Add new generic variables at the end if necessary
            newvarias_S = [R.S[i] for i in vcat(setdiff(1:N, ivarias_gen), ivarias_gen)]
            append!(newvarias_S, [:A,:B][nb_lf:-1:1] )

            # We change the polynomial ring and add linear form(s) if necessary
            Rold, Nold = R, N
            R, varias = polynomial_ring(base_ring(Rold), newvarias_S)
            N = length(varias)

            index_permut = Dict((v, i) for (i, v) in enumerate(newvarias_S))
            permut_varias = [ index_permut[v] for v in Rold.S ]
            # TODO: permut exponents and use polynomial constructor 
            F = [ evaluate(f, [ varias[permut_varias[i]] for i in 1:Nold ]) for f in F ] 

            if nb_lf > 0
                lf_cfs = [ rand(ZZ(-100):ZZ(100), N) for _ in 1:nb_lf]
                append!(F, [ transpose(lf)*varias  for lf in lf_cfs ])
            else
                lf_cfs = Vector{Vector{ZZRingElem}}(undef,0)
            end
        else
            lf_cfs = Vector{Vector{ZZRingElem}}(undef,0)
        end
    else
        # Add new generic variables at the end if necessary
        newvarias_S = vcat(R.S, [:B,:A])

        # We change the polynomial ring and add linear form(s) if necessary
        Rold, Nold = R, N
        R, varias = polynomial_ring(base_ring(Rold), newvarias_S)
        N = length(varias)

        F = [ evaluate(f, [ varias[i] for i in 1:Nold ]) for f in F ] 
        lf_cfs = [ rand(ZZ(-100):ZZ(100), N) for _ in 1:2]
        append!(F, [ transpose(lf)*varias  for lf in lf_cfs ])
    end
    # Compute DEG+1 evaluations of the param (whose deg is bounded by DEG)
    PARAM  = Vector{Vector{AlgebraicSolving.QQPolyRingElem}}(undef,0)
    _values = Vector{QQFieldElem}(undef,0)
    i = 1
    while length(_values) < DEG+2
        if i > 2*(DEG+2)
            error("Too many bad specializations")
        end

        #try
            Feval = [ evaluate(f, [N-1], [QQ(i)]) for f in F ]
            Feval = [ rem_var(f, N-1) for f in Feval ]
            r = rational_parametrization(Ideal(Feval))
            # For lifting: the same variable must be chosen for the param
            if  r.vars == rem_ind(R.S, N-1)
                lc = leading_coefficient(r.elim)
                # TODO: why dividing by lc ?
                rr = [ p/lc for p in vcat(r.elim, r.denom, r.param) ]
                push!(PARAM, rr)
                push!(_values, QQ(i))
            else
                println("bad specialization: ", i)
            end
        #catch
        #    error("bad specialization: ", i)
        #end
        i += 1
    end

    # Interpolate each coefficient of each poly in the param 
    POLY_PARAM = Vector{QQMPolyRingElem}(undef,0)
    T, (x,y) = polynomial_ring(QQ, [:x,:y])
    for count in 1:N
        COEFFS = []
        for deg in 0:DEG
            _evals = [ coeff(PARAM[i][count], deg) for i in 1:length(PARAM)]
            push!(COEFFS, interpolate(polynomial_ring(QQ,"u")[1], _values, _evals))
        end

        C = [ collect(coefficients(c)) for c in COEFFS ]
        POL_term = [C[i][j]*x^(j-1)*y^(i-1) for i in 1:length(C) for j in 1:length(C[i])]
        POL = length(POL_term) > 0 ? sum(POL_term) : T(0) 
        
        push!(POLY_PARAM, POL)
    end

    return [R.S, lf_cfs, POLY_PARAM[1], POLY_PARAM[2], POLY_PARAM[3:end]]
end

#=
# Tests
R, (a,b,c) = polynomial_ring(QQ, [:a,:b,:c])
F = [-62*a^2-24*a*b+83*a*c-46*b^2-45*b*c+16*c^2+13*a-65*b-84*c-19, 
92*a^2-a*b+39*a*c-16*b^2+61*b*c+78*c^2+17*a+49*b-97*c-56
]

R, (x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])
F = [
x1^4 - 12*x1^3*x2 - 40*x1^3*x3 + 12*x1^3*x4 + 56*x1^2*x2^2 + 320*x1^2*x2*x3 - 148*x1^2*x2*x4 + 802*x1^2*x3^2 + 56*x1^2*x3*x4 + 288*x1^2*x4^2 - 20*x1^2 - 120*x1*x2^3 - 880*x1*x2^2*x3 + 576*x1*x2^2*x4 - 4012*x1*x2*x3^2 + 224*x1*x2*x3*x4 - 1968*x1*x2*x4^2 + 120*x1*x2 - 8040*x1*x3^3 - 3508*x1*x3^2*x4 - 3264*x1*x3*x4^2 + 400*x1*x3 + 1512*x1*x4^3 - 120*x1*x4 + 100*x2^4 + 800*x2^3*x3 - 760*x2^3*x4 + 5620*x2^2*x3^2 - 80*x2^2*x3*x4 + 3964*x2^2*x4^2 - 200*x2^2 + 16080*x2*x3^3 - 3436*x2*x3^2*x4 - 1168*x2*x3*x4^2 - 800*x2*x3 - 9576*x2*x4^3 + 760*x2*x4 + 40401*x3^4 + 59496*x3^3*x4 + 72556*x3^2*x4^2 - 4020*x3^2 + 37296*x3*x4^3 - 2960*x3*x4 + 15876*x4^4 - 2484*x4^2 + 64,
 -40*x1^3 + 320*x1^2*x2 + 1604*x1^2*x3 + 56*x1^2*x4 - 880*x1*x2^2 - 8024*x1*x2*x3 + 224*x1*x2*x4 - 24120*x1*x3^2 - 7016*x1*x3*x4 - 3264*x1*x4^2 + 400*x1 + 800*x2^3 + 11240*x2^2*x3 - 80*x2^2*x4 + 48240*x2*x3^2 - 6872*x2*x3*x4 - 1168*x2*x4^2 - 800*x2 + 161604*x3^3 + 178488*x3^2*x4 + 145112*x3*x4^2 - 8040*x3 + 37296*x4^3 - 2960*x4,
 12*x1^3 - 148*x1^2*x2 + 56*x1^2*x3 + 576*x1^2*x4 + 576*x1*x2^2 + 224*x1*x2*x3 - 3936*x1*x2*x4 - 3508*x1*x3^2 - 6528*x1*x3*x4 + 4536*x1*x4^2 - 120*x1 - 760*x2^3 - 80*x2^2*x3 + 7928*x2^2*x4 - 3436*x2*x3^2 - 2336*x2*x3*x4 - 28728*x2*x4^2 + 760*x2 + 59496*x3^3 + 145112*x3^2*x4 + 111888*x3*x4^2 - 2960*x3 + 63504*x4^3 - 4968*x4
 ]

 @time begin
 compute_param(F)
 end
=#