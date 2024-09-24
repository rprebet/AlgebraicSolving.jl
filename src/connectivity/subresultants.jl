# Some tools
function poly_to_array(P::MPolyRingElem)
    # return P as two lists ci and vi of resp. coeffs and exponents
    return [collect(coefficients(P)), collect(exponent_vectors(P))]
end

function rem_ind(L,i)
    # remove element at index i
    return [ L[1:i-1]; L[i+1:end] ]
end

function add_ind(L,i,x)
    # add element x at index i
    return [ L[1:i-1]; [x]; L[i:end] ]
end

function parray_asvar(LP, idx)
    # takes the above representation of a poly P  
    # and outputs a representation of P seen as a univariate poly in the i-th variable
    deg = sort(LP[2], by= x->x[idx])[end][idx]
    NLP = [ [[],[]] for _ in 1:deg+1 ]
    for i in eachindex(LP[2])
        di = LP[2][i][idx]+1
        push!(NLP[di][1], LP[1][i])
        push!(NLP[di][2], rem_ind(LP[2][i], idx)) 
    end
    return NLP
end

function array_to_poly(LP, A)
    C = MPolyBuildCtx(A)
    R = base_ring(A)
    for i in eachindex(LP[1])
        push_term!(C, R(LP[1][i]), LP[2][i]);
    end
    return finish(C)
end

# Univariate resultant
function subresultants(P::PolyRingElem{T}, Q::PolyRingElem{T}) where T <: RingElement
    if degree(P) < degree(Q)
        P, Q = Q, P
    end
    S = [Q]
    s = leading_coefficient(Q)^(degree(P)-degree(Q))
    A = Q
    B = pseudorem(P,-Q)
    ring = parent(P)
    while true
        d = degree(A)
        e = degree(B)
        #println("($d,$e)")
        if iszero(B)
            return S
        end
        pushfirst!(S, copy(B))
        delta = d - e
        if delta > 1
            if length(S) > 1
                n = degree(S[2]) - degree(S[1]) - 1
                if n == 0
                    C = copy(S[1])
                else
                    x = leading_coefficient(S[1])
                    y = leading_coefficient(S[2])
                    a = 1 << (length(bits(ZZ(n))) - 1)
                    c = x
                    n = n - a
                    while a > 1
                        a /= 2
                        c = c^2 / y
                        if n >= a
                            c = c * x / y
                            n = n - a
                        end
                    end
                    C = c * S[1] / y
                end
            else
                C = leading_coefficient(B)^(delta-1) * B / s^(delta-1)
            end
            pushfirst!(S, copy(C))
        else
            C = copy(B)
        end
        if e == 0
            return S
        end
        #@time 
        B = pseudorem(A,-B) / (s^delta * leading_coefficient(A))
        A = copy(C)
        s = leading_coefficient(A)
    end
end

#Bivariate subresultants
function subresultants(P::MPolyRingElem{T}, Q::MPolyRingElem{T}, idx; list=false) where T <: RingElement
    LPQ = map(poly_to_array, [P,Q])
    ULPQ = [ parray_asvar(lpq, idx) for lpq in LPQ ]

    R, x = polynomial_ring(base_ring(parent(P)), "x")
    S, y = polynomial_ring(R, "y")

    UP, UQ = [ S([ array_to_poly(l, R) for l in lpq ]) for lpq in ULPQ ]

    sr = subresultants(UP, UQ)

    # Get it back to initial polynomial ring
    Lsr = [ [collect(coefficients(csr)) for csr in coefficients(sr)] for sr in sr]
    newsr = []
    if list
        for lsr in Lsr
            mlsr = []
            for i in 1:length(lsr)
                tmp = [[],[]]
                for j in 1:length(lsr[i])
                    push!(tmp[1], lsr[i][j])
                    push!(tmp[2],add_ind([j-1],idx, 0))
                end
                push!(mlsr, tmp)
            end
            mlsr = [ array_to_poly(tmp, parent(P)) for tmp in mlsr ]
            push!(newsr, mlsr)
        end
    else
        for lsr in Lsr
            mlsr = [[],[]]
            for i in 1:length(lsr)
                for j in 1:length(lsr[i])
                    push!(mlsr[1], lsr[i][j])
                    push!(mlsr[2],add_ind([j-1],idx, i-1))
                end
            end
            push!(newsr, array_to_poly(mlsr, parent(P)))
        end
    end

    return newsr
end

function interp_subresultants(P::MPolyRingElem{T}, Q::MPolyRingElem{T}, iv; list=false) where T <: RingElement
	R = parent(P)
    A = coefficient_ring(R)
    iu = iv%2+1
    Vu, Vv = [:x, :y]
	# Modern computer Algebra Th 6.51.
	degs = degrees(P)
	n, d = degs[iv], degs[iu]
	degsmax = [ (2*n-1-2*k)*d for k=0:d-1 ] 
	#Dmax = degsmax[end]>>1+1
	##
    B,t =polynomial_ring(A, Vv)
    xs = A.(0:degsmax[1])
    Srev = Vector{Vector{typeof(t)}}(undef, length(xs))
    evPQ = [ [B(coefficients_of_univariate(evaluate(pq,[iu],[xx]))) for pq in [P,Q] ] for xx in xs]
    Threads.@threads for i in 1:length(xs)
	    Srev[i] = subresultants(evPQ[i][1], evPQ[i][2])
    end
    Srev = [ [ collect(coefficients(ss)) for ss in s] for s in Srev ]
	C, = Nemo.polynomial_ring(A,Vu)
	Sr = [ [ C(0) for _=1:i] for i=1:n ]
    Lys = [ [ [ (try A(Srev[l][i][j]) catch; A(0) end) for l in eachindex(Srev) ] for j in 1:i ] for i in 1:n ]
	for i = 1:n
        # For the i-th subresultant
        Threads.@threads for j in 1:i
            # Each coefficient is interpolated
            Sr[i][j] = interpolate(C, xs, Lys[i][j])
        end
	end
	return Sr
end

function mmod_subresultants(P::QQMPolyRingElem, Q::QQMPolyRingElem, idx; list=false, n_threads=Threads.nthreads(), n_ssr=-1,v=0)
    PZ, QZ = change_coefficient_ring.(Ref(ZZ), int_coeffs.([P,Q]))
    sr = mmod_subresultants(PZ, QZ, idx; list=list, n_threads=n_threads, n_ssr=n_ssr, v=v)
    return sr
    if list
        return [ change_coefficient_ring.(Ref(QQ),srp) for srp in sr ]
    else
        return change_coefficient_ring.(Ref(QQ),sr)
    end
end

function mmod_subresultants(P::ZZMPolyRingElem, Q::ZZMPolyRingElem, idx; list=false, n_threads=Threads.nthreads(), n_ssr=0, v=0)
    prim = ZZ(1)<<(61)
    L1, primprod = [], ZZ(1)
    Ltemp = Vector{Vector{Vector{Vector{ZZRingElem}}}}(undef,n_threads+1)
    lcpq =  [ evaluate(leading_coefficient(l, 2),[0,0]) for l in [P,Q] ]
    compt=n_threads
    while true
        v > 0 && print("$compt primes ")
        LP, LQ, Lprim = MPolyRingElem[], MPolyRingElem[], ZZRingElem[]
        while length(Lprim) < n_threads
            prim = next_prime(prim)
            if !any(getindex.(divides.(lcpq, Ref(prim)),1)) 
                Pprim, Qprim = [ change_coefficient_ring(residue_ring(ZZ, prim)[1], poly) for poly in [P,Q] ]
                #Pprim, Qprim = [ change_coefficient_ring(GF(prim), poly) for poly in [P,Q] ]
                push!.([LP,LQ,Lprim], [Pprim,Qprim,prim])
            else
                return prim, lcpq
            end
        end
        Threads.@threads for j in 1:n_threads
            sr = subresultants(LP[j], LQ[j], idx, list=true)
            Ltemp[j] = [ [ lift.(Ref(ZZ), coefficients_of_univariate(srii)) for srii in sri ] for sri in sr ]
            #sr = interp_subresultants(LP[j], LQ[j], idx)[1][1]
            #Ltemp[j] = lift.(Ref(ZZ), coefficients(sr))
        end

        if n_ssr > 0
            for j in 1:n_threads
                for k in eachindex(Ltemp[j])
                    Ltemp[j][k] = Ltemp[j][k][end-min(k,n_ssr)+1:end]
                end
            end
        end
        
        if L1 != []
            Ltemp[n_threads+1] = copy(L1)
            push!(Lprim, primprod)
        end
        #@assert allequal(length.(Ltemp[1:nthreads])) "Specialization problem"  
        Nssr = n_ssr<0 ? length(last(first(Ltemp))) : n_ssr 
        #println(Ltemp)
        L2 = deepcopy(L1)
        L1 = Vector{Vector{Vector{ZZRingElem}}}(undef, length(first(Ltemp)))
        Threads.@threads for i in eachindex(L1)
            L1[i] = Vector{Vector{ZZRingElem}}(undef, min(i, Nssr))
            Threads.@threads for j in eachindex(L1[i])
                #println("($i,$j)")
                L1[i][j] = Vector{ZZRingElem}(undef, length(Ltemp[1][i][j]))
                Threads.@threads for k in eachindex(L1[i][j])
                    cfs = [ try Ltemp[p][i][j][k] catch; ZZ(0) end for p in eachindex(Lprim) ]
                    L1[i][j][k] = crt(cfs, Lprim, true)
                end
            end
        end
        #println(L1)
        primprod = prod(Lprim)
        compt+=n_threads
        v > 0 && println(ndigits(primprod, 10), " digits")
        L1 != L2 || break
        #compt < 20 || break
    end

    uvar = parent(P).S[(idx+1)%2]
    A, = polynomial_ring(ZZ, uvar)
    sr = [ [ change_ringvar(A(ss),[uvar]) for ss in s ] for s in L1 ]
    return sr
end

function fact_gcd(delta::T, LP::Vector{T}) where (T <:PolyRingElem)
    # delta : poly to factor w.r.t the polynomials in LP
    Lphi = [gcd(delta, LP[1])]
    Ldelta = [delta/Lphi[1]]
    i = 2
    while degree(Lphi[end])>0
        push!(Lphi, gcd(Lphi[i-1], LP[i]))
        push!(Ldelta, Lphi[i-1]/Lphi[i])
        i+=1
    end
    return filter(s->degree(s[2])>0, Dict(enumerate(Ldelta)))
end

function fact_gcd(delta::T, LP::Vector{T}) where (T <:MPolyRingElem)
    @assert is_univariate(delta) && all(is_univariate.(LP)) "Not univariate polynomial"
    R = parent(delta)
    A, = polynomial_ring(coefficient_ring(R), first(R.S))
    out = fact_gcd(A(coefficients_of_univariate(delta)), A.(coefficients_of_univariate.(LP)))
    return Dict([ (i, change_ringvar(f, [first(R.S)]) ) for (i,f) in out ])
end

function param_crit_split(f)
    # Compute subresultants and factor the first subresultant
    if total_degree(f) > 30
        @time sr = mmod_subresultants(f, derivative(f, 2), 2, list=true, n_ssr=2)
    else
        @time sr = subresultants(f, derivative(f, 2), 2, list=true)
    end
    sqr = collect(factor_squarefree(sr[1][1]))

    # Filter out factors with zero multiplicity and sort by multiplicity
    sqr = sort!(filter(t -> t[2] > 0, sqr), by = t -> t[2])
    # Group factors by multiplicity
    sqrmult = unique(getindex.(sqr, 2))
    group_sqr = Dict(m => [r[1] for r in sqr if r[2] == m] for m in sqrmult)
    
    # Initalization
    singmult = filter(p->p*(p-1)<=sqrmult[end], 2:sqrmult[end])
    param_crit = Dict(p => [QQMPolyRingElem[], -sr[p][end-1], (p-1)*sr[p][end]]  for p in singmult)
    lsr = [ sr[p][end] for p in singmult ]

    # Critical points : multiplicity 1 in res
    (1 in sqrmult) && push!(param_crit, 1=>[group_sqr[1], -sr[2][end-1], sr[2][end]])
    # Nodes : multiplicity 2 in res
    (2 in sqrmult) && append!(param_crit[2][1], group_sqr[2])
    # Other sing
    filter!(m->!(m in [1,2]), sqrmult)
    #TODO: simpler criterion for mult=p*(p-1)?
    Ld = Vector{Vector{Dict{Int, QQMPolyRingElem}}}(undef, length(sqrmult))
    #Threads.@threads 
    for k in eachindex(sqrmult)
        for l in eachindex(group_sqr[sqrmult[k]])
            Ld[k][l] = fact_gcd(group_sqr[sqrmult[k]][l], lsr)
        end
    end
    for k in eachindex(Ld)
        for l in eachindex(Ld[k])
            for (i, dji) in Ld[k][l]
                push!(param_crit[i+1][1], dji)
            end
        end
    end

    filter!(p->length(p[2][1])>0, param_crit)
    return param_crit
end