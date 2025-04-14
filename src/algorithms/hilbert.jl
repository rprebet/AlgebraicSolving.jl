export affine_hilbert_series, hilbert_series, hilbert_dimension, hilbert_degree, hilbert_polynomial

function hilbert_series(I; variant::Int=0)
    gb = get(I.gb, 0, groebner_basis(I, complete_reduction = true, nr_thrds=Threads.nthreads()))
    lexps = (_drl_lead_exp).(gb)
    return _hilbert_series_mono(lexps, variant=variant)
end

function affine_hilbert_series(I; variant::Int=0)
    gb = get(I.gb, 0, groebner_basis(I, complete_reduction = true, nr_thrds=Threads.nthreads()))
    lexps = (_drl_lead_exp).(homogenize(gb))
    return _hilbert_series_mono(lexps, variant=variant)
end

function hilbert_degree(I)
    return numerator(hilbert_series(I))(1) |> abs
end

function hilbert_dimension(I)
    return denominator(hilbert_series(I)) |> degree
end

function hilbert_polynomial(I)
    H = hilbert_series(I)
    num, dim = numerator(H), degree(denominator(H))
    num = iseven(dim) ? numerator(H) : -numerator(H)
    t = gen(parent(num))
    La = Vector{ZZPolyRingElem}(undef, dim)
    while dim>0
        num, La[dim] = divrem(num, 1-t)
        dim -= 1
    end
    println(La, num)
    Hpolyfct = d->sum(La[i](0)*binomial(i+d, i) for i in 1:length(La))
    dim = degree(denominator(H))
    A, = polynomial_ring(QQ, :d)
    Hpoly = interpolate(A, QQ.(0:dim+1), [QQ(Hpolyfct(d)) for d in 0:dim+1])
    @assert(degree(Hpoly)==dim, "Degree of poly does not match the dimension")
    # Hilbert poly, index of regularity
    return Hpoly, degree(num)+1
end


function _hilbert_series_mono(exps::Vector{Vector{Int}}; variant::Int=0)
    h = _num_hilbert_series_mono(exps, variant=variant)
    t = gen(parent(h))
    return h//(1-t)^length(first(exps))
end

function _num_hilbert_series_mono(exps::Vector{Vector{Int}}; variant::Int=0)
    A, t = polynomial_ring(ZZ, 't')
    r = length(exps)
    r == 0 && return one(A)
    N = length(first(exps))
    ## Base cases ##
    r == 1 && return (1-t^sum(first(exps)))
    supp = findall.(Ref(!iszero), exps)
    pow_supp = findall(s->length(s)==1, supp)
    # If exps is a product of simple powers
    if length(pow_supp) == r
        #println("Simple power")
        return prod(1-t^(exps[i][supp[i][1]]) for i in pow_supp)
    # Only one non-simple power P
    elseif length(pow_supp) == r-1
        #println("Mixed pow")
        inpow = setdiff(eachindex(exps), pow_supp) |> first
        # P has disjoint support with other powers
        if all(iszero(exps[inpow][supp[i][1]]) for i in pow_supp)
            return (1-t^sum(exps[inpow]))*prod(1-t^(exps[i][supp[i][1]]) for i in pow_supp)
        else
            return prod(1-t^(exps[i][supp[i][1]]) for i in pow_supp) - t^sum(exps[inpow]) *
            prod(1-t^(exps[i][supp[i][1]]-exps[inpow][supp[i][1]]) for i in pow_supp)
        end
    end

    # Variable index occuring the most in exps
    counts = sum(x->x .> 0, eachcol(reduce(hcat, exps)))
    ivarmax = argmax(counts)

    ## Splitting recursive cases ##
    # Monomials have disjoint supports
    if counts[ivarmax] == 1
        return prod(1-t^sum(mono) for mono in exps)
    # Heuristic where general splitting is useful
    elseif 8 <= r <= N
        # Finest partition of monomial supports
        LV, h = _monomial_support_partition(exps), one(A)
        rem_mon = collect(1:r)
        for V in LV
            JV, iJV = Vector{Vector{Int}}(), Int[]
            for (k, i) in enumerate(rem_mon)
                mono = exps[i]
                if any(mono[j] != 0 for j in V)
                    push!(iJV, k)
                    push!(JV, mono)
                end
            end
            # Interreduce JV
            JV = [JV[j] for j in eachindex(JV) if
                      !any(all(JV[k] .<= JV[j]) for k in eachindex(JV) if k!=j)]
            h *= _num_hilbert_series_mono(JV, variant=variant)
            # Avoid re-check monomials (LV partitions sat)
            deleteat!(rem_mon, iJV)
        end
        return h
    end

    ## Pivot recursive case ##
    C = 0
    while C < r
        C +=1
        # Exponent of ivarmax in gcd of two random generators
        pivexp = max(1, minimum(mon[ivarmax] for mon in rand(exps, 2)))
        # Compute interreduced generators for (exps):pivot
        sat = Vector{Vector{Int64}}(undef, r)
        for j in 1:r
            sat[j] = [exps[j][1:ivarmax-1]; max(exps[j][ivarmax]-pivexp, 0); exps[j][ivarmax+1:end]]
        end
        sat = [sat[j] for j in eachindex(sat) if !iszero(sat[j]) &&
                         !any(all(sat[k] .<= sat[j]) for k in eachindex(sat) if k!=j)]
        isempty(sat) && continue # We must split the ideal

        # Interreduce exps + pivot
        filter!(e->(pivexp > e[ivarmax]), exps)
        push!(exps,[zeros(Int64,ivarmax-1); pivexp; zeros(Int64,N-ivarmax)])

        a = _num_hilbert_series_mono(exps)
        b = _num_hilbert_series_mono(sat)
        return a+b*t^pivexp
    end
    error("Impossible to find a suitable pivot")
end

function _monomial_support_partition(L::Vector{Vector{Int}})
    # Build adjacency graph: connect variables that appear together in a monomial
    n = length(first(L))
    adj = [Set{Int}() for _ in 1:n]
    supp_tot = (!).(trues(n))
    for mono in L
        support = findall(!=(0), mono)
        for i in support
            supp_tot[i] = supp_tot[i] || true
        end
        for i in support, j in support
            if i != j
                push!(adj[i], j)
            end
        end
    end

    # DFS to extract connected components
    visited = falses(n)
    components = Vector{Vector{Int}}()

    function dfs(u, comp)
        visited[u] = true
        push!(comp, u)
        for v in adj[u]
            if !visited[v]
                dfs(v, comp)
            end
        end
    end

    for v in [i for i in 1:n if supp_tot[i]]
        if !visited[v]
            comp = Int[]
            dfs(v, comp)
            push!(components, comp)
        end
    end

    return components
end