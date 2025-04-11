include("dimension.jl")

export affine_hilbert_series, hilbert_series, hilbert_series_mono, hilbert_dimension, hilbert_degree

function hilbert_series_mono(exps::Vector{Vector{Int}}; variant::Int=0)
    h = _num_hilbert_series_mono(exps, variant=variant)
    t = gen(parent(h))
    return h//(1-t)^length(first(exps))
end

function _num_hilbert_series_mono(exps::Vector{Vector{Int}}; variant::Int=0)
    r = sort!(exps, by=reverse) |> length

    A, t = polynomial_ring(ZZ, 't')
    h = iszero(r) ? one(A) : 1-t^(sum(exps[1]))
    for i in 2:r
        # Compute generators for (x^a1,...,x^a{i-1}):x^ai
        sat = [ [ max(a[j]-exps[i][j], 0) for j in eachindex(a)]
                                                for a in exps[1:i-1] ]
        # Reduce to minimal generators
        sat = [sat[j] for j in eachindex(sat) if
                      !any(all(sat[k] .<= sat[j]) for k in eachindex(sat) if k!=j)]

        if variant==0
            # Trivial partition of monomial supports
            nolin_sat = [u for u in sat if sum(u)>1]
            hsat = (1-t)^(length(sat)-length(nolin_sat))*_num_hilbert_series_mono(nolin_sat, variant=variant)
        elseif variant==1
            # Finest partition of monomial supports
            LV, hsat = _monomial_support_partition(sat), one(A)
            rem_mon = collect(1:length(sat))
            for V in LV
                JV, iJV = Vector{Vector{Int}}(), Int[]
                for (k, i) in enumerate(rem_mon)
                    mono = sat[i]
                    if any(mono[j] != 0 for j in V)
                        push!(iJV, k)
                        push!(JV, mono)
                    end
                end
                hsat *= _num_hilbert_series_mono(JV, variant=variant)
                # Avoid re-check monomials (LV partitions sat)
                deleteat!(rem_mon, iJV)
            end
        else
            # No partition of monomial supports
            hsat = _num_hilbert_series_mono(sat, variant=variant)
        end
        h -= t^(sum(exps[i]))*hsat
    end
    return h
end

function hilbert_series(I; variant::Int=0)
    gb = get(I.gb, 0, groebner_basis(I, complete_reduction = true))
    lexps = (_drl_lead_exp).(gb)
    return hilbert_series_mono(lexps, variant=variant)
end

function affine_hilbert_series(I; variant::Int=0)
    gb = get(I.gb, 0, groebner_basis(I, complete_reduction = true))
    lexps = (_drl_lead_exp).(homogenize(gb))
    return hilbert_series_mono(lexps, variant=variant)
end

function homogenize(F::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), ["x$i" for i in 1:nvars(R)+1], internal_ordering=:degrevlex)
    res = typeof(first(F))[]
    for f in F
        ctx = MPolyBuildCtx(S)
        d = total_degree(f)
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            enew = push!(e, d - sum(e))
            push_term!(ctx, c, e)
        end
        push!(res, finish(ctx))
    end
    return res
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

function hilbert_degree(I)
    return numerator(hilbert_series(I))(1) |> abs
end

function hilbert_dimension(I)
    return denominator(hilbert_series(I)) |> degree
end