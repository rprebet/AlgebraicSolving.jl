@doc Markdown.doc"""
    curve_rational_parametrization(I::Ideal{<:MPolyRingElem}, <keyword arguments>)

Given a **radical** ideal `I` with solution set X being of dimension 1 over the complex numbers,
return a rational curve parametrization of the one-dimensional irreducible components of X.

**Important**: In the output, the variables x and y correspond respectively to the last and second-to-last entries of the vars attribute.

**Note**: At the moment only QQ is supported as ground field. If the dimension of the ideal
is not one an ErrorException is thrown.

# Arguments
- `I::Ideal{<:QQMPolyRingElem}`: input generators.
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).
- `cfs_lfs::Vector{Vector{ZZRingElem}} = []`: coefficients for the above linear forms
- `nr_thrds::Int=1`: number of threads for msolve
- `extra_form::Bool=false`: if `cfs_lfs` is not given, search for a 3rd linear
form in addition to the usual "x, y" pair, found (or, if `cfs_lfs` has 3
entries, verified) AFTER them, distinct from both, and never referenced by
either.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x1,x2,x3) = polynomial_ring(QQ, ["x1","x2","x3"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> I = Ideal([x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1])
QQMPolyRingElem[x1 + 2*x2 + 2*x3 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2]

julia> curve_rational_parametrization(I)
AlgebraicSolving.CurveRationalParametrization([:x1, :x2, :x3, :_Z2, :_Z1], Vector{ZZRingElem}[[0, 0, 1, 0, -1], [0, 1, 0, -1, 0]], x^2 + 4//3*x*y - 1//3*x + y^2 - 1//3*y, 4//3*x + 2*y - 1//3, QQMPolyRingElem[4//3*x^2 - 4//3*x*y + 2//3*x + 4//3*y - 1//3, -2*x^2 - 4//3*x*y + 2//3*x + 1//3*y, 4//3*x^2 + 2*x*y - 1//3*x])

julia> curve_rational_parametrization(I, cfs_lfs=[[-3,2,2], [1,4,-3]])
AlgebraicSolving.CurveRationalParametrization([:x1, :x2, :x3, :_Z2, :_Z1], Vector{ZZRingElem}[[-3, 2, 2, 0, -1], [1, 4, -3, -1, 0]], 127//128*x^2 + 3//8*x*y + 161//64*x + y^2 - 7//8*y - 49//128, 3//8*x + 2*y - 7//8, QQMPolyRingElem[-3//32*x^2 - 1//2*x*y + 5//16*x + 1//2*y - 7//32, -1//4*x^2 + 1//8*x*y - 3//4*x + 3//8*y, 19//64*x^2 + 1//8*x*y + 25//32*x + 3//8*y - 21//64])

julia> curve_rational_parametrization(I, cfs_lfs=[[-3,2,2,-1,-2], [1,4,-3,2,-1]])
AlgebraicSolving.CurveRationalParametrization([:x1, :x2, :x3, :_Z2, :_Z1], Vector{ZZRingElem}[[-3, 2, 2, -1, -2], [1, 4, -3, 2, -1]], 244//181*x^2 - 148//543*x*y + 532//543*x + y^2 + 182//181*y - 49//543, -148//543*x + 2*y + 182//181, QQMPolyRingElem[440//543*x^2 - 580//543*x*y - 44//543*x + 136//181*y + 112//543, 80//181*x^2 + 320//543*x*y + 122//181*x + 81//181*y + 49//543, -460//543*x^2 - 10//181*x*y - 418//543*x + 32//181*y + 56//181])
```
"""
function curve_rational_parametrization(I::Ideal{<:QQMPolyRingElem}; kwargs...)
    return curve_rational_parametrization([I]; kwargs...)[1]
end

@doc Markdown.doc"""
    curve_rational_parametrization(curves::Vector{<:Ideal{<:QQMPolyRingElem}}, <keyword arguments>)

Vector-of-curves version of `curve_rational_parametrization`: finds a SINGLE
set of linear forms (`cfs_lfs`, either given or searched for jointly, as
with a single curve) valid for EVERY curve in `curves` at once.
"""
function curve_rational_parametrization(
        curves::Vector{<:Ideal{<:QQMPolyRingElem}};                                    # input generators, sharing the SAME linear forms
        info_level::Int=0,                                                             # info level for print outs
        cfs_lfs::Union{Vector{<:Vector{<:Union{Int,ZZRingElem}}}, Nothing} = nothing,  # coeffs of linear forms
        nr_thrds::Int=1,                                                               # number of threads (msolve)
        check_cfs::Bool = true,
        extra_form::Bool = false                                                      # search a 3rd, distinct form (see `_add_genvars_extra`)
    )
    @assert !isempty(curves) "Must provide at least one curve."
    @assert all(nvars(parent(I)) >= 2 for I in curves) "Each curve must be defined in a ring with at least 2 variables"
    if !isnothing(cfs_lfs)
        @assert length(cfs_lfs) >= 2 "When specified, at least two linear forms must be provided"
        cfs_lfs = Vector{ZZRingElem}[[ ZZRingElem(c) for c in cfs_lf] for cfs_lf in cfs_lfs] # Convert Int64 into ZZRingElem
    end

    info_level > 0 && println("Compute generic linear forms...")
    n_gen = isnothing(cfs_lfs) ? (extra_form ? 3 : 2) : length(cfs_lfs)
    news, cfs_lfs = _add_genvars(curves, n_gen, cfs_lfs, check_cfs = check_cfs)

    return [ _curve_rational_parametrization_ext(I, Inew, cfs_lfs; info_level, nr_thrds)
             for (I, Inew) in zip(curves, news) ]
end

# Given the original ideal `I` and its already-embedded/pinned version
# `Inew` (`n_gen` new variables, with the linear-form equations
# already added, whose coeficients are in `cfs_lfs`), computes and
# caches (`I.rat_param`) the actual bivariate rational parametrization
# of `I`.
# This is the actual core function
function _curve_rational_parametrization_ext(I::Ideal, Inew::Ideal, cfs_lfs; info_level::Int=0, nr_thrds::Int=1)
    if Inew.dim == -1
        T = polynomial_ring(QQ, [:x, :y])[1]
        I.dim = -1
        I.rat_param = CurveRationalParametrization(Symbol[], Vector{ZZRingElem}[], T(-1), T(-1), QQMPolyRingElem[])
        return I.rat_param
    end
    @assert Inew.dim == 1 "Input ideal(s) must define a curve or an empty set"

    R = parent(Inew)
    N = nvars(R)
    DEG, F = Inew.deg, Inew.gens

    # Compute DEG+2 evaluations of x in the param (whose total deg is bounded by DEG)
    PARAM  = Vector{Vector{QQPolyRingElem}}(undef, DEG+2)
    _values = Vector{ZZRingElem}(undef, DEG+2)

    i = 1
    free_ind = collect(1:DEG+2)
    used_ind = falses(DEG+2)
    lc = nothing

    while length(free_ind) > 0
        if i > 2 * (DEG + 2)
            error("Too many bad specializations. Check radicality and generic linear forms.")
        end

        # Determine values to evaluate at to keep bitsize low
        curr_values = ZZ.([-(i - 1 + (length(free_ind) + 1) ÷ 2):-i; i:(i - 1 + length(free_ind) ÷ 2)])
        LFeval = _evalvar(F, N, curr_values)

        # Compute parametrization of each evaluation
        Lr = Vector{RationalParametrization}(undef, length(free_ind))
        for j in eachindex(free_ind)
            info_level > 0 && print("Evaluated parametrizations: $(j)/$(length(free_ind))\r")
            Lr[j] = rational_parametrization(Ideal(LFeval[j]), nr_thrds=nr_thrds, info_level=0)

            # Specialization checks: same vars order, generic degree
            if Lr[j].vars == symbols(R)[1:N-1] && degree(Lr[j].elim) == DEG
                if isnothing(lc)
                    lc = leading_coefficient(Lr[j].elim)
                    rr = vcat([Lr[j].elim, Lr[j].denom], Lr[j].param)
                else
                    # Adjust when the rat_param is multiplied by some constant factor
                    fact = lc / leading_coefficient(Lr[j].elim)
                    rr = vcat([Lr[j].elim * fact, Lr[j].denom * fact], Lr[j].param .* fact)
                end
                PARAM[free_ind[j]] = rr
                _values[free_ind[j]] = curr_values[j]
                used_ind[j] = true
            end
        end

        # Update range, free indices and used indices
        i += length(free_ind)
        free_ind = free_ind[.!used_ind]
        used_ind = falses(length(free_ind))

        if info_level * length(free_ind) != 0
            println("bad specialization(s): ", curr_values[free_ind])
        end
    end

    # Interpolate each coefficient of each poly in the param
    T, = polynomial_ring(QQ, [:x, :y])
    A, = polynomial_ring(QQ)

    POLY_PARAM = Vector{QQMPolyRingElem}(undef, N)
    for count in 1:N
        info_level > 0 && print("Interpolate parametrizations: $count/$N\r")
        COEFFS = Vector{QQPolyRingElem}(undef, DEG + 1)

        for deg in 0:DEG
            _evals = [coeff(PARAM[i][count], deg) for i in eachindex(PARAM)]
            # Remove denominators for faster interpolation with FLINT
            # TODO: remove dens mult when interface's ready in Nemo
            den = reduce(lcm, (denominator(e) for e in _evals), init=one(ZZ))
            scaled_evals = [ZZ(_evals[i] * den) for i in eachindex(_evals)]
            COEFFS[deg + 1] = interpolate(A, _values, scaled_evals) / (lc * den)
        end

        ctx = MPolyBuildCtx(T)
        for (idx, c) in enumerate(COEFFS)
            for (j, coef) in enumerate(coefficients(c))
                !iszero(coef) && push_term!(ctx, coef, [j - 1, idx - 1])
            end
        end
        POLY_PARAM[count] = finish(ctx)
    end
    info_level > 0 && println()

    I.deg, I.dim = Inew.deg, Inew.dim
    I.rat_param = CurveRationalParametrization(symbols(R), cfs_lfs, POLY_PARAM[1], POLY_PARAM[2], POLY_PARAM[3:end])
    return I.rat_param
end


# Embeds I into a ring with n_gen extra free variables
function _embed_extra_vars(
    I::Ideal{T},
    n_gen::Int,
    genS::Vector{Symbol} = Symbol[]
) where T <: MPolyRingElem
    F = I.gens
    R = parent(I)
    K, n = base_ring(R), nvars(R)

    # Add new variables (reverse index order)
    @assert isempty(genS) || (length(unique(genS)) == length(genS) == n_gen) "Bad provided names for generic variables"
    genS = isempty(genS) ? [Symbol("_Z", i) for i in n_gen:-1:1] : genS
    newS = vcat(symbols(R), genS)
    R_ext, all_vars = polynomial_ring(K, newS)

    # Inject F in this new ring efficiently by reconstruction
    F_ext = Vector{T}(undef, length(F))
    ctx = MPolyBuildCtx(R_ext)
    new_e = zeros(Int, n + n_gen) # Pre-allocated buffer

    for i in eachindex(F)
        for (e, c) in zip(exponent_vectors(F[i]), coefficients(F[i]))
            new_e[1:n] .= e
            push_term!(ctx, c, new_e)
        end
        F_ext[i] = finish(ctx)
    end

    return Ideal(F_ext)
end

# Embeds I with n_gen new variables (`_embed_extra_vars`) and adds the
# n_gen equations Li(X) - Zi = 0 pinning each new variable to a generic linear
# form: either given via `cfs_lfs`, or found via incremental process.
# cfs_lfs is checked to be generic iff `check_cfs==true`.
function _add_genvars(
    ideals::Vector{Ideal{T}},
    n_gen::Int,
    cfs_lfs::Union{Vector{<:Vector{<:RingElem}}, Nothing},
    genS::Vector{Symbol} = Symbol[];
    check_cfs::Bool = true,
    excluded::Vector{Vector{ZZRingElem}} = Vector{ZZRingElem}[],
    points::Vector{<:Ideal} = Ideal[]
) where T <:MPolyRingElem
    # `_add_genvars_extra` splits into "the 2 standard ones, then 1 extra" instead
    if n_gen == 3 && check_cfs && isempty(genS)
        return _add_genvars_extra(ideals, cfs_lfs; points=points)
    end

    n = nvars(parent(ideals[1]))
    ideals_ext = [_embed_extra_vars(I, n_gen, genS) for I in ideals]

    # Find generic linear forms
    if !isnothing(cfs_lfs)
        @assert length(cfs_lfs) == n_gen "Expected $n_gen linear forms, got $(length(cfs_lfs))"
        @assert all(length(c) in [n, n + n_gen] for c in cfs_lfs) "Linear forms must have $n or $(n + n_gen) coefficients"
        if length(first(cfs_lfs)) == n
            cfs_lfs = [vcat(c, [-ZZ(j == i) for j in n_gen:-1:1]) for (i, c) in enumerate(cfs_lfs)]
        end
    end

    if check_cfs
        (DEGs, DIMs), cfs_lfs = _find_generic_linear_forms(ideals_ext, n_gen, cfs_lfs; excluded=excluded, points=points)
    else
        DEGs = Vector{Int64}(undef, length(ideals))
        DIMs = Vector{Int64}(undef, length(ideals))
        uZ = one(ZZRingElem)
        for (i, I) in enumerate(ideals_ext)
            F = I.gens
            lucky_prime = first(_generate_lucky_primes(F, uZ<<30, (uZ<<31)-1, 1))
            Itest = Ideal(change_base_ring.(Ref(GF(lucky_prime)), F))
            DEGs[i], DIMs[i] = hilbert_degree(Itest), dimension(Itest)
        end
    end

    # Add equations Li(X) - Zi = 0
    all_vars = gens(parent(ideals_ext[1]))
    lf = [transpose(c) * all_vars for c in cfs_lfs]
    ideals_new = [ Ideal(vcat(I.gens, lf)) for I in ideals_ext ]
    for i in 1:length(ideals_new)
        ideals_new[i].deg = DEGs[i]
        ideals_new[i].dim = max(DIMs[i] - n_gen, -1)
    end

    return ideals_new, cfs_lfs
end

# Single-ideal version
function _add_genvars(I::Ideal{T}, args...; kwargs...) where T <:MPolyRingElem
    Inew, cfs_lfs = _add_genvars([I], args...; kwargs...)
    return Inew[1], cfs_lfs
end

# Handles n_gen == 3 specially: the first two forms are found (or
# verified) exactly as the standard n_gen = 2 case and then the third
# is chosen (or verified) **distinct** using the exact same genericity
# test as the very first form
function _add_genvars_extra(
        ideals::Vector{<:Ideal{T}},
        cfs_lfs::Union{Vector{<:Vector{<:RingElem}}, Nothing};
        points::Vector{<:Ideal} = Ideal[]
    ) where T <: MPolyRingElem
    n = nvars(parent(ideals[1]))

    # Normalize to raw (length n) coefficients
    if !isnothing(cfs_lfs)
        @assert length(cfs_lfs) == 3 "Expected 3 linear forms, got $(length(cfs_lfs))"
        @assert all(length(c) in [n, n + 3] for c in cfs_lfs) "Linear forms must have $n or $(n + 3) coefficients"
        cfs_lfs = [ ZZRingElem.(c[1:n]) for c in cfs_lfs ]
    end

    # First two forms
    cfs_lfs12 = isnothing(cfs_lfs) ? nothing : cfs_lfs[1:2]
    ideals_xy, cfs_lfs12 = _add_genvars(ideals, 2, cfs_lfs12; check_cfs=true, points=points)

    # Third one
    cfs_lfs3 = isnothing(cfs_lfs) ? nothing : [cfs_lfs[3]]
    excluded = Vector{ZZRingElem}[ c[1:n] for c in cfs_lfs12 ] # distinct from the first two
    _, cfs_lfs3 = _add_genvars(ideals, 1, cfs_lfs3, [:_Z3]; check_cfs=true, excluded=excluded)
    z = only(cfs_lfs3)[1:n]

    # Re-embed the ORIGINAL ideals with all 3 (now validated) forms
    cfs_lfs_out = [cfs_lfs12[1][1:n], cfs_lfs12[2][1:n], z]
    ideals_ext = [ _embed_extra_vars(I, 3, [:_Z3, :_Z2, :_Z1]) for I in ideals ]
    all_vars = gens(parent(ideals_ext[1]))
    padded = [ vcat(c, [-ZZ(j == i) for j in 3:-1:1]) for (i, c) in enumerate(cfs_lfs_out) ]
    lf = [ transpose(c) * all_vars for c in padded ]
    ideals_new = [ Ideal(vcat(I.gens, lf)) for I in ideals_ext ]
    for i in eachindex(ideals_new)
        ideals_new[i].deg, ideals_new[i].dim = ideals_xy[i].deg, ideals_xy[i].dim
    end

    return ideals_new, padded
end


# Computes/tests n_gen sequential generic linear forms so that the last
# n_gen variables are in generic position (w.r.t. each other, starting from
# the last one) SIMULTANEOUSLY for every ideal in `ideals`
# -- each of which must already live in a ring with n_gen extra (free,
# unconstrained) variables appended (see `_embed_extra_vars`)
#-- ALL IN THE SAME RING (same total number of variables).

# Two modes:
# - search (`cfs_lfs === nothing`): draws candidates from `_candidate_stream`
#   and accepts the first one valid for EVERY ideal, one linear form at a
#   time;
# - verification (`cfs_lfs` given): checks that the already-fixed forms are
#   generic enough for EVERY ideal.
function _find_generic_linear_forms(
        ideals::Vector{<:Ideal{<:MPolyRingElem}},
        n_gen::Int,
        cfs_lfs::Union{Vector{Vector{ZZRingElem}}, Nothing} = nothing;
        excluded::Vector{Vector{ZZRingElem}} = Vector{ZZRingElem}[],
        points::Vector{<:Ideal} = Ideal[]
    )
    n = nvars(parent(ideals[1]))
    n_nogen = n - n_gen
    uZ = one(ZZRingElem)
    is_verif = !isnothing(cfs_lfs)
    max_iter = is_verif ? 1 : 10000

    # Per-ideal context: its own degree/dimension, bifurcation bound
    # and a mutable "probe" system (`:F`) that accumulates one
    # evaluation constraint per accepted form.
    ctxs = map(ideals) do I
        F = I.gens
        prime, Itest, DEG, DIM = _reference_deg_dim(F)

        if DIM > 0
            @assert DIM == n_gen + 1 "Input ideal(s) must define curves or empty sets"
            @assert 2*n_gen < DIM + 2 "Too many generic linear forms asked > dim + 1"
            _radical_evidence(Itest, DIM, DEG) || error("Input ideal(s) are not radical")
        end
        # Bound on bifurcation set degree (e.g., Jelonek & Kurdyka, 2005),
        # with a 2^20 margin
        max_deg = maximum(f -> total_degree(f), F; init=1)
        bif_bound = uZ << (n * ceil(Int, log2(max_deg)) + 21)
        # The data
        Dict(:F0 => F, :F => copy(F), :lucky_prime => prime,
             :DEG => DEG, :DIM => DIM, :bif_bound => bif_bound)
    end

    # 0-dim ideals (control/intersection points)
    pctxs = Tuple{Ideal, Int, ZZRingElem}[]
    for I in points
        q, Imod, DEGp, DIMp = _reference_deg_dim(I.gens)
        @assert DIMp <= 0 "Input points must be finitely many"
        if DIMp == 0
            _radical_evidence(Imod, 0, DEGp) || error("Input point ideal(s) are not radical")
            push!(pctxs, (_embed_extra_vars(I, n_gen), DEGp, q))
        end
    end

    cfs_lfs_out = Vector{Vector{ZZRingElem}}()
    # A candidate rejected for one k is rejected for further k
    stream = _candidate_stream(n_nogen)

    for k in 1:n_gen
        vals = Vector{Vector{ZZRingElem}}(undef, length(ctxs))
        # checks that cand is a valid lf for all ideals, and for the points
        function valid_for_all(cand)
            cvals = map(_draw_val, ctxs)
            ok = all(eachindex(ctxs)) do idx
                     c = ctxs[idx]
                     _is_valid_linear_form(c[:F], cand, cvals[idx], c[:DEG], c[:DIM], k, c[:F0], cfs_lfs_out, c[:lucky_prime])
                 end &&
                 (k != 1 ||
                  all(pc -> _is_valid_points_form(pc[1], pc[2], vcat(cfs_lfs_out, [cand]), pc[3]), pctxs))
            ok && (vals .= cvals)
            return ok
        end

        if is_verif
            coeffs = cfs_lfs[k]
            coeffs[1:n_nogen] in excluded && error("Provided linear form number $k is not distinct from an excluded one.")
            valid_for_all(coeffs) || error("Provided linear form number $k, failed the genericity test.")
        else
            coeffs = nothing
            for (attempt, base) in enumerate(stream)
                attempt > max_iter && break
                base in excluded && continue # no excluded choice
                cand = vcat(base, [-ZZ(j == k) for j in n_gen:-1:1])
                cand in cfs_lfs_out && continue # no redundant choice
                if valid_for_all(cand)
                    coeffs = cand
                    break
                end
            end
            isnothing(coeffs) && error("Failed to find a generic linear form after $max_iter tests.")
        end
        push!(cfs_lfs_out, coeffs)

        # Update every ideal's own running probe system with this form.
        for (idx, c) in enumerate(ctxs)
            cvars = gens(parent(first(c[:F])))
            L = sum(coeffs[i] * cvars[i] for i in eachindex(coeffs))
            push!(c[:F], L, vals[idx][1] * cvars[n - k + 1] + vals[idx][2])
        end
    end

    return ([c[:DEG] for c in ctxs], [c[:DIM] for c in ctxs]), cfs_lfs_out
end

# (reduction, degree, dimension) of `F` modulo `p`
function _deg_dim(F, p)
    Imod = Ideal(change_base_ring.(Ref(GF(p)), F))
    return Imod, hilbert_degree(Imod), dimension(Imod)
end

# Modular degree/dimension of `F` on *at most* three primes
function _reference_deg_dim(F)
    uZ = one(ZZRingElem)
    ps = _generate_lucky_primes(F, uZ<<30, (uZ<<31)-1, 3)
    I1, D1, M1 = _deg_dim(F, ps[1])
    I2, D2, M2 = _deg_dim(F, ps[2])
    (D1, M1) == (D2, M2) && return ps[1], I1, D1, M1
    _, D3, M3 = _deg_dim(F, ps[3])
    (D3, M3) == (D1, M1) && return ps[1], I1, D1, M1
    (D3, M3) == (D2, M2) && return ps[2], I2, D2, M2
    error("Degree/dimension differ modulo three distinct primes")
end

# Cheap (generic) evidence of non-radicality
# if the tests fails (e.g. unlucky fiber/prime)
# just returns true as an inconclusive test (will be catch by another test)
function _radical_evidence(Imod::Ideal{<:MPolyRingElem}, dim::Int, DEG::Int)
    dim < 0 && return true
    R = parent(Imod)
    n, vars, K = nvars(R), gens(R), base_ring(R)
    rnd() = K(rand(-100:100))
    cuts = [ sum(rnd() * v for v in vars) + rnd() for _ in 1:dim ]
    Isliced = Ideal(vcat(Imod.gens, cuts))

    dimension(Isliced) == 0 || return true
    El = eliminate(Isliced, n - 1)
    (isempty(El) || length(El) > 1) && return true
    Ielim = Ideal(El); Ielim.gb[0] = El
    hilbert_degree(Ielim) == DEG || return true
    # Actual radicality test (above detect degenerate situations)
    return is_squarefree(El[1])
end

# Checks, w.r.t. the degree reverse lexicographical order (msolve's default,
# mon_order = 0), that x_i^2 lies OUTSIDE the staircase for every variable x_i except `last_idx`.
# This is a LITERAL port of msolve's own pre-filter.
function _is_staircase_generic(I::Ideal{T} where T <: MPolyRingElem, last_idx::Int)
    n = nvars(parent(I))
    gb = groebner_basis(I, complete_reduction = false)

    # Degenerate case: `I` is already the whole ring
    any(g -> !iszero(g) && total_degree(g) == 0, gb) && return true

    lead_exps = [_lead_exp_ord(g, :degrevlex) for g in gb if !iszero(g)]
    for i in 1:n
        i == last_idx && continue
        e1 = zeros(Int, n); e1[i] = 1
        e2 = zeros(Int, n); e2[i] = 2
        any(e -> e == e1 || e == e2, lead_exps) || return false
    end
    return true
end

# Final validation of a linear form by computing a generic fiber with msolve
# **Note**: this evaluation is intentionally NOT reused later to avoid collision
# and unecessary technicalities
function _validate_against_msolve(real_F, n, DEG, val)
    R = parent(first(real_F))
    LFeval = _evalvar(real_F, n, [-val[2]//val[1]])[1]
    real_param = rational_parametrization(Ideal(LFeval))
    # Degenerate/empty fiber: original ideal is empty
    isempty(real_param.vars) && return true
    degree(real_param.elim) != DEG && return false

    return real_param.vars == symbols(R)[1:n-1]
end

# On a 0-dim ideal `Iext` of `DEG` points, embedded and pinned
# by the forms `cfs`: the first form (i.e. the last ring variable, the one
# `param_newvars` parametrizes by) must separate the points.
function _is_valid_points_form(Iext::Ideal, DEG::Int, cfs::Vector{Vector{ZZRingElem}}, prime)
    vars = gens(parent(Iext))
    n = length(vars)
    F = vcat(Iext.gens, transpose(cfs[1]) * vars, vars[n-1])
    Imod = Ideal(change_base_ring.(Ref(GF(prime)), F))

    _is_staircase_generic(Imod, n) || return false
    Imod_elim = Ideal(eliminate(Imod, n-1))
    Imod_elim.gb[0] = Imod_elim.gens
    return hilbert_degree(Imod_elim) == DEG
end

# The 3 linear forms used by the curve graph routines: found/validated on
# every curve and every 0-dim ideal of `points`, with no parametrization.
function _graph_linear_forms(curves::Vector{<:Ideal{T}},
                             points::Vector{<:Ideal{T}} = Ideal{T}[]) where T <: MPolyRingElem
    _, cfs_lfs = _add_genvars(collect(Ideal{T}, curves), 3, nothing; points=points)
    return cfs_lfs
end

# A generic evaluation point for the ideal `c`, avoiding 0 and multiples of
# its working prime
function _draw_val(c)
    val = [ZZ(), ZZ()]
    while iszero(val[1]) || is_divisible_by(val[1], c[:lucky_prime]) || is_divisible_by(val[2], c[:lucky_prime])
        val = rand(-c[:bif_bound]:c[:bif_bound], 2)
    end
    return val
end

# Check whether `coeffs` is a generic enough linear form (the k-th one
# to be fixed) for the ideal whose running probe system is `F`.
function _is_valid_linear_form(F, coeffs, val, DEG, DIM, k, F_orig, cfs_lfs_out, prime)
    R = parent(first(F))
    n, vars = nvars(R), gens(R)
    L = sum(coeffs[i] * vars[i] for i in 1:n)
    probe = val[1] * vars[n - k + 1] + val[2]
    FL = vcat(F, L)
    Feval = vcat(FL, probe)

    modK = GF(prime)
    Imod = Ideal(change_base_ring.(Ref(modK), 2*k <= DIM ? Feval : FL))

    if 2*k <= DIM
        # --- Case A: Projection Form ---
        return dimension(Imod) == DIM - 2*k && hilbert_degree(Imod) == DEG
    else
        # --- Case B: Generic staircase + Separating Form + msolve run (2*k == DIM + 1) ---
        # Intersects transversally the original variety
        if DIM > 0
            Imod0 = Ideal(change_base_ring.(Ref(modK), vcat(F_orig, L, probe)))
            dimension(Imod0) == DIM - 2 && hilbert_degree(Imod0) == DEG || return false
        end

        _is_staircase_generic(Imod, n - k + 1) || return false

        Imod_new = Ideal(change_ringvar(Imod.gens, symbols(R)[vcat(1:n-k, n-k+2:n, n-k+1)]))

        Imod_elim = Ideal(eliminate(Imod_new, n-1))
        Imod_elim.gb[0] = Imod_elim.gens

        hilbert_degree(Imod_elim) == DEG || return false

        # Both mod p tests passed, we need a final check on the exact system in QQ
        real_F = vcat(F_orig, [transpose(c) * vars for c in vcat(cfs_lfs_out, [coeffs])])
        return _validate_against_msolve(real_F, n, DEG, val) == true
    end
end

# A stateful, lazy generator for candidate linear forms with n vars
function _candidate_stream(n::Int)
    Channel{Vector{ZZRingElem}}() do ch
        # 1. Quick coordinate projection check backward from the last variable
        for i in n:-1:1
            coeffs = zeros(ZZRingElem, n)
            coeffs[i] = 1
            put!(ch, coeffs)
        end

        # 2. Bitmask Layering Search
        sorted_masks = sort(collect(1:(1 << (n-1)) - 1), by=count_ones)
        queue = [ones(ZZRingElem, n)]
        tested = Set{Vector{ZZRingElem}}([queue[1]])

        while true
            # Yield the current layer's candidates one by one
            for coeffs in queue
                put!(ch, coeffs)
            end

            # Generate the next layer
            new_L = [
                l .+ [ZZ((mask >> (k-1)) & 1) for k in 1:n]
                for mask in sorted_masks
                for l in queue
            ]

            queue = Vector{Vector{ZZRingElem}}()
            for coeffs in new_L
                if !(coeffs in tested)
                    push!(tested, coeffs)
                    push!(queue, coeffs)
                end
            end
        end
    end
end

function _evalvar(
    F::Vector{<:MPolyRingElem},
    i::Int,
    La::Vector{<:RingElem}
    )
    R = parent(first(F))
    indnewvars = setdiff(1:nvars(R), i)
    C, = polynomial_ring(base_ring(R), symbols(R)[indnewvars])

    LFeval = Vector{Vector{elem_type(C)}}()
    ctx = MPolyBuildCtx(C)

    max_deg = isempty(F) ? 0 : maximum(f -> degree(f, i), F)

    for a in La
        pow_a = [one(parent(a))]
        for d in 1:max_deg
            push!(pow_a, pow_a[end] * a)
        end

        push!(LFeval, elem_type(C)[])
        for f in F
            for (e, c) in zip(exponent_vectors(f), coefficients(f))
                aei = pow_a[e[i] + 1]
                push_term!(ctx, c * aei, [e[j] for j in indnewvars])
            end
            push!(LFeval[end], finish(ctx))
        end
    end
    return LFeval
end

# Generate N random primes between low and up
# that do not divide any numerator/denominator
# of any coefficient in polynomials from LP
function _generate_lucky_primes(
    LF::Vector{<:MPolyRingElem},
    low::ZZRingElem,
    up::ZZRingElem,
    N::Int64
    )
    # Using a Set avoids resizing and `unique!` shifting overhead
    CF_set = Set{ZZRingElem}()
    for f in LF, c in coefficients(f)
        !isone(numerator(c)) && push!(CF_set, numerator(c))
        !isone(denominator(c)) && push!(CF_set, denominator(c))
    end

    CF = sort!(collect(CF_set), rev=true)
    Lprim = ZZRingElem[]

    while length(Lprim) < N
        cur_prim = next_prime(rand(low:up))
        is_lucky = !(cur_prim in Lprim)
        idx = firstindex(CF)
        # Exploit decreasing order of CF
        while is_lucky && idx <= lastindex(CF) && CF[idx] > cur_prim
            is_lucky = !is_divisible_by(CF[idx], cur_prim)
            idx += 1
        end
        is_lucky && push!(Lprim, cur_prim)
    end
    return Lprim
end