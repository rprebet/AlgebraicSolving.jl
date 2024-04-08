# TODO: is this needed?
function num_eqns(X::LocClosedSet)
    return length(X.seq)
end

# for displaying locally closed sets
function Base.show(io::IO, lc::LocClosedSet)
    string_rep_seq = variety_string_rep(lc.seq)
    for (i, ineqn_set) in enumerate(lc.ineqns)
        ineqn_string = variety_string_rep(ineqn_set,
                                          sep = "*", lpar = "(", rpar = ")")
        if isone(i)
            string_rep_seq *= "\\" * ineqn_string
        else
            string_rep_seq *= (" ∪ " * ineqn_string)
        end
    end
    print(io, string_rep_seq)
end

function variety_string_rep(F::Vector{<:MPolyRingElem};
                            sep = ", ", lpar = "", rpar = "")
    string_rep = "V("
    for (i, f) in enumerate(F)
        if isone(i)
            string_rep *= (lpar * "$f" * rpar)
        else
            string_rep *= (sep * lpar * "$f" * rpar)
        end
    end
    string_rep *= ")"
    return string_rep
end

function ring(X::LocClosedSet)
    return parent(first(X.seq))
end

function is_lci(X::LocClosedSet)
    neqns = num_eqns(X)
    return all(gb -> codim(gb) == neqns, X.gbs)
end

function is_empty_set(X::LocClosedSet)
    R = ring(X)
    if isempty(X.gbs)
        return true
    elseif all(gb -> one(R) in gb, X.gbs)
        return true
    elseif all(isempty, X.gbs)
        return true
    end
    return false
end

function add_inequation!(X::LocClosedSet, h::P;
                             known_zds=P[]::Vector{P}) where P
    R = ring(X)
    for (i, gb) in enumerate(X.gbs)
        X.gbs[i] = saturate(vcat(gb, known_zds), h)
        push!(X.ineqns[i], h)
    end
end

function add_inequation(X::LocClosedSet, h::MPolyRingElem)
    Y = deepcopy(X)
    if isone(h)
        return Y
    end
    add_inequation!(Y, h)
    return Y
end

# which equations are hull equations needs to be managed outside of this function
function split(X::LocClosedSet, g::MPolyRingElem)
    tim = @elapsed X_min_g = add_inequation(X, g)
    @info "initial saturation time $(tim)"

    if is_empty_set(X_min_g)
        @info "equation vanishes"
    end

    X_hull_g = deepcopy(X)

    col_gbs = X_min_g.gbs
    hull_gbs = Vector{typeof(g)}[]
    new_ineqns = Vector{typeof(g)}[]
    R = ring(X)
    for (i, (X_gb, col_gb)) in enumerate(zip(X.gbs, col_gbs))
        if one(R) in col_gb
            @info "equation vanishes on one GB"
            tim = @elapsed new_gb = saturate(vcat(X_gb, [g]), X.ineqns[i])
            @info "adding equation time $(tim)"
            push!(hull_gbs, new_gb)
            push!(new_ineqns, X.ineqns[i])
            continue
        end
        sort(col_gb, by = p -> total_degree(p))
        H_rand = filter(!iszero, my_normal_form(random_lin_combs(col_gb), X_gb))
        new_ineqns_and_gbs = remove(X_gb, H_rand, known_eqns = [g])
        append!(new_ineqns, [vcat(X.ineqns[i], [f]) for f in (x -> x[1]).(new_ineqns_and_gbs)])
        append!(hull_gbs, (x -> x[2]).(new_ineqns_and_gbs))
    end
    X_hull_g.gbs = hull_gbs
    X_hull_g.ineqns = new_ineqns

    return X_hull_g, X_min_g
end

function remove(gb::Vector{P},
                H::Vector{P};
                known_eqns::Vector{P}=P[]) where P

    res = Tuple{P, Vector{P}}[]
    isempty(H) && return res

    R = base_ring(first(gb))
    h = first(H)
    tim = @elapsed gb1 = saturate(vcat(gb, known_eqns), h)
    @info "remove time $(tim) for degree $(total_degree(h))"
    if one(R) in gb1
        @info "is empty"
        return remove(gb, H[2:end], known_eqns=known_eqns)
    end
    push!(res, (h, gb1))
    tim = @elapsed G = filter(!iszero, my_normal_form(random_lin_combs(gb1), gb))
    @info "normal form time $(tim)"
    for htil in H[2:end]
        tim = @elapsed Grem = filter(!iszero, my_normal_form(G .* htil, gb))
        @info "constructing new equations time $(tim)"
        new_gbs = remove(gb, Grem, known_eqns=known_eqns)
        append!(res, new_gbs)
    end
    return res
end

# ------------------------ #
# --- helper functions --- #
# ------------------------ #

function saturate(F::Vector{P}, nz::P) where {P <: MPolyRingElem}
    return saturate(F, [nz])
end

function saturate(F::Vector{P}, nzs::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), vcat(["t$i" for i in 1:length(nzs)], ["x$i" for i in 1:nvars(R)]),
                              ordering = :degrevlex)
    Fconv = [convert_poly_to_t_ring(f, S) for f in F]

    for (i, h) in enumerate(nzs)
        ti = vars[i]
        push!(Fconv, ti*convert_poly_to_t_ring(h, S)-1)
    end

    gb = groebner_basis(Ideal(Fconv), complete_reduction = true,
                        eliminate = length(nzs))

    # convert back to original ring
    res = Vector{P}(undef, length(gb))

    for (i, p) in enumerate(gb)
        res[i] = convert_to_orig_ring(p, R) 
    end
    return res
end

function quotient(F::Vector{P}, nz::P) where {P <: MPolyRingElem}
    return quotient(F, [nz])
end

function quotient(F::Vector{P}, nzs::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), vcat(["t$i" for i in 1:length(nzs)], ["x$i" for i in 1:nvars(R)]),
                              ordering = :degrevlex)
    Fconv = [convert_poly_to_t_ring(f, S) for f in F]

    for (i, h) in enumerate(nzs)
        ti = vars[i]
        Fconv .*= ti
        push!(Fconv, (ti-1)*convert_poly_to_t_ring(h, S))
    end

    gb = groebner_basis(Ideal(Fconv), complete_reduction = true,
                        eliminate = length(nzs))

    # convert back to original ring
    res = Vector{P}(undef, length(gb))

    for (i, p) in enumerate(gb)
        res[i] = divides(convert_to_orig_ring(p, R), prod(nzs))[2]
    end
    return res
end

# assumes H is sorted by degree
function random_lin_combs(H::Vector{P}) where {P <: MPolyRingElem}
    res = P[]
    chr = characteristic(base_ring(first(H)))
    curr_deg = total_degree(first(H))
    curr_pol = zero(parent(first(H)))
    for h in H
        if total_degree(h) == curr_deg
            curr_pol += rand(1:chr-1)*h
        else
            push!(res, curr_pol)
            curr_pol = h
            curr_deg = total_degree(h)
        end
    end
    push!(res, curr_pol)
    return res
end

function random_lin_comb(F::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    res = zero(R)
    chr = characteristic(R)
    for f in F
        res += rand(1:chr-1)*f
    end
    return res
end

function convert_poly_to_t_ring(f::P, S::MPolyRing) where {P <: MPolyRingElem}
    ctx = MPolyBuildCtx(S)
    R = parent(f)
    nts = nvars(S) - nvars(R)
    for (e, c) in zip(exponent_vectors(f), coefficients(f))
        enew = vcat(zeros(Int, nts), e)
        push_term!(ctx, c, enew)
    end
    return finish(ctx)
end

function convert_to_orig_ring(f::P, R::MPolyRing) where {P <: MPolyRingElem}
    ctx = MPolyBuildCtx(R)
    S = parent(f)
    nts = nvars(S) - nvars(R)
    for (e, c) in zip(exponent_vectors(f), coefficients(f))
        push_term!(ctx, c, e[nts+1:end])
    end
    return finish(ctx)
end

function max_ind_sets(gb::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(gb))
    res = [trues(ngens(R))]

    lms = (Nemo.leading_monomial).(gb)
    for lm in lms
        to_del = Int[]
        new_miss = BitVector[]
        for (i, mis) in enumerate(res)
            nz_exps_inds = findall(e -> !iszero(e),
                                   first(Nemo.exponent_vectors(lm)))
            ind_var_inds = findall(mis)
            if issubset(nz_exps_inds, ind_var_inds)
                for j in nz_exps_inds
                    new_mis = copy(mis)
                    new_mis[j] = false
                    push!(new_miss, new_mis)
                end
                push!(to_del, i)
            end
        end
        deleteat!(res, to_del)
        append!(res, new_miss)
        unique!(res)
    end

    max_length = maximum(mis -> length(findall(mis)), res)
    filter!(mis -> length(findall(mis)) == max_length, res)
    return res
end

function codim(gb::Vector{P}) where P
    miss = max_ind_sets(gb)
    cd = length(findall(b -> !b, first(miss)))
    return cd
end

# for debugging
function check_decomp(F::Vector{P}, Xs::Vector{<:LocClosedSet}) where {P <: MPolyRingElem}
    println("checking decomp")
    gb_ch = F
    for X in Xs
        for gb in X.gbs
            g = random_lin_comb(gb)
            gb_ch = saturate(gb_ch, g)
        end
    end
    R = parent(first(F))
    gb_ch = saturate(gb_ch, last(gens(R)))
    return one(parent(first(F))) in gb_ch
end
