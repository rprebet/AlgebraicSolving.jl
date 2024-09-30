# Return the polynomials in F, but injected in the polynomial ring with newvarias_S as new variables
function change_ringvar(F::Vector{P}, newvarias_S::Vector{Symbol}) where {P <: MPolyRingElem}
    R = parent(first(F))
    # Locate variables of R in newvarias
    to_varias = Vector{Int}(undef,0)
    for v in newvarias_S
        ind = findfirst(x->x==v, R.S)
        push!(to_varias, typeof(ind)==Nothing ? length(R.S)+1 : ind)
    end

    ind_novarias = setdiff(eachindex(R.S), to_varias)
    newR, newvarias = polynomial_ring(base_ring(R), newvarias_S)

    res = typeof(first(F))[]
    ctx = MPolyBuildCtx(newR)

    for f in F
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            @assert(all([ e[i]==0 for i in ind_novarias ]), "Occurence of old variable.s found!")
            push!(e, 0)
            push_term!(ctx, c, [e[i] for i in to_varias ])
        end
        push!(res, finish(ctx))
    end

    return res
end

function change_ringvar(F::Vector{P}, newvarias_S::Vector{Symbol}) where {P <: PolyRingElem}
    R = parent(first(F))
    A, (x,) = polynomial_ring(base_ring(R), [R.S])
    return change_ringvar([ evaluate(f, x) for f in F ], newvarias_S)
end

function change_ringvar(f::Union{MPolyRingElem, PolyRingElem}, newvarias_S::Vector{Symbol})
    return first(change_ringvar([f], newvarias_S))
end

function change_ringvar(F::Vector{P}, ind_newvarias::Union{Vector{I}, UnitRange{I}}) where {P <: MPolyRingElem, I <: Int64}
    R = parent(first(F))
    return change_ringvar(F, [R.S[i] for i in ind_newvarias])
end

function change_ringvar(f::MPolyRingElem, ind_newvarias::Union{Vector{I}, UnitRange{I}}) where {I <: Int64}
    R = parent(f)
    return first(change_ringvar([f], [R.S[i] for i in ind_newvarias]))
end

# Return the polynomials in F, but injected in the polynomial ring with the variables occuring in F
function change_ringvar(F::Vector{P}) where {P <: MPolyRingElem}
    union_varias = Set{Symbol}()
    for f in F
        union!(union_varias, map(Symbol, vars(f)) )
    end
    return change_ringvar(F, collect(union_varias))
end

function change_ringvar(f::MPolyRingElem)
    return first(change_ringvar([f]))
end

# Return the polynomials in F, but injected in the polynomial ring with newvarias_S as new variables
function change_ringvar_mod(F::Vector{P}, newvarias_S::Vector{Symbol}, oldvarias_S::Vector{Symbol}) where {P <: MPolyRingElem}
    R = parent(first(F))
    # Locate variables of R in newvarias
    to_varias = Vector{Int}(undef,0)
    for v in newvarias_S
        ind = findfirst(x->x==v, oldvarias_S)
        push!(to_varias, typeof(ind)==Nothing ? length(oldvarias_S)+1 : ind)
    end

    ind_novarias = setdiff(eachindex(oldvarias_S), to_varias)
    newR, newvarias = polynomial_ring(base_ring(R), newvarias_S)

    res = typeof(first(F))[]
    ctx = MPolyBuildCtx(newR)

    for f in F
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            @assert(all([ e[i]==0 for i in ind_novarias ]), "Occurence of old variable.s found!")
            push!(e, 0)
            push_term!(ctx, c, [e[i] for i in to_varias ])
        end
        push!(res, finish(ctx))
    end

    return res
end

function change_ringvar_mod(F::Vector{P}, newvarias_S::Vector{Symbol}, oldvarias_S::Vector{Symbol}) where {P <: PolyRingElem}
    R = parent(first(F))
    A, (x,) = polynomial_ring(base_ring(R), oldvarias_S)
    return change_ringvar_mod([ evaluate(f, x) for f in F ], newvarias_S, oldvarias_S)
end

function change_ringvar_mod(f::Union{MPolyRingElem, PolyRingElem}, newvarias_S::Vector{Symbol}, oldvarias_S::Vector{Symbol})
    return first(change_ringvar_mod([f], newvarias_S, oldvarias_S))
end

#function change_ringvar_mod(f::Union{MPolyRingElem, PolyRingElem}, newvarias_S::Vector{Symbol})
#    return first(change_ringvar_mod([f], newvarias_S, oldvarias_S))
#end

function computepolarproj(j::Int, V::AlgebraicSolving.Ideal, dimV::Int, varbs; dimproj=j-1, characteristic=0, output="minors", verb=0, nr_thrds=1)
    # Compute the set of points x where pi_j(T_x(V)) has dimension < dimproj
    @assert output in ["minors", "groebner", "real", "parametric", "interval"] "Wrong output parameter"
    R = parent(V)
    n, hV = R.nvars, V.gens
    c = n - dimV

    JW = transpose([ derivative(f, k) for k=j+1:n, f in hV ])
    sizeminors = c + dimproj - (j-1)
    hW = vcat(hV, compute_minors(sizeminors, JW, R))
     output_functions = Dict(
        "minors" => x -> x,
        "groebner" => x -> groebner_basis(x, info_level=verb,nr_thrds=nr_thrds),
        "real" => x -> real_solutions(x, info_level=verb,nr_thrds=nr_thrds),
        "interval" => x -> real_solutions(x, interval=true, info_level=verb,nr_thrds=nr_thrds),
        "parametric" => x -> rational_parametrization(x, info_level=verb,nr_thrds=nr_thrds)
    )

    return output_functions[output](AlgebraicSolving.Ideal(hW))
end

function computepolarphi(j::Int, V::AlgebraicSolving.Ideal, phi::MPolyRingElem, dimV::Int, varbs; dimproj=j-1, characteristic=0, output="minors", verb=0, nr_thrds=1)
    # Compute the set of points x where pi_j(T_x(V)) has dimension < dimproj
    @assert output in ["minors", "groebner", "real", "parametric", "interval"] "Wrong output parameter"
    R = parent(V)
    n, hV = R.nvars, V.gens
    c = n - dimV

    if j==0
        JW = transpose([ derivative(f, k) for k=1:n, f in hV ])
        sizeminors = c
    else
        JW = transpose([ derivative(f, k) for k=j:n, f in [hV;[phi]] ])
        sizeminors = c + 1 + dimproj - (j-1)
    end
    hW = vcat(hV, compute_minors(sizeminors, JW, R))
     output_functions = Dict(
        "minors" => x -> x,
        "groebner" => x -> groebner_basis(x, info_level=verb,nr_thrds=nr_thrds),
        "real" => x -> real_solutions(x, info_level=verb,nr_thrds=nr_thrds),
        "interval" => x -> real_solutions(x, interval=true, info_level=verb,nr_thrds=nr_thrds),
        "parametric" => x -> rational_parametrization(x, info_level=verb,nr_thrds=nr_thrds)
    )

    return output_functions[output](AlgebraicSolving.Ideal(hW))
end

function compute_minors(p, A, R)
    #Computes the p-minors of a matrix A
    n, m = size(A)
    rowsmins = collect(combinations(1:n, p))
    colsmins = collect(combinations(1:m, p))
    mins = Vector{eltype(A)}(undef, length(rowsmins) * length(colsmins))
    k = 1
    for rowsmin in rowsmins
        for colsmin in colsmins
            mins[k] = detmpoly(A[rowsmin, colsmin], R)
            k += 1
        end
    end

    return mins
end

function combinations(a, n)
    # Helper function to recursively generate combinations
    function _combinations(a, n, start, chosen)
        if length(chosen) == n
            return [chosen]
        elseif start > length(a)
            return Vector{Int}([])
        else
            # Include the current element and recurse
            include_current = _combinations(a, n, start + 1, [chosen; a[start]])
            # Exclude the current element and recurse
            exclude_current = _combinations(a, n, start + 1, chosen)
            return vcat(include_current, exclude_current)
        end
    end
    return _combinations(a, n, 1, Vector{Int}([]))
end

function detmpoly(A::Matrix{QQMPolyRingElem}, R)
    # Get the size of the matrix
    n = size(A, 1)
    if n != size(A, 2)
        throw(ArgumentError("Matrix must be square"))
    end

    if n == 1
        return A[1, 1]
    end

    # Initialize the determinant polynomial
    detA = QQMPolyRingElem(R,0)

    # Compute the determinant polynomial
    for j = 1:n
        submatrix = A[2:end, [i for i = 1:n if i != j]]
        detA += (-1)^(1+j)*A[1, j] * detmpoly(submatrix, R)
    end

    return detA
end

function contfrac_convergents(x::Rational{Int})
    q = Rational{Int}[]
    fpart, ipart = modf(x)
    push!(q, ipart)
    na, da, nb, db = Int(1), Int(0), Int(ipart), Int(1)
    while true
        fpart == 0 && return q
        x = inv(fpart)
        fpart, ipart = modf(x)
        na, da, nb, db = nb, db, na+ipart*nb, da+ipart*db
        push!(q, nb//db)
    end
end

function modfQQ(x::QQFieldElem)
    ipart, fpart = [ f(numerator(x), denominator(x)) for f in [div, rem]]
    return QQFieldElem(fpart, denominator(x)), ipart
end

function small_mid_point(a::QQFieldElem,b::QQFieldElem)
    a == b && return a
    if b < a
        a, b = b, a
    end
    x = (a+b)/ QQFieldElem(2)
    fpart, ipart = modfQQ(x)
    println("\n",fpart," , ", ipart)
    q1, q2 = QQFieldElem(1//0), QQFieldElem(ipart,1)
    while true
        #println("(",q1," ; ",q2,")")
        if numerator(fpart) == 0 || (q2 > a && q2 < b)
            return q2
        end
        x = inv(fpart)
        fpart, ipart = modfQQ(x)
        println(fpart," , ", ipart)
        q1, q2 = q2, QQFieldElem(numerator(q1) + ipart*numerator(q2), denominator(q1) + ipart*denominator(q2))
    end
end

function MidRationalPoints(S::Vector{Vector{QQFieldElem}})
    # S is a list of [ [l_1,r_1], ..., [l_n, r_n] ]
    # such that the [l_i, r_i] are rational and disjoint open intervals.
    #
    # It orders the [l_i,r_i] and computes the simplest
    # rational points between these intervals.
    n = length(S)
    if n==0
        return Vector{QQFieldElem}(undef, 0)
    end
    S1 = sort(S, lt=(x, y) -> x[2] <= y[1])
    ratioP = Vector{QQFieldElem}(undef, n-1)
    for i in 1:(n- 1)
        ri, li1 = S1[i][2], S1[i+1][1]
        @assert ri < li1 "The intervals are not disjoint."
        #ratioP[i] = small_mid_point(ri, li1)
        ratioP[i] = simplest_between(ri, li1)
    end
    return ratioP
end

