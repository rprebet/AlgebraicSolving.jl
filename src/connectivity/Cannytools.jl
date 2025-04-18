# Compute the set of points x where phi_j(T_x(V)) has dimension < dimproj
function computepolar(
        j::Int,                 # j-th first coordinate images of phi
        V::Ideal{P};            # input ideal
        phi::Vector{P} = P[],   # polynomial map in consideration (completed by sufficiently many projections)
        dimproj = j-1,          # maximum dimension of tangent space of phi
        v=0,                    # verbosity level
    ) where (P <: MPolyRingElem)
    V.dim == -1 && dimension(V)
    R = parent(V)
    n = nvars(R)
    c = n - V.dim
    nphi = length(phi)

    JW = transpose([ derivative(f, k) for k=max(j+1-nphi,0):n, f in vcat(V.gens, phi[1:min(j,nphi)])])
    sizeminors = c + min(nphi,j) + min(dimproj,j-1) - (j-1)
    minors = compute_minors(sizeminors, JW, R)

    return Ideal(vcat(V.gens, minors))
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

function detmpoly(A::Matrix{T} where T<:MPolyRingElem, R)
    # Get the size of the matrix
    n = size(A, 1)
    if n != size(A, 2)
        throw(ArgumentError("Matrix must be square"))
    end

    if n == 1
        return A[1, 1]
    end

    # Initialize the determinant polynomial
    detA = zero(R)

    # Compute the determinant polynomial
    for j = 1:n
        submatrix = A[2:end, [i for i = 1:n if i != j]]
        detA += (-1)^(1+j)*A[1, j] * detmpoly(submatrix, R)
    end

    return detA
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
    # It orders the [l_i,r_i], removes repetitions and computes the simplest
    # rational points between these intervals.
    if isempty(S)
        return QQFieldElem[]
    end
    S1 = sort(S, lt=(x, y) -> x[2] <= y[1])
    n = unique!(S1) |> length
    ratioP = Vector{QQFieldElem}(undef, n-1)
    for i in 1:(n- 1)
        ri, li1 = S1[i][2], S1[i+1][1]
        @assert ri < li1 "The intervals are not disjoint."
        #ratioP[i] = small_mid_point(ri, li1)
        ratioP[i] = simplest_between(ri, li1)
    end
    return ratioP
end

