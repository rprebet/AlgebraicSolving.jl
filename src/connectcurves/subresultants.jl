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
    # takes a sparse representation of a poly P  
    # and outputs a sparse rpz of P seen as a univariate poly in the i-th variable
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
    S = []
    s = leading_coefficient(Q)^(degree(P)-degree(Q))
    A = Q
    B = pseudorem(P,-Q)
    ring = parent(P)
    while true
        d = degree(A)
        e = degree(B)
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
                        a ÷= 2
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
        B = pseudorem(A,-B) / (s^delta * leading_coefficient(A))
        A = copy(C)
        s = leading_coefficient(A)
    end
end

#Bivariate subresultants
function subresultants(P::MPolyRingElem{T}, Q::MPolyRingElem{T}, idx) where T <: RingElement
    LPQ = map(poly_to_array, [P,Q])
    ULPQ = [ parray_asvar(lpq, idx) for lpq in LPQ ]

    R, x = polynomial_ring(base_ring(parent(P)), "x")
    S, y = polynomial_ring(R, "y")

    UP, UQ = [ S([ array_to_poly(l, R) for l in lpq ]) for lpq in ULPQ ]

    sr = subresultants(UP, UQ)

    # Get it back to initial polynomial ring
    Lsr = [ [collect(coefficients(csr)) for csr in coefficients(sr)] for sr in sr]
    newsr = []
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

    return newsr
end