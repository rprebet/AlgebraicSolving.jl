function convert_to_msolve(
        F::Vector{T}) where T <: MPolyElem

    R = first(F).parent
    
    if R.ord != :degrevlex
        error("Only applicablefor DRL monomial ordering.")
    end

    nr_vars    = R.num_vars
    nr_gens    = length(F)
    lens       = Int32[F[i].length for i in 1:nr_gens]
    nr_terms   = sum(lens)
    field_char = characteristic(R)

    if field_char > 2^31
        error("At the moment we only support finite fields up to prime characteristic < 2^31.")
    end
    # get coefficients
    if field_char == 0
        cfs = BigInt[]
    else
        cfs = Int32[]
    end
    if field_char == 0
        for i in 1:nr_gens
            for cf in coefficients(F[i])
                push!(cfs, BigInt(numerator(cf)))
                push!(cfs, BigInt(denominator(cf)))
            end
        end
    else
        for i in 1:nr_gens
            for cf in coefficients(F[i])
                push!(cfs, Int32(data(cf)))
            end
        end
    end

    # get exponent vectors
    exps  = Int32[]
    for i in 1:nr_gens
        for ev in exponent_vectors(F[i])
            append!(exps, convert(Vector{Int32}, ev))
        end
    end

    return lens, cfs, exps
end

function convert_finite_field_gb_to_abstract_algebra(
        bld::Int32,
        blen::Vector{Int32},
        bcf::Vector{Int32},
        bexp::Vector{Int32},
        R::MPolyRing
        )

    if characteristic(R) == 0
        error("At the moment we only support finite fields.")
    end

    nr_gens = bld
    nr_vars = nvars(R)
    CR      = coefficient_ring(R)

    basis = Vector{MPolyElem}(undef, nr_gens)

    len   = 0

    for i in 1:nr_gens
        g  = MPolyBuildCtx(R)
        for j in 1:blen[i]
            push_term!(g, CR(bcf[len+j]),
                       convert(Vector{Int}, bexp[(len+j-1)*nr_vars+1:(len+j)*nr_vars]))
        end
        len +=  blen[i]
        basis[i]  = g.poly
    end

    return basis
end