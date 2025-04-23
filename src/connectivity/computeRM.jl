#using Pkg
#Pkg.activate("AlgebraicSolving.jl")
#using Revise
#using AlgebraicSolving
#using Nemo

export roadmap, computepolar, MidRationalPoints
include("Cannytools.jl")

function roadmap(
        V::Ideal{T} where T <: QQMPolyRingElem;                 # input ideal
        Q::Vector{Vector{QQFieldElem}}=Vector{QQFieldElem}[],   # base points with rational coefficients
        C::Vector{Vector{QQFieldElem}}=Vector{QQFieldElem}[],   # query points with rational coefficients
        v::Int=0,                                               # verbosity level
        checks::Bool=false                                      # perform checks (dimension, regularity, etc.)
    )
    # Some base cases
    A = parent(V)
    varias = gens(A)
    if length(varias)<=2
        return [V]
    end
    # Some preprocessing
    V.dim == -1 && dimension(V)
    isempty(Q) && push!(Q,[])
    # Base points
    e = length(first(Q))

    RM = Vector{Ideal{QQMPolyRingElem}}(undef,0)
    for q in Q
        ## Fq ##
        # Genericity assumption (can be checked)
        Fq = fbr(V, q)
        if checks
            @assert(dimension(Fq) == V.dim - e, "Non-generic polar variety")
        else
            Fq.dim = V.dim - e
        end

        # Terminal case (dim <=1)
        if Fq.dim <= 1
            push!(RM, Fq)
            continue
        end

        ## sing(Fq) ##
        if checks
            v>0 && println("Check real quasi-smoothness")
            singFq = fbr(computepolar(e, V, v=max(v-1,0)), q)
            @assert(isempty(real_solutions(singFq, info_level=max(v-1,0), nr_thrds=Threads.nthreads())),
                    "Non-empty real sing locus!")
        end

        ## K(pi_1,Fq) ##
        v>0 && println("V-critical points")
        K1Fq = fbr(computepolar(e+1, V, v=max(v-1,0)), q)
        K1Fq = real_solutions(K1Fq, info_level=max(v-1,0), nr_thrds=Threads.nthreads(), interval=true)

        ## K(pi_2, Fq) ##
        v>0 && println("Polar variety")
        K2Fqmins = computepolar(e+2, V, v=max(v-1,0))
        K2Fq = fbr(K2Fqmins, q)
        if checks
            @assert(isone(dimension(K2Fq)), "Non-generic polar variety")
        else
            K2Fq.dim = 1
        end
        push!(RM, K2Fq)

        ## Points with vertical tg in K(pi_2, Fq) ##
        v>0 && println("W-critical points with vertical tangent")
        K1WmFq = fbr(computepolar(e+2, K2Fqmins, dim=e+1, dimproj=e+0, v=max(v-1,0)), q)
        K1WmFq = real_solutions(K1WmFq, info_level=max(v-1,0), nr_thrds=Threads.nthreads(), interval=true)

        ## New base and query points ##
        Cq = isempty(q) ? C : [ c for c in C if c[e] == q[e]]
        K1W = vcat(K1Fq, K1WmFq)
        # Heuristic to be proven (Reeb's th)
        #K1W = K1W[2:end-1]
        ##########
        K1WRat = MidRationalPoints(getindex.(K1W,e+1), unique(getindex.(Cq, e+1)))
        newQ = vcat.(Ref(q), K1WRat)

        # Recursively compute roadmap of possible fibers
        if !isempty(newQ)
            RMFq = roadmap(V, Q=newQ, C=Cq)
            append!(RM, RMFq)
        end
    end

    return RM
end

function roadmap(
    V::Ideal{P},                # input ideal
    C::Ideal{P};                # ideal defining query points
    v::Int=0,                   # verbosity level
    checks::Bool=false          # perform checks (dimension, regularity, etc.)
) where (P <: QQMPolyRingElem)
    @assert(parent(V)==parent(C), "Equations for variety and query points must live the same ring")
    CQ = real_solutions(C, info_level=max(v-1,0), nr_thrds=Threads.nthreads())
    return roadmap(V, C=CQ, v=v, checks=checks)
end

function fbr(I::Ideal{P} where P <: QQMPolyRingElem, Q::Vector{QQFieldElem})
    vars = gens(parent(I))
    return Ideal(vcat(I.gens, [vars[i] - Q[i] for i in 1:min(length(vars),length(Q))]))
end

#=
## Test ##
R,(x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])
#h = x1^2+x2^2+x3^2+x4^2-1
h = (x1^2+x2^2+x3^2+x4^2+9-1)^2-4*9*(x1^2+x2^2+x3^2)
h = evaluate(h,[x1,x2,x3],[x1+rand(-10:10)*x2+rand(-10:10)*x3+rand(-10:10)*x4,x2+rand(-10:10)*x3+rand(-10:10)*x4,x3+rand(-10:10)*x4])
V = Ideal([h])

carte = (computeRM(V, 3, [Vector{QQFieldElem}(undef,0)]))
#println(carte)
=#


