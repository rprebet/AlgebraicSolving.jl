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
    fixvarias = varias[1:e]

    RM = Vector{Ideal{QQMPolyRingElem}}(undef,0)
    for q in Q
        ## Fq ##
        Fq = e > 0 ? Ideal(change_ringvar([evaluate(h, fixvarias, q) for h in V.gens], A.S[e+1:end])) : V
        # Genericity assumption (can be checked)
        if checks
            @assert(dimension(Fq) == V.dim - e, "Non-generic polar variety")
        else
            Fq.dim = V.dim - e
        end

        # Terminal case (dim <=1)
        if Fq.dim <= 1
            curve = change_ringvar(Fq.gens, A.S)
            push!(RM, Ideal(vcat(curve, [fixvarias[j] - q[j] for j in 1:e])))
            continue
        end

        ## sing(Fq) ##
        if checks
            v>0 && println("Check real quasi-smoothness")
            singFq = computepolar(0, Fq, v=max(v-1,0))
            @assert(isempty(real_solutions(singFq, info_level=max(v-1,0), nr_thrds=Threads.nthreads())),
                    "Non-empty real sing locus!")
        end

        ## K(pi_1,Fq) ##
        v>0 && println("V-critical points")
        K1Fq = computepolar(1, Fq, v=max(v-1,0))
        K1Fq = real_solutions(K1Fq, info_level=max(v-1,0), nr_thrds=Threads.nthreads(), interval=true)

        ## K(pi_2, Fq) ##
        v>0 && println("Polar variety")
        K2Fq = computepolar(2, Fq, v=max(v-1,0))
        if checks
            @assert(isone(dimension(K2Fq)), "Non-generic polar variety")
        else
            K2Fq.dim = 1
        end
        polar = change_ringvar(K2Fq.gens, A.S)
        push!(RM, Ideal(vcat(polar, [fixvarias[j] - q[j] for j in 1:e])))

        ## Points with vertical tg in K(pi_2, Fq) ##
        v>0 && println("W-critical points with vertical tangent")
        K1WmFq = computepolar(2, K2Fq, dimproj=0, v=max(v-1,0))
        K1WmFq = real_solutions(K1WmFq, info_level=max(v-1,0), nr_thrds=Threads.nthreads(), interval=true)

        ## New base and query points ##
        Cq = isempty(q) ? C : [ c[2:end] for c in C if c[1] == q[e]]
        K1W = vcat(K1Fq, K1WmFq)
        # Heuristic to be proven (Reeb's th)
        #K1W = K1W[2:end-1]
        ##########
        K1WRat = MidRationalPoints(first.(K1W), unique(first.(Cq)))
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


