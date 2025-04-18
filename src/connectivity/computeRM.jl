#using Pkg
#Pkg.activate("AlgebraicSolving.jl")
#using Revise
#using AlgebraicSolving
#using Nemo

export computeRM, computepolarproj, computepolarphi, computepolar
# DEBUG
export compute_minors, detmpoly, change_ringvar_mod, MPolyBuild

include("Cannytools.jl")

function computeRM(
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
    println(Q)
    println(C)
    println()
    isempty(Q) && push!(Q,[])
    # Base points
    e = length(first(Q))
    fixvarias = varias[1:e]

    RM = Vector{Ideal{QQMPolyRingElem}}(undef,0)
    for q in Q
        ## Fq ##
        if e > 0
            hFq = change_ringvar([evaluate(h, fixvarias, q) for h in V.gens], A.S[e+1:end])
            Fq = Ideal(hFq)
        else
            Fq = V
        end
        # Genericity assumption (can be checked)
        Fq.dim = V.dim - e

        if V.dim - e <= 1
            curve = change_ringvar(Fq.gens, A.S)
            push!(RM, Ideal(vcat(curve, [fixvarias[j] - q[j] for j in 1:e])))
        else
            ## sing(Fq) ##
            if checks
                v>0 && println("Compute first the singular points")
                singFq = computepolar(0, Fq, v=max(v-1,0))
                @assert(isempty(real_solutions(singFq, info_level=max(v-1,0), nr_thrds=Threads.nthreads())), "Non-emtpy real sing locus!")
            end

            ## K(pi_1,Fq) ##
            v>0 && println("First critical points")
            K1Fq = computepolar(1, Fq, v=max(v-1,0))
            K1Fq = real_solutions(K1Fq, info_level=max(v-1,0), nr_thrds=Threads.nthreads(), interval=true)

            ## K(pi_2, Fq) ##
            v>0 && println("Second critical points")
            K2Fq = computepolar(2, Fq, v=max(v-1,0))
            if checks
                @assert(isone(dimension(K2Fq)), "Non-generic polar variety")
            end
            K2Fq.dim = 1
            polar = change_ringvar(K2Fq.gens, A.S)
            push!(RM, Ideal(vcat(polar, [fixvarias[j] - q[j] for j in 1:e])))

            ## Points with vertical tg in K(pi_2, Fq) ##
            v>0 && println("Vertical tg points")
            K1WmFq = computepolar(2, K2Fq, dimproj=0, v=max(v-1,0))
            K1WmFq = real_solutions(K1WmFq, info_level=max(v-1,0), nr_thrds=Threads.nthreads(), interval=true)

            ## New base points ##
            K1W = vcat(K1Fq, K1WmFq)
            K1WRat = MidRationalPoints(getindex(K1W, 1))
            # Heuristic to be proven
            #K1WRat = K1WRat[2:end-1]
            ##########
            Cq = [c for c in C if c[1:e]==q]
            append!(K1WRat, getindex.(Cq, e+1)) |> sort!
            newQ = [ vcat(q, [kv]) for kv in K1WRat ]

            if !isempty(newQ)
                RMFq = computeRM(V, Q=newQ, C=Cq)
                append!(RM, RMFq)
            end
        end
    end

    return RM
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


