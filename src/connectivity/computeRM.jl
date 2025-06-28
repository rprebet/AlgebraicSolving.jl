#using Pkg
#Pkg.activate("AlgebraicSolving.jl")
#using Revise
#using AlgebraicSolving
#using Nemo

export roadmap, computepolar, MidRationalPoints, fbr, all_eqs, all_base_pts
include("Cannytools.jl")

function roadmap(
    V::Ideal{T} where T <: QQMPolyRingElem;                 # input ideal
    q::Vector{QQFieldElem}=QQFieldElem[],                   # single base point with rational coefficients
    C::Vector{Vector{QQFieldElem}}=Vector{QQFieldElem}[],   # query points with rational coefficients
    v::Int=0,                                               # verbosity level
    checks::Bool=false                                      # perform checks (dimension, regularity, etc.)
)
    # Some base cases
    A = parent(V)
    varias = gens(A)
    if length(varias)<=2
        return [V.gens]
    end
    # Some preprocessing
    isnothing(V.dim) && dimension(V)
    # Base points
    e = length(q)

    ## Fq ##
    # Genericity assumption (can be checked)
    if checks
        @assert(dimension(fbr(V,q)) == V.dim - e, "Non-generic coordinates")
    end

    # Terminal case (dim <=1)
    if V.dim - e <= 1
        return RMnode([], q, RMnode[])
    end

    ## sing(Fq) ##
    if checks
        v>0 && println("Check real quasi-smoothness")
        singFq = computepolar(e, V, v=max(v-1,0)) |> Ideal
        @assert(isempty(real_solutions(fbr(singFq, q), info_level=max(v-1,0), nr_thrds=Threads.nthreads())),
                "Non-empty real sing locus!")
    end

    ## K(pi_1,Fq) ##
    v>0 && println("V-critical points")
    K1Fq = computepolar(e+1, V, v=max(v-1,0)) |> Ideal
    K1Fq = real_solutions(fbr(K1Fq,q), info_level=max(v-1,0), nr_thrds=Threads.nthreads(), interval=true)

    ## K(pi_2, Fq) ##
    v>0 && println("Polar variety")
    K2Fqmins = computepolar(e+2, V, v=max(v-1,0), only_mins=true)
    K2Fq = vcat(V.gens, K2Fqmins) |> Ideal
    if checks
        @assert(fbr(K2Fq, q) |> dimension == 1, "Non-generic polar variety")
    else
        K2Fq.dim = e + 1
    end
    RM = RMnode(K2Fqmins, q, RMnode[])

    ## Points with vertical tg in K(pi_2, Fq) ##
    v>0 && println("W-critical points with vertical tangent")
    K1WmFq = computepolar(e+2, K2Fq, dimproj=e, v=max(v-1,0)) |> Ideal
    K1WmFq = real_solutions(fbr(K1WmFq,q), info_level=max(v-1,0), nr_thrds=Threads.nthreads(), interval=true)

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
        for newq in newQ
            RMFq = roadmap(V, q=newq, C=Cq)
            push!(RM.children, RMFq)
        end
    end

    if e == 0
        return Roadmap(V, RM)
    else
        return RM
    end
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


function fbr(F::Vector{P} where P <: QQMPolyRingElem, Q::Vector{QQFieldElem})
    @assert(!isempty(F), "Empty polynomial vector")
    vars = gens(parent(first(F)))
    return Ideal(vcat(F, [vars[i] - Q[i] for i in 1:min(length(vars),length(Q))]))
end

function fbr(I::Ideal{P} where P <: QQMPolyRingElem, Q::Vector{QQFieldElem})
    return fbr(I.gens, Q)
end

mutable struct RMnode
    polar_eqs::Vector{QQMPolyRingElem}
    base_pt::Vector{QQFieldElem}
    children::Vector{RMnode}
end

mutable struct Roadmap
    initial_ideal::Ideal{QQMPolyRingElem}
    root::RMnode
end

function all_eqs(RM::Roadmap)
    function all_eqs_rec(RMn::RMnode)
        eqs = [fbr(vcat(RM.initial_ideal.gens, RMn.polar_eqs), RMn.base_pt)]
        for child in RMn.children
            append!(eqs, all_eqs_rec(child))
        end
        return eqs
    end
    return all_eqs_rec(RM.root)
end

function all_base_pts(RM::Roadmap)
    function all_base_pts_rec(RMn::RMnode)
        eqs = [RMn.base_pt]
        for child in RMn.children
            append!(eqs, all_base_pts_rec(child))
        end
        return eqs
    end
    return all_base_pts_rec(RM.root)
end

function nb_nodes(RM::Roadmap)
    function nb_nodes_rec(RMn::RMnode)
        nb = 1
        for child in RMn.children
            nb += nb_nodes_rec(child)
        end
        return nb
    end
    return nb_nodes_rec(RM.root)
end

Base.show(io::IO, RM::Roadmap) = print(io, all_base_pts(RM))
Base.getindex(RM::Roadmap, idx::Union{Int, UnitRange}) = all_eqs(RM)[idx]
Base.lastindex(RM::Roadmap) = nb_nodes(RM)

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


