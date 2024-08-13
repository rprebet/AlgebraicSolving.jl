#using Pkg
#Pkg.activate("AlgebraicSolving.jl")
#using Revise
#using AlgebraicSolving
#using Nemo

export computeRM, computepolarproj

include("Cannytools.jl")

function computeRM(V::Ideal{T} where T <: QQMPolyRingElem, dimV::Int, Q=Vector{Vector{QQFieldElem}}([[]]) )
	# V is a variety of dimension dimV
	# Q are base points with rational coefficients
	# C are query points with rational coefficients

	A = parent(V)
	varias, nvarias = gens(A), nvars(A)
  hV = V.gens
  #println(Q)
  e = length(first(Q)) 
  fixvarias, newvarias = varias[1:e], varias[e+1:nvarias]

  R = Vector{AlgebraicSolving.Ideal{QQMPolyRingElem}}(undef,0)

  for q in Q
      ## Fq ##
      if e > 0
        hFq = change_ringvar([evaluate(h, fixvarias, q) for h in V.gens], collect(1:e))
        Fq = AlgebraicSolving.Ideal(hFq)
      else
        Fq = V
      end
      
      if dimV - e <= 1
        push!(R, AlgebraicSolving.Ideal(vcat(Fq.gens, [fixvarias[j] - q[j] for j in 1:e])))
      
      else
        ## sing(Fq) ##
        #println("Compute first the singular points")
        singFq = computepolarproj(0, Fq, dimV-e, newvarias, output="real", verb=0)
        @assert(length(singFq)==0, "Non-emtpy real sing locus!")
        #@assert(degree(singFq.elim)==0, "Non-emtpy sing locus!")

        ## K(pi_1,Fq) ##
        K1Fq = computepolarproj(1, Fq, dimV-e, newvarias, output="interval", verb=0)

        ## K(pi_2, Fq) ##
        K2Fq = computepolarproj(2, Fq, dimV-e, newvarias, output="minors")
        push!(R, AlgebraicSolving.Ideal(vcat(K2Fq.gens, [fixvarias[j] - q[j] for j in 1:e])))

        ## Points with vertical tg in K(pi_2, Fq) ##
        K1WmFq = computepolarproj(2, K2Fq, 1, newvarias, output="interval", dimproj=0, verb=0)

        ## New base points ##
        K1W = vcat(K1Fq, K1WmFq)
        K1WRat = MidRationalPoints(getindex.(K1W,1))
        newQ = [ vcat(q, [kv]) for kv in K1WRat ]

        RFq = computeRM(V, dimV, newQ)

        R = vcat(R, RFq)
      end
  end

  return R
end

#=
## Test ##
R,(x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])
#h = x1^2+x2^2+x3^2+x4^2-1
h = (x1^2+x2^2+x3^2+x4^2+9-1)^2-4*9*(x1^2+x2^2+x3^2)
h = evaluate(h,[x1,x2,x3],[x1+rand(-10:10)*x2+rand(-10:10)*x3+rand(-10:10)*x4,x2+rand(-10:10)*x3+rand(-10:10)*x4,x3+rand(-10:10)*x4])
V = AlgebraicSolving.Ideal([h])

carte = (computeRM(V, 3, [Vector{QQFieldElem}(undef,0)]))
#println(carte)
=#

