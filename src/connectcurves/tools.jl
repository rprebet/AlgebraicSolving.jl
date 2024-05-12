#using Nemo
using Singular

function randQQ(size=3)
	fac = 10^(size)
	return Nemo.QQ((-1)^(trunc(Int,rand()*10))*trunc(Int,rand()*fac),trunc(Int,rand()*fac))
end

function _mpoly_to_singular(f)
	@assert Nemo.is_univariate(f) "Not univariate polynomial"
	R, x = Singular.polynomial_ring(Singular.QQ, [:x])
	x = first(x)
	cf = Nemo.coefficients_of_univariate(f)
	return sum([ cf[j]*x^(j-1) for j = eachindex(cf) ])
end

function nb_roots_interval(f,a,b)
	g = _mpoly_to_singular(f)
	return Singular.LibRootsur.sturm(g, Singular.QQ(a), Singular.QQ(b))
end

"""
	nb_roots_interval_partition(f, a, b)

Compute the number of real roots of f(x) in the intervals:
	(-infinity, a), (a,b), (b, +infinity)

# Input:
* f: polynomial with involving one variable, with coeffients in QQ
* a,b:  rational elements of QQ
# Output:
* [ n(-infinity, a), n(a,b), n(b, +infinity) ]

"""
function nb_roots_interval_partition(f, a, b)
	g = _mpoly_to_singular(f)
	M = Singular.LibRootsur.maxabs(g)
	a1, b1 = Singular.QQ(a), Singular.QQ(b)
	M1, M2 = min(-M, a1), max(b1, M)
	B = [ M1-1, a1, b1, M2+1 ]
	return [ Singular.LibRootsur.sturm(g, B[i], B[i+1]) for i=1:3 ]
end

function sturm_habicht_seq(f)
	g = _mpoly_to_singular(f)
	SH = map(Singular.coefficients_of_univariate, Singular.LibRootsur.sturmhaseq(g))
	return [ vcat(SH[k], [n_Q(0) for _=1:(k-length(SH[k]))]) for k in eachindex(SH) ]
end

function subresultants_seq(f)
	SH = sturm_habicht_seq(f)
	d(k) = (-1)^((k*(k+1))>>1)
	return [ d(length(SH)-1-k)*SH[k] for k=1:length(SH)-2 ]
end

function subresultants_seq_biv(f, iv)
	iu = iv%2+1
	# Modern computer Algebra Th 6.51.
	degs = degrees(f)
	n, d = degs[iv], degs[iu]
	degsmax = [ (2*n-1-2*k)*d for k=0:d-2 ] 
	Dmax = degsmax[end]>>1+1
	##
	xs = [ randQQ() for _=1:(degsmax[1]+1) ]
	@time begin
	Srev = [ subresultants_seq(evaluate(f,[iu],[xx])) for xx in xs ]
	end
	S,t = Nemo.polynomial_ring(Nemo.QQ,"t")
	Sr = [ [ S(0) for _=1:i] for i=1:n-2 ]
	for i = 1:n-2
		for j = 1:i
			ys = [ Nemo.QQ(Srev[l][i][j]) for l=eachindex(Srev) ]
			Sr[i][j] = Nemo.interpolate(S, xs, ys)
		end
	end
	return Sr
end

#R, (x,y) = Nemo.polynomial_ring(Nemo.QQ, [:x,:y])
#f = x^7+y^2*x^5-2*x*y^3+1+y^6
#f = (x^2+y^2-1)*(2*x-y)
#Sr = subresultants_seq_biv(f, 1)
