function Arb_to_rat(x)
	r = radius(x)
	return map(simplest_rational_inside, [x-r, x+r])
end

function rat_to_Arb(x, prec)
    x1,x2 = x
    xm, xd = ArbField(prec)((x1+x2)/2), ArbField(prec)(x2-x1)
    return ball(xm,xd)
end

function evaluate_Arb(f, x::ArbFieldElem)
	if is_zero(f)
		return zero(parent(x))
	else
		cf = coefficients_of_univariate(f)
		return evalpoly(x, cf)
	end
end

function evaluate_Arb(f, g, x::ArbFieldElem)
	@assert !is_zero(g) "Denominator must be non-zero"
	if is_zero(f)
		return zero(parent(x))
	else
		cf, cg = coefficients_of_univariate.([f, g])
		return evalpoly(x, cf) / evalpoly(x, cg)
	end
end

function Arb_eval(f, B, prec)
	BA = [ rat_to_Arb(b, prec) for b in B ]
	return evaluate(change_base_ring(ArbField(prec), f), BA)
end