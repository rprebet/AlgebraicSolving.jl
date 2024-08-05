function PolyToBRstring(f)
    fc, fe = map(collect, [coefficients(f), exponent_vectors(f)])
    s = string(length(fc), " 2\n")
    M = ["$(fc[i]) $(fe[i][2]) $(fe[i][1])" for i in eachindex(fc)]
    s *= join(M, '\n')
    return s
end

function BRstringToPoly(s)
    Ls = split(s, "\n")
    Ls = [split(ls, " ") for ls in Ls]
    Ls = [[parse(Int, l) for l in ls] for ls in Ls]

    F, i = [], 1
    while i <= length(Ls)
        nterm = Ls[i][1]
        j, f = 1, 0
        while j <= nterm
            c = Ls[i + j]
            f += c[1] * y^(c[2]) * x^(c[3])
            j += 1
        end
        push!(F, f)
        i += nterm + 1
    end

    return length(F) == 1 ? first(F) : F
end

function ToBresultant(f, g, finput="/tmp/in.bs")
    @assert parent(f) == parent(g) "The two polynomials must belong to the same polynomial ring."
    @assert characteristic(parent(f)) == 0 "Positive characteristic is not handled"

    s = PolyToBRstring(f) * "\n" * PolyToBRstring(g)

    ff = open(finput, "w")
    write(ff, s)
    close(ff)
end

function FromBresultant(output="p", foutput="/tmp/out.bs")
    #x, y = symbols("x y")
    #A, (x, y) = ZZ[], PolynomialRing(ZZ, ["x", "y"])
    ff = open(foutput, "r")
    s = read(ff, String)
    s = replace(s, "\n" => "", ":" => "", ";" => "")
    L = eval(Meta.parse(s))
    close(ff)
    return L
end

function Bresultant(f, g; fname1="/tmp/in.bs", fname2="/tmp/out.bs", bspath="./bresultant", v=0)
    ToBresultant(f, g, fname1)
    out = run(pipeline(`$bspath -f $fname1 -o $fname2`, stderr="/tmp/err.log.bs"))
    if v > 0
    	ff = open("/tmp/err.log.bs", "r")
    	println(read(ff, String))
    	close(ff)
    end
    #if v == 1
    #    println(split(stderr(out), '\n')[end - 1] |> x -> split(x, ": ")[end] * "s")
    #elseif v == 2
    #    println(stderr(out))
    #end
    return FromBresultant(fname2)
end

#using Nemo
#R, (x,y) = Nemo.polynomial_ring(Nemo.ZZ, [:x,:y])
#f = x^2+4*x+x*y-y^2+3*y-4
#g = derivative(f,y)
