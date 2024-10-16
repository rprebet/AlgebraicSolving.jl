export parse_poly

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

function parse_poly(poly_str)
    # Regular expression to match terms of the form "ax^n", "bx", or constants
    term_pattern = r"([+-]?[ ]*\d*)[\*]?x\^(\d+)|([+-]?[ ]*\d*)[\*]?x|([+-]?[ ]*\d+)"

    # Initialize a dictionary to store the coefficients with their exponents
    coeffs = Dict{Int, BigInt}()

    # Use the regular expression to find all matches
    for m in eachmatch(term_pattern, poly_str)
        #println(m)
        if !isnothing(m.captures[1]) && !isnothing(m.captures[2])  # Matches ax^n
            s = m.captures[1] * (m.captures[1] in ["","-","+"] ? "1" : "")
            coeffs[parse(Int, m.captures[2])] = parse(BigInt, s)
        elseif !isnothing(m.captures[3])  # Matches bx (degree 1 term)
            s = m.captures[3] * (m.captures[3] in ["","-","+"] ? "1" : "")
            coeffs[1] = parse(BigInt, s)
        elseif !isnothing(m.captures[4])  # Matches constant term
            coeffs[0] = parse(BigInt, m.captures[4])
        end
    end

    # Get the highest degree
    max_degree = maximum(keys(coeffs))
    # Create the list of coefficients starting from the highest degree
    result = [get(coeffs, i, 0) for i in 0:max_degree]

    return result
end

function FromBresultant(output="p", foutput="/tmp/out.bs")
    #x, y = symbols("x y")
    A, x = polynomial_ring(QQ, :x)
    ff = open(foutput, "r+")
    s = read(ff, String)
    #println(s)
    Ls = split.(split(replace(s, r":|;\n|\[\n\[|\]\n\]" => ""), "],\n["), Ref(",\n"))
    Lcs = [ parse_poly.(ls) for ls in Ls ]
    #println(Lcs)
    return [ A.(ls) for ls in Lcs ]
end

function Bresultant(f, g; fname1="/tmp/in.bs", fname2="/tmp/out.bs", bspath="./bresultant", v=0)
    @time ToBresultant(f, g, fname1)
    #err = Pipe()
    @time process = run(pipeline(`$bspath -f $fname1 -o $fname2`))

    #=
    #close(err.in)
    if v > 0
        c=0
        while c <100
            out = @async String(read(err))
            println(out)
            c +=1
        end
    end
    #if v == 1
    #    println(split(stderr(out), '\n')[end - 1] |> x -> split(x, ": ")[end] * "s")
    #elseif v == 2
    #    println(stderr(out))
    #end
    wait(process)
    =#

    return @time FromBresultant(fname2)
end

#using Nemo
#R, (x,y) = Nemo.polynomial_ring(Nemo.ZZ, [:x,:y])
#f = x^2+4*x+x*y-y^2+3*y-4
#g = derivative(f,y)
