function poly_to_usolve_in(f, finput="/tmp/in.us")
    @assert characteristic(parent(f)) == 0 "Positive characteristic is not handled"

		fc = coefficients_of_univariate(f)
    s = string(length(fc)-1, "\n")
    M = ["$(ZZ(fc[i]))" for i in 1:length(fc)]
    s *= join(M, '\n')
   
    ff = open(finput, "w")
    write(ff, s)
    close(ff)
end

function usolve_out_to_inter(foutput="/tmp/out.us")
	Ls = usolve_out_to_tuple(foutput)
	return tuple_to_inter(Ls)
end

function tuple_to_inter(Ls)
	return [ [ls[1]//QQ(2)^(ls[2]), (ls[1]+1)//QQ(2)^(ls[2])] for ls in Ls ]
end
	
function usolve_out_to_tuple(foutput="/tmp/out.us")
	ff = open(foutput, "r")
	s = read(ff, String)
	Ls = split(s, "\n")[1:end-1]
	Ls = [ split(ls, " ") for ls in Ls[2:end] ]
	return [ [eval(Meta.parse("ZZ($l)")) for l in ls] for ls in Ls ]
end

function tuple_to_usolve_out(Ls, foutput="/tmp/out.us")
	s = string(length(Ls), "\n")
  M = ["$(Ls[i][1]) $(Ls[i][2]) $(Ls[i][3])" for i in 1:length(Ls)]
  s *= join(M, '\n')
	
	ff = open(finput, "w")
  write(ff, s)
  close(ff) 
end

function usolve(f; precision=16,
fname1="/tmp/in.us", fname2="/tmp/out.us", uspath="./usolve", v=0, output="")
    poly_to_usolve_in(f, fname1)
    out = run(pipeline(`$uspath -f $fname1 -o $fname2 -v $v -p $precision`, stderr="/tmp/err.log.us", stdout="/tmp/out.log.us"))
    #if v == 1
    #    println(split(stderr(out), '\n')[end - 1] |> x -> split(x, ": ")[end] * "s")
    #elseif v == 2
    #    println(stderr(out))
    #end
    Ls = usolve_out_to_tuple(fname2)
    Is = tuple_to_inter(Ls)
    if output == "inter"
	    return Is
	  else
	  	ff = open(fname2, "r")
    	return Ls, Is
    end
end
	
	
#using Nemo
#R, (x,y) = polynomial_ring(ZZ, [:x,:y])
#f = x^2+5*x-4
