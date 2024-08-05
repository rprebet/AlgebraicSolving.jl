function in_inter(I, J)
    return J[1] <= I[1] && I[2] <= J[2]
end

function overlap_inter(I,J)
    return max(I[1], J[1]) <= min(I[2], J[2])
end

# To try/do : isolate with usolve, call msolve with only one variable
function intersect_box(f, B; prec=100)
    L = Array{Any}(undef, 4)
    for i in 1:2
        # Lxi
        L[i] = Array{Any}(undef,2)
        while true
            flag = false
            L[i][1] = isolate_eval(f, 2, B[2][i], prec=prec)
            L[i][2] = []
            for (j, l) in pairs(L[i][1])
                if in_inter(l, B[1])
                    push!(L[i][2], j)
                elseif overlap_inter(l, B[1])
                    prec *= 2
                    println("Increase precision to ", prec)
                    flag = true
                    break
                end
            end
            flag || break
        end
        # Lyi
        L[i+2] = Array{Any}(undef,2)
        while true
            flag = false
            L[i+2][1] = isolate_eval(f, 1, B[1][i], prec=prec)
            L[i+2][2] = []
            for (j, l) in pairs(L[i+2][1])
                if in_inter(l, B[2])
                    push!(L[i+2][2], j)
                elseif overlap_inter(l, B[2])
                    prec *= 2
                    #println("Increase precision to ", prec)
                    flag = true
                    break
                end
            end
            flag || break
        end
    end
    return L
end

function refine_xboxes(f, LB, prec)
    # Refine LB along first axis, being roots of f
    xnew = isolate(f, prec=prec)
    for i in eachindex(LB)
		  LB[i] = [ xnew[i], LB[i][2] ]
    end
end

function refine_ith_xboxes(f, LB, i, prec)
    # Refine only LB[i] along first axis, being i-th root of f
    xnew = isolate(f, prec=prec)
	LB[i] = [ xnew[i], LB[i][2] ]
end
