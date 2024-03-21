println("Loading libraries and data..")
using Oscar
include("data.jl")
include("src/usolve/usolve.jl")
#include("tools.jl")

#using Arblib

function order_permut2d(L)
    # Create a list of tuples with elements and their corresponding indices
    LL = [(L[i][j], (i, j)) for i in eachindex(L) for j in eachindex(L[i])]
    # Sort the enumerated list based on the values
    sorted_LL = sort(LL, by = x -> x[1])
    # Extract the sorted values and their corresponding indices
    sorted_L = [pair[1] for pair in sorted_LL]
    sorted_ind = [pair[2] for pair in sorted_LL]
    
    return sorted_L, sorted_ind
end

function isolate(f; prec = 32, software="usolve")
	if software == "usolve"
		return usolve(f, precision = prec, uspath="src/usolve/usolve")
	end
end

function Arb_to_rat(x)
	r = radius(x)
	return map(simplest_rational_inside, [x-2*r, x+2*r])
end

function evaluate_Arb(f, x)
	cf = coefficients_of_univariate(f)
	return evalpoly(RR(x), cf) 
end

println("Isolating critical values..")
xcrit = [ isolate(first(p), prec=32) for p in params ]
xcrit_usolve = getindex.(xcrit, 1)
xcrit = getindex.(xcrit, 2)
xcritorder, xcritpermut = order_permut2d(xcrit)

println("Computing isolating critical boxes using Arb..")
#RR(x) = Arb(x, prec=120)
RR = ArbField(120)
Pcrit = [ [ [xc, evaluate_Arb(params[i][2], xc[1])/evaluate_Arb(params[i][3],xc[1])] for xc in xcrit[i]] for i in eachindex(xcrit) ]

LBcrit = [ [ [ pc[1], Arb_to_rat(pc[2]) ]  for pc in pcrit] for pcrit in Pcrit ]

printstep("Graph computation")
# Would be nice to have only one intermediate fiber (take the average of abscissa and ordinates) for plot
# And even remove this fiber for the graph
Vert = []
Edg = []
Corr = [[[[], [[], [], []], []] for _ in xcrit[i] ] for _ in xcrit ]
Viso = []
ite = 1

for ind in 1:length(xcritpermut)
    i, j = xcritpermut[ind]
    if ind > 1
        i1, j1 = xcritpermut[ind - 1]
    end
    if ind < length(xcritpermut)
        i2, j2 = xcritpermut[ind + 1]
        I2L, nI2L = LPCside[i2][j2][2]
    end

    PCside, nPCside = LPCside[i][j], LnPCside[i][j]
    I, nI = unzip(PCside[2:end])
    ymincrit = min(nI[1][1] + nI[1][2] + [length(I[1][1])])

    # Construct vertices
    ###########################
    # On the vertical left side
    if ind > 1
        for k in 1:length(I[1][1])
            push!(Corr[i][j][1][1], Corr[i1][j1][3][k])
        end
    else
        for k in 1:length(I[1][1])
            push!(Vert, (xcrit[i][j].lower() - exp, I[1][1][k]))
            push!(Corr[i][j][1][1], ite)
            ite += 1
        end
    end
    ###########################
    # On the vertical right side
    if ind < length(xcritpermut)
        for k in 1:length(I[1][2])
            push!(Vert, ((xcrit[i][j] + xcrit[i2][j2]).lower() / 2 + exp, (I[1][2][k] + I2L[k]).lower() / 2))
            push!(Corr[i][j][3], ite)
            ite += 1
        end
    else
        for k in 1:length(I[1][2])
            push!(Vert, (xcrit[i][j].upper() + exp, I[1][2][k]))
            push!(Corr[i][j][3], ite)
            ite += 1
        end
    end
    ###########################
    # Below the critical point
    for k in 1:ymincrit
        push!(Vert, (xcrit[i][j], (I[1][1][k] + I[1][2][k]) / 2))
        push!(Corr[i][j][2][1], ite)
        push!(Edg, (Corr[i][j][1][1][k], ite))  # left
        push!(Edg, (ite, Corr[i][j][3][k]))  # right
        ite += 1
    end
    ###########################
    # The critical point
    ##########################
    # if it is an isolated point
    if isempty(nI[1][1]) && isempty(nI[1][2])
        #pass
        # We can add the isolated  vertex
        # push!(Vert, Pcrit[i][j])
        # push!(Corr[i][j][2][1], ite)
        # We will subsequently add the vertex in the graph
        # push!(Viso, ite)
        # ite += 1
    ############################################
    ## TO BE REPLACED BY APPSING IDENTIFICATOR ##
    ## works for space curves without nodes   ##
    ############################################
    # If we are dealing with a node
    elseif i == 1
        # We connect the pairwise opposite branches nI[1][1][i] and nI[1][2][i+1 mod 2], i=1,2
        push!(Edg, (Corr[i][j][1][1][nI[1][1][1]], Corr[i][j][3][nI[1][2][2]]))
        push!(Edg, (Corr[i][j][1][1][nI[1][1][2]], Corr[i][j][3][nI[1][2][1]]))
    else
        # We can add the vertex
        push!(Vert, Pcrit[i][j])
        push!(Corr[i][j][2][1], ite)
        # We connect to the vertical left side of the critical box
        for k in nI[1][1]
            push!(Edg, (Corr[i][j][1][1][k], ite))
        end
        # We connect to the vertical right side of the critical box
        for k in nI[1][2]
            push!(Edg, (ite, Corr[i][j][3][k]))
        end
        ite += 1
    end
    ###########################
    # Above the critical point
    for k in length(I[1][1]) - length(nI[1][1]) - ymincrit:-1:1
        push!(Vert, (xcrit[i][j], (I[1][1][end - k + 1] + I[1][2][end - k + 1]) / 2))
        push!(Corr[i][j][2][3], ite)
        push!(Edg, (Corr[i][j][1][1][end - k + 1], ite))  # left
        push!(Edg, (ite, Corr[i][j][3][end - k + 1]))  # right
        ite += 1
    end
end