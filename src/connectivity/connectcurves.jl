#println("\nLoading libraries and data..")
#using Nemo
#using Plots, Colors
#pythonplot()

export compute_graph, connected_components, number_connected_components, group_by_component,
 plot_graph, plot_graphs, plot_graph_comp, compute_param

 # DEBUG
 export subresultants, diff, diff_list

include("tools.jl")
include("subresultants.jl")
include("isolate.jl")
include("boxes.jl")
include("graph.jl")
include("plots.jl")
include("arbtools.jl")

function compute_graph(F, C=[]; param=true, generic=true, precx=150, v=0, arb=true, int_coeff=false)
    if !(param)
        println("Compute rational parametrization...")
        @time begin
            Fparam = compute_param(F)[3]
        end
    else
        Fparam = F[1]
    end

    return compute_graph_param(Fparam, C, generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff)
end

function compute_graph_param(f, C=[]; generic=true, precx = 150, v=0, arb=true, int_coeff=false)
    println("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    println("!! Careful: this is a WIP version !!")
    println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    R = parent(f)
    x, y = gens(R)

    if int_coeff
        CD = lcm(map(denominator, collect(coefficients(f))))
        f *= CD
    end
    #println(f)
    # Generic change of variables
    changemat = [1 0; 0 1]
    if  !generic
        changemat = trimat_rand(QQ, 2, range=-100:100)
    end
    f = evaluate(f, collect(changemat*[x; y]));
    

    println("\nCompute parametrization of critical pts...")
    @time begin
    sr = subresultants(f, derivative(f,y), 2, list=true)
    # Take sqfree factors of the resultant
    sqr = collect(factor_squarefree(sr[1][1]))
    # Keep only deg>0 factors 
    filter!(t->t[2]>0, sqr)
    # Order by multiplicity
    sort!(sqr, by=t->t[2])
    # Group by multiplicity
    group_sqr = [ [R(1),i] for i in 1:sqr[end][2] ]
    for r in sqr
        group_sqr[r[2]][1] *= r[1]
    end
    # Construct the parametrization of the critical points
    params = [ [ q[1], -sr[2][1], sr[2][2] ] for q in group_sqr ];
    append!(params, C)
    end

    if arb
        # TODO : check that no overlap between different isolations
        println("\nIsolating critical values with precision ", precx,"..")
        @time begin
        xcrit = [ isolate(first(p), prec=precx) for p in params ]
        for i in eachindex(xcrit)
            for j in eachindex(xcrit[i])
                if xcrit[i][j][1]==xcrit[i][j][2]
                    xcrit[i][j] = [xcrit[i][j][1]-1//ZZ(1)<<precx, xcrit[i][j][1]+1//ZZ(1)<<precx]
                end
            end
        end
        xcritpermut = order_permut2d(xcrit);
        end

        println("\nComputing isolating critical boxes using Arb with precision ",max(precx,150),"..")
        @time begin
        #RR = ArbField(max(precx,150))
        precArb =precx
        Pcrit = [ [ [xc, evaluate_Arb(params[i][2], xc[1], precArb)/evaluate_Arb(params[i][3],xc[1],precArb)] for xc in xcrit[i]] for i in eachindex(xcrit) ]
        LBcrit = [ [ [ map(QQ, pc[1]), map(QQ, Arb_to_rat(pc[2])) ]  for pc in pcrit] for pcrit in Pcrit ]
        end
    else
        println("\nCompute critical boxes with msolve with precision ", precx,"..")
        @time begin
        LBcrit = [ sort(real_solutions(AlgebraicSolving.Ideal([p[1],  p[3]*y-p[2]]), precision=precx, interval=true),by=t->t[1]) for p in params ]
        for i in eachindex(LBcrit)
            for j in eachindex(LBcrit[i])
                for k in eachindex(LBcrit[i][j])
                    if LBcrit[i][j][k][1]==LBcrit[i][j][k][2]
                        LBcrit[i][j][k] = [LBcrit[i][j][k][1]-1//ZZ(1)<<precx, LBcrit[i][j][k][1]+1//ZZ(1)<<precx]
                    end
                end
            end
        end
        xcrit = [ [ B[1] for B in lbcrit ] for lbcrit in LBcrit ]
        xcritpermut = order_permut2d(xcrit);
        end
    end

    @time begin
    println("\nTest for identifying singular boxes");ts=time();
    ########################################################
    # For each mult of sing pts (keys) give the corresponding factors of the resultant (values)
    ### DEBUG: for now, only nodes for singular points, hardcode below ###
    mult_to_factor = Dict(2=>2:length(group_sqr))
    #########################################################
    mults = keys(mult_to_factor)
    Lfyk = diff_list(f, 2, maximum(mults))
    for m in mults
        for k in mult_to_factor[m]
            while true
                flag = false
                for j in eachindex(LBcrit[k])
                    pcrit = [ rat_to_Arb(c, precx) for c in LBcrit[k][j] ]
                    if contains_zero(evaluate(Lfyk[m+1], pcrit))
                        println("Refine singular boxes of multiplicity ", k)
                        # TODO refine_boxes(params[k], LBcrit[k])
                        #flag = true
                        break
                    end
                end
                flag || break
            end
        end
    end
    end

    @time begin
    # Could be improved by handling nodes as extreme boxes:
    # when npcside = [2,2,0,0] just take nearest below and above
    # intersections b with the curves on the vertical sides
    # and change into npcside = [0,0,2,2]
    ## TODO : Refine only the intervals that need to be refined
    println("\nCompute intersections with critical boxes..")
    LPCside = Array{Any}(undef,length(LBcrit))
    ndig = maximum([Int(floor(log10(max(1,length(LB))))) for LB in LBcrit])
    for i in eachindex(LBcrit)
        ndigi = Int(floor(log10(max(1,length(LBcrit[i])))))
        LPCside[i] = Array{Any}(undef, length(LBcrit[i]))
        precxtmp = precx 
        while true
            flag = false
            for j in eachindex(LBcrit[i])
                print("mult=$i ; $(j)/$(length(LBcrit[i]))$(repeat(" ", ndig-ndigi+1))pts","\r")
                pcside = intersect_box(f, LBcrit[i][j], prec=precxtmp)
                npcside = [length(n) for (I, n) in pcside]
                if (i == 1 && sum(npcside) > 2) || (i > 1 && sum(npcside[1:2]) != 0)
                    precxtmp *= 2
                    println("\nRefine boxes along x-axis to precision ", precxtmp)
                    refine_xboxes(params[i][1], LBcrit[i], precxtmp)
                    flag = true
                    break
                end
                LPCside[i][j] = pcside
            end
            flag || break
        end
        println("")
    end
    LnPCside = [ [[length(indI) for (L, indI) in PB] for PB in lpcside] for lpcside in LPCside ] 

    # Update extreme boxes
    for j in eachindex(LBcrit[1])
        # If the curve does not intersect the box only on vertical sides
        if !(LnPCside[1][j][1:2] == [0, 0])
            PCside, nPCside = LPCside[1][j], LnPCside[1][j]
            I = [ l[1] for l in PCside[3:end] ]
            nI = [ l[2] for l in PCside[3:end] ]
            # Locate the orientation of the extreme point
            # s is the index on the side where there are more branches
            # s=1: left; s=2: right
            s = argmax([length(I[1]), length(I[2])])
            # Ordinate range of the extreme point
            ycrit = LBcrit[1][j][2]
            # If it intersects on the bottom side
            if nPCside[1] == 1
                # yinf: the intersection with the vertical side just below the extreme point
                yinf = maximum([i for (i, yy) in pairs(I[s]) if yy[1] < ycrit[1]])
                # We vertically enlarge the box until it intersects on the horizontal side
                push!(LPCside[1][j][s + 2][2], yinf)
                LPCside[1][j][1][2] = []
                # We update the intersection numbers
                LnPCside[1][j][s + 2] += 1
                LnPCside[1][j][1] = 0
            end
            # If it intersects on the top side
            if nPCside[2] == 1 
                # ymax: the intersection with the vertical side just above the extreme point
                ymax = minimum([i for (i, yy) in pairs(I[s]) if yy[2] > ycrit[2]])
                # We vertically enlarge the box until it intersects on the horizontal side
                push!(LPCside[1][j][s + 2][2], ymax)
                LPCside[1][j][2][2] = []
                # We update the intersection numbers
                LnPCside[1][j][s + 2] += 1
                LnPCside[1][j][2] = 0
            end
        end
    end
    end

    println("\nGraph computation")
    # Would be nice to have only one intermediate fiber (take the average of abscissa and ordinates) for plot
    # And even remove this fiber for the graph
    Vert = []
    Edg = []
    Corr = [[[[], [[], [], []], []] for j in xcrit[i] ] for i in eachindex(xcrit) ]
    Viso = []
    Vcon = [ [] for _ in 1:length(C) ]

    Lapp = [[],[]]

    for ind in 1:length(xcritpermut)
        i, j = xcritpermut[ind]

        if ind > 1
            i1, j1 = xcritpermut[ind - 1]
        end
        if ind < length(xcritpermut)
            i2, j2 = xcritpermut[ind + 1]
            I2L, nI2L = LPCside[i2][j2][3]
        end

        PCside, nPCside = LPCside[i][j], LnPCside[i][j]
        I = [ l[1] for l in PCside[3:end] ]
        nI = [ l[2] for l in PCside[3:end] ]

        xcmid = sum(LBcrit[i][j][1])//2
        ycmid = sum(LBcrit[i][j][2])//2

        ymincrit = minimum(vcat(nI[1], nI[2], [length(I[1])+1]))
        # Construct vertices
        ###########################
        # On the vertical left side
        if ind > 1
            for k in 1:length(I[1])
                push!(Corr[i][j][1], Corr[i1][j1][3][k])
            end
        else
            for k in 1:length(I[1])
                push!(Vert, [xcrit[i][j][1], sum(I[1][k])//2])
                push!(Corr[i][j][1], length(Vert))
            end
        end
        ###########################
        # On the vertical right side
        if ind < length(xcritpermut)
            for k in 1:length(I[2])
                push!(Vert, [(xcrit[i][j][2] + xcrit[i2][j2][1])//2, sum(I[2][k] + I2L[k])//4])
                push!(Corr[i][j][3], length(Vert))
            end
        else
            for k in 1:length(I[2])
                push!(Vert, [xcrit[i][j][2], sum(I[2][k])//2])
                push!(Corr[i][j][3], length(Vert))
            end
        end
        ###########################
        # Below the critical point
        #println()
        #println(map(length,I))
        #println(nI[1],", ", nI[2], ", ", [length(I[1])+1])
        for k in 1:ymincrit-1
            push!(Vert, [xcmid, sum(I[1][k] + I[2][k])// 4])
            push!(Corr[i][j][2][1], length(Vert))
            push!(Edg, [Corr[i][j][1][k], length(Vert)])  # left
            push!(Edg, [length(Vert), Corr[i][j][3][k]])  # right
        end
        ###########################
        # The critical point
        ##########################
        # If we are dealing with a node
        if length(group_sqr)>1 && i == 2
            # if it is an isolated point
            if isempty(nI[1]) && isempty(nI[2])
                push!(Lapp[1], (i,j))
                #pass
                # We can add the isolated  vertex
                # push!(Vert, [xcmid, ycmid])
                # push!(Corr[i][j][2][1], length(Vert))
                # We will subsequently add the vertex in the graph
                # push!(Viso, length(Vert))
            ############################################
            ## TO BE REPLACED BY APPSING IDENTIFICATOR ##
            ## works for space curves without nodes   ##
            ############################################
            else
                # We connect the pairwise opposite branches nI[1][1][i] and nI[1][2][i+1 mod 2], i=1,2
                push!(Edg, [Corr[i][j][1][nI[1][1]], Corr[i][j][3][nI[2][2]]])
                push!(Edg, [Corr[i][j][1][nI[1][2]], Corr[i][j][3][nI[2][1]]])
                push!(Lapp[2], (i,j))
            end
        else
            # We can add the vertex
            push!(Vert, [xcmid, ycmid])
            push!(Corr[i][j][2][2], length(Vert))
            # We connect to the vertical left side of the critical box
            for k in nI[1]
                push!(Edg, [Corr[i][j][1][k], length(Vert)])
            end
            # We connect to the vertical right side of the critical box
            for k in nI[2]
                push!(Edg, [length(Vert), Corr[i][j][3][k]])
            end
            if i > length(params)-length(C)
                # If this is a control point
                push!(Vcon[i + length(C) - length(params)], length(Vert))
            end
        end
        ###########################
        # Above the critical point
        for k=(length(I[1]) - length(nI[1]) - ymincrit+1):-1:1
            push!(Vert, [xcmid, sum(I[1][end - k + 1] + I[2][end - k + 1])//4])
            push!(Corr[i][j][2][3], length(Vert))
            push!(Edg, [Corr[i][j][1][end - k + 1], length(Vert)])  # left
            push!(Edg, [length(Vert), Corr[i][j][3][end - k + 1]])  # right
        end
    end

    if Lapp != [[],[]]
        println("Removed apparent singularities: $(length(Lapp[1])) isolated  $(length(Lapp[2])) nodes")
    end
    #EdgPlot = [[Vert[k] for k in [i, j]] for (i, j) in Edg]
    #plot_graph(Vert, EdgPlot)
    #gui()

    #plot_graph_comp(Vert,CEdg)
    # Operate inverse change of variable if necessary
    if !(generic)
        Vert = [ changemat*v for v in Vert ]
    end

    if length(C)>0
        return (Vert, Edg), Vcon
    else
        return Vert, Edg
    end
end

function group_by_component(F::Vector{P}, C::Vector{P}; param=true, generic=true, precx=150, v=0, arb=true, int_coeff=false) where {P <: MPolyRingElem}
    # Compute a graph homeomorphic to Z(F) and return the vertices identified by Z(C)
    G, Vcon = compute_graph(F, C, param=param, generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff)
    # Compute the partition of the vertices in Vcon according to the connected component in G
    sort!(Vcon)
    CVcon = group_by_component(G, Vcon)
    
    # Convert the partition into abscissa order 
    index_map = Dict((val, idx) for (idx, val) in enumerate(Vcon))
    return [map(v -> index_map[v], C) for C in CVcon]
end