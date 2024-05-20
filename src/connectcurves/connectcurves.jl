#println("\nLoading libraries and data..")
#using Nemo
#using Plots, Colors
#pythonplot()

export compute_graph, connected_components, number_connected_components,
 plot_graph, plot_graphs, plot_graph_comp, compute_param

include("tools.jl")
include("subresultants.jl")
include("isolate.jl")
include("boxes.jl")
include("graph.jl")
include("plots.jl")
include("arbtools.jl")
include("param-curve.jl")

function compute_graph(F; param=false, generic=false, precx=150, v=0, arb=false)

    if !(param)
        println("Compute rational parametrization...")
        @time begin
            Fparam = compute_param(F)
        end
    end

    return compute_graph_param(Fparam[1], generic=generic, precx=precx, v=v, arb=arb)
end

function compute_graph_param(f; generic=false, precx = 150,v=0, arb=false)
    println("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    println("!! Careful: this is a WIP version !!")
    println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

@time begin
    R = parent(f)
    x, y = gens(R)
    # Generic change of variables
    changemat = [1 0; 0 1]
    if  !generic
        changemat = map(QQ,rand(-500:500, 2, 2))
    end
    f = evaluate(f, collect(changemat*[x; y]));
    

    println("\nCompute parametrization of critical pts...")
    @time begin
    sr = subresultants(f, derivative(f,y), 2, list=true);
    # Take sqfree factors of the resultant
    sqr = collect(factor_squarefree(sr[1][1]))
    # Keep only deg>0 factors 
    filter(t->t[2]>0, sqr)
    # Order by multiplicity
    sqr = sort(sqr, by=t->t[2])
    # Group by multiplicity
    group_sqr = [ [R(1),i] for i in 1:sqr[end][2] ]
    for r in sqr
        group_sqr[r[2]][1] *= r[1]
    end
    # Construct the parametrization of the critical points
    params = [ [ q[1], -sr[2][1], sr[2][2] ] for q in group_sqr ];
    end

    if arb
        # TODO : check that no overlap between different isolations
        println("\nIsolating critical values with precision ", precx,"..")
        @time begin
        xcrit = [ isolate(first(p), prec=precx) for p in params ]
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
        LBcrit = [ sort(inter_solutions(AlgebraicSolving.Ideal([p[1],  p[3]*y-p[2]]), precision=precx),by=t->t[1]) for p in params ]
        xcrit = [ [ B[1] for B in lbcrit ] for lbcrit in LBcrit ]
        xcritpermut = order_permut2d(xcrit);
        end
    end

    print("\nTest for identifying singular boxes");ts=time();
    ########################################################
    ### TODO ###############################################
    ########################################################

    # Could be improved by handling nodes as extreme boxes:
    # when npcside = [2,2,0,0] just take nearest below and above
    # intersections b with the curves on the vertical sides
    # and change into npcside = [0,0,2,2]
    ## TODO : Refine only the intervals that need to be refined
    println("\nCompute intersections with critical boxes..")
    @time begin
    LPCside = Array{Any}(undef,length(LBcrit))
    ndig = maximum([Int(floor(log10(length(LB)))) for LB in LBcrit])
    for i in eachindex(LBcrit)
        ndigi = Int(floor(log10(length(LBcrit[i]))))
        LPCside[i] = Array{Any}(undef, length(LBcrit[i]))
        precxtmp = precx 
        while true
            flag = false
            for j in eachindex(LBcrit[i])
                print("mult=$i ; $(j)/$(length(LBcrit[i]))$(repeat(" ", ndig-ndigi+1))pts","\r")
                pcside = intersect_box(f, LBcrit[i][j], prec=precxtmp)
                npcside = [length(n) for (I, n) in pcside]
                if i == 1 && sum(npcside) > 2
                    precxtmp *= 2
                    println("\nRefine extreme boxes along x-axis to precision ", precxtmp)
                    refine_xboxes(params[1][1], LBcrit[1], precxtmp)
                    flag = true
                    break
                elseif i > 1 && sum(npcside[1:2]) != 0
                    precxtmp *= 2
                    println("\nRefine singular boxes along x-axis to precision ", precxtmp)
                    refine_xboxes(params[2][1], LBcrit[2], precxtmp)
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
                yinf = maximum([i for (i, yy) in pairs(I[s]) if yy[1] < ycrit[2]])
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
                ymax = minimum([i for (i, yy) in pairs(I[s]) if yy[2] > ycrit[1]])
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

    println("Graph computation")
    # Would be nice to have only one intermediate fiber (take the average of abscissa and ordinates) for plot
    # And even remove this fiber for the graph
    Vert = []
    Edg = []
    Corr = [[[[], [[], [], []], []] for j in xcrit[i] ] for i in eachindex(xcrit) ]
    Viso = []

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
        #println(vcat(nI[1], nI[2], [length(I[1])+1]))
        for k in 1:ymincrit-1
            push!(Vert, [xcmid, sum(I[1][k] + I[2][k])// 4])
            push!(Corr[i][j][2][1], length(Vert))
            push!(Edg, [Corr[i][j][1][k], length(Vert)])  # left
            push!(Edg, [length(Vert), Corr[i][j][3][k]])  # right
        end
        ###########################
        # The critical point
        ##########################
        # if it is an isolated point
        if isempty(nI[1]) && isempty(nI[2])
            #pass
            # We can add the isolated  vertex
            # push!(Vert, Pcrit[i][j])
            # push!(Corr[i][j][2][1], length(Vert))
            # We will subsequently add the vertex in the graph
            # push!(Viso, length(Vert))
        ############################################
        ## TO BE REPLACED BY APPSING IDENTIFICATOR ##
        ## works for space curves without nodes   ##
        ############################################
        # If we are dealing with a node
        elseif i == 2
            # We connect the pairwise opposite branches nI[1][1][i] and nI[1][2][i+1 mod 2], i=1,2
            push!(Edg, [Corr[i][j][1][nI[1][1]], Corr[i][j][3][nI[2][2]]])
            push!(Edg, [Corr[i][j][1][nI[1][2]], Corr[i][j][3][nI[2][1]]])
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

    #EdgPlot = [[Vert[k] for k in [i, j]] for (i, j) in Edg]
    #plot_graph(Vert, EdgPlot)
    #gui()

    #plot_graph_comp(Vert,CEdg)
    # Operate inverse change of variable if necessary
    if !(generic)
        Vert = [ changemat*v for v in Vert ]
    end
println("Total time:")
end
    return Vert, Edg
end