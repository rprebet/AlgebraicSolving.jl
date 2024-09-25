
## Plot functions
function plot_graph(G, Vemph::Vector{Vector{T}} where T<:Any; color="red", width=3, vert=true, subplt=false)
    if !subplt 
        println("Plotting the graph")
        plot(legend=false)
    end
    V, E = G
    #col = distinguishable_colors(length(Vemph)+2)
    for e in E
        v1, v2 = [ map(Float64, V[ee]) for ee in e ]
        plot!([v1[1], v2[1]], [v1[2], v2[2]], lc=color, lw=width)
    end
    if vert
        scatter!( map(Float64, [v[1] for v in V]),  map(Float64, [v[2] for v in V]), mc="black", m=:diamond)
    end
    for i in 1:length(Vemph)
        scatter!( map(Float64, [V[j][1] for j in Vemph[i]]),  map(Float64, [V[j][2] for j in Vemph[i]]), mc=color, m=:diamond)
    end
    subplt || gui()
end

function plot_graph(G; color="red", width=3, vert=true, subplt=false)
    plot_graph(G, Vector{Vector{QQMPolyRingElem}}(); color=color, width=width, vert=vert, subplt=subplt)
end

function plot_graph(G, Vemph::Vector{T} where T<:Any; color="red", width=3, vert=true, subplt=false)
    plot_graph(G, [Vemph], color=color, width=width, vert=vert, subplt=subplt)
end

function plot_graph(G, Vemph::Dict{Int, Vector{T}} where T<:Any; color="red", width=3, vert=true, subplt=false)
    plot_graph(G, Vemph|>values|>collect, color=color, width=width, vert=vert, subplt=subplt)
end

function plot_graphs(CG; width=3, vert=true, subplt=false)
    if !subplt
        println("Plotting the graph")
        plot(legend=false)
    end
    col = distinguishable_colors(length(CG)+2)
    for j in eachindex(CG)
        G, CVemph = length(CG[j][1])==2 ? CG[j] : (CG[j], Vector{Int}())
        plot_graph(G, CVemph, color=col[j+2], vert=vert, subplt=true)
    end
    subplt || gui()
end

function plot_graph_comp(G, Vemph=[]; width=3, vert=true, subplt=false)
    println("Plotting the graph")
    plot(legend=false)
    CG = connected_components(G, Vemph)
    plot_graphs(CG, width=width, vert=vert, subplt=true)
    subplt || gui()#savefig("/tmp/test.html")
end