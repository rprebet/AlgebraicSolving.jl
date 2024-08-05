
## Plot functions
function plot_graph(G; color="red", width=3, vert = true)
    println("Plotting the graph")
    V, E = G
    plot(legend=false)
    for e in E
        v1, v2 = [ map(Float64, V[ee]) for ee in e ]
        plot!([v1[1], v2[1]], [v1[2], v2[2]], lc=color, lw=width)
    end
    if vert
        scatter!( map(Float64, [v[1] for v in V]),  map(Float64, [v[2] for v in V]), mc="black", m=:diamond)
    end
    gui()
end

function plot_graphs(CG; width=3, vert=true)
    println("Plotting the graph")
    plot(legend=false)
    col = distinguishable_colors(length(CG)+2)
    for j in eachindex(CG)
        V, E = CG[j]
        for e in E
            v1, v2 = [ map(Float64, V[ee]) for ee in e ]
            plot!([v1[1], v2[1]], [v1[2], v2[2]], lc=col[j+2], lw=width)
        end
        if vert
            scatter!( map(Float64, [v[1] for v in V]),  map(Float64, [v[2] for v in V]), mc="black", m=:diamond)
        end
    end
    gui()
end

function plot_graph_comp(G; width=3, vert=true)
    println("Plotting the graph")
    plot(legend=false)
    CG = connected_components(G)
    col = distinguishable_colors(length(CG)+2)
    for j in eachindex(CG)
        V, E = CG[j]
        for e in E
            v1, v2 = [ map(Float64, V[ee]) for ee in e ]
            plot!([v1[1], v2[1]], [v1[2], v2[2]], lc=col[j+2], lw=width)
        end
        if vert
            scatter!( map(Float64, [v[1] for v in V]),  map(Float64, [v[2] for v in V]), mc="black", m=:diamond)
        end
    end
    gui()
end