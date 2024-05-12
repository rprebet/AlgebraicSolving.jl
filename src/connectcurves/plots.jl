
## Plot functions
function plot_graph(V, E)
    println("Plotting the graph")
    plot(legend=false)
    for e in E
        xe, ye = map(Float64,[e[1][1], e[2][1]]), map(Float64,[e[1][2], e[2][2]])
        plot!(xe, ye, lc="red", lw=3)
    end
    scatter!( map(Float64, [v[1] for v in V]),  map(Float64, [v[2] for v in V]), mc="black", m=:diamond)
    #plt.show(block=true)
    plot!()
end

function plot_graph_comp(V, CE)
    println("Plotting the graph")
    plot(legend=false)
    col = distinguishable_colors(length(CE)+2)
    for j in eachindex(CE)
        for e in CE[j]
            v1, v2 = map(Float64, V[e[1]]), map(Float64, V[e[2]])
            plot!([v1[1], v2[1]], [v1[2], v2[2]], lc=col[j+2], lw=3)
        end
    end
    #plt.show(block=true)
    plot!()
end