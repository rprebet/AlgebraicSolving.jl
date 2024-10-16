
## Plot functions
function plot_graph(G, Vemph::Vector{Vector{T}} where T<:Any; color="red", width=3, vert=true, subplt=false)
    if !subplt
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
        plot(legend=false)
    end
    col = distinguishable_colors(length(CG)+2)
    for j in eachindex(CG)
        G, CVemph = length(CG[j][1])==2 ? CG[j] : (CG[j], Vector{Int}())
        plot_graph(G, CVemph, color=col[j+2], vert=vert, subplt=true)
    end
    #subplt || gui()
    savefig("/home/remi/Documents/gittravail/test.html")
    replace_width_height_in_file("/home/remi/Documents/gittravail/test.html", 1000, 750)
end

function plot_graph_comp(G, Vemph=[]; width=3, vert=true, subplt=false)
    plot(legend=false)
    CG = connected_components(G, Vemph)
    println(Vemph)
    plot_graphs(CG, width=width, vert=vert, subplt=true)
    subplt || gui()#savefig("/tmp/test.html")
end

function replace_width_height_in_file(filename, x, y)
    # Open the file and read its content as a string
    content = read(filename, String)

    # Build the patterns to search for and the replacements
    width_pattern = "\"width\": 600"
    height_pattern = "\"height\": 400"

    # Create the replacement strings with new values of x and y
    new_width = "\"width\": $x"
    new_height = "\"height\": $y"

    # Replace the last occurrence of "width": 600
    last_width_pos = findlast(occursin(width_pattern), split(content, "\n"))
    if last_width_pos !== nothing
        content = replace(content, width_pattern => new_width; count=1)
    end

    # Replace the last occurrence of "height": 400
    last_height_pos = findlast(occursin(height_pattern), split(content, "\n"))
    if last_height_pos !== nothing
        content = replace(content, height_pattern => new_height; count=1)
    end

    # Write the modified content back to the file
    write(filename, content)
end