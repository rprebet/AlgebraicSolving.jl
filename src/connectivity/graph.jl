## TODO : deal with isolated vertices (not appearing in the edges)
############################

function index_of(x, L)
   return findfirst(t->t==x,L)
end

function connected_components(G)
    # Outputs subgraphs of the connected components
    ###########
    vert, edges = G

    adj_list = Dict{Int, Vector{Int}}()
    visited = Set{Int}()
    components = Vector{Vector{Int}}()

    # Build adjacency list
    for (i, j) in edges
        push!(get!(adj_list, i, []), j)
        push!(get!(adj_list, j, []), i)
    end

    # Depth-first search
    function dfs(node, component)
        push!(visited, node)
        push!(component, node)
        for neighbor in adj_list[node]
            if !(neighbor in visited)
                dfs(neighbor, component)
            end
        end
    end

    # Find connected components
    for node in keys(adj_list)
        if !(node in visited)
            component = Vector{Int}()
            dfs(node, component)
            push!(components, component)
        end
    end

    # Group edges by connected component
    grouped_edges = [Tuple{Int, Int}[] for _ in components]
    for (idx, component) in enumerate(components)
        for node in component
            for neighbor in adj_list[node]
                if neighbor > node && neighbor in component
                    push!(grouped_edges[idx], (index_of(node,component), index_of(neighbor,component)))
                end
            end
        end
    end

    return [ [[ vert[cv] for cv in components[i] ], grouped_edges[i]] for i in eachindex(components) ]
end

function number_connected_components(G)
    return length(connected_components(G))
end