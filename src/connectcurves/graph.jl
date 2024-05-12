
## Graph functions
function connected_components(vert, edges)
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
                    push!(grouped_edges[idx], (node, neighbor))
                end
            end
        end
    end

    grouped_vert = [ [ vert[cv] for cv in c ] for c in components ]
    return grouped_vert, grouped_edges
end