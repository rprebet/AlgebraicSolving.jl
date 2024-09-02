## TODO : deal with isolated vertices (not appearing in the edges)
############################

function index_of(x, L)
   return findfirst(t->t==x,L)
end

function connected_components(G; Vemph=[], ind_abs=false)
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

    if length(Vemph)>0
        CVemph = [ [] for _ in eachindex(components) ]
        for i in eachindex(components)
            for j in eachindex(components[i])
                if components[i][j] in Vemph
                    if ind_abs
                        push!(CVemph[i], components[i][j])
                    else
                        push!(CVemph[i], j)
                    end
                end
            end
        end
    end

    if length(Vemph)>0
        return [ [[[ vert[cv] for cv in components[i] ], grouped_edges[i]], CVemph[i]] for i in eachindex(components) ]
    else
        return [ [[ vert[cv] for cv in components[i] ], grouped_edges[i]] for i in eachindex(components) ]
    end
end

function group_by_component(G, V)
    CG = connected_components(G, Vemph=V, ind_abs=true)
    return filter(c->length(c)>0, [ C[2] for C in CG ])
end

function number_connected_components(G)
    return length(connected_components(G))
end


function merge_graphs(LG, VC)
    # Initialize merged graph with the vertices and edges of the first graph
    V_merged, E_merged = deepcopy(LG[1])

    # To keep track of the mapping of vertices from each graph to the merged graph
    all_vertex_maps = [Dict(i => i for i in 1:length(V_merged))]

    # Iterate over each graph G_i in LG starting from the second graph
    for i in 2:length(LG)
        V_i, E_i = LG[i]  # Vertices and edges of the current graph
        VC_i = VC[i]      # Dictionary of common vertices for the current graph

        # Dictionary to map vertices of G_i to their new indices in V_merged
        vertex_index_map = Dict()

        # Step 1: Map common vertices from G_i to existing vertices in V_merged
        for (k, common_indices) in VC_i
            if k < i
                # Map each common vertex in G_i to its corresponding vertex in G_k (already merged)
                for (j, idx) in enumerate(common_indices)
                    vertex_index_map[idx] = all_vertex_maps[k][VC[k][i][j]]
                end
            end
        end

        # Step 2: Add non-common vertices from G_i to V_merged
        for idx in 1:length(V_i)
            if !haskey(vertex_index_map, idx)  # If not already mapped, it's a new vertex
                vertex_index_map[idx] = push!(V_merged, V_i[idx]) |> length
            end
        end

        # Step 3: Merge edges of G_i into the merged graph
        new_edges = [(vertex_index_map[i], vertex_index_map[j]) for (i, j) in E_i]
        append!(E_merged, new_edges)

        # Store the vertex index mapping for G_i for future merges
        push!(all_vertex_maps, vertex_index_map)
    end

    return V_merged, E_merged
end