"""
randommodel(species::Int64, connectance::Float64)

Return a network of randomly assembled interactions according to
the cascade model.

This is essentially a wrapper for the `erdos_renyi` function from
Graphs.jl and just packes it into a SpeciesInteraction network

#### References

Erdős, Paul, and Alfréd Rényi. 1959. “On Random Graphs I.” 
Publicationes Mathematicae. https://doi.org/10.5486/PMD.1959.6.3-4.12.

Graphs.jl TODO
"""
function randommodel(species::Int64, links::Int64)
    
    N = erdos_renyi(species, links)

    adj = Graphs.SimpleGraphs.adj(N)

    edges = zeros(Bool, (species, species))

    for i in eachindex(adj)
        for j in eachindex(adj[i])
            edges[i,adj[i][j]] = 1
        end
    end

    edges = Binary(edges)
    nodes =  Unipartite(edges)
    return SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
end 
