using Graphs

include("internals.jl")

"""
randommodel(S::Int64, L::Int64)

    Return a network of randomly assembled interactions according to
    the Erdős-Rényi model.

    This is essentially a wrapper for the `erdos_renyi` function from
    Graphs.jl and just packes it into a SpeciesInteraction network

    #### References

    Erdős, Paul, and Alfréd Rényi. 1959. “On Random Graphs I.” 
    Publicationes Mathematicae. https://doi.org/10.5486/PMD.1959.6.3-4.12.

    Graphs.jl TODO
"""
function randommodel(S::Int64, L::Int64)

    N = erdos_renyi(S, L)

    adj = Graphs.SimpleGraphs.adj(N)

    edges = _get_matrix(adj)

    edges = Binary(edges)
    nodes = Unipartite(edges)
    return SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
end
