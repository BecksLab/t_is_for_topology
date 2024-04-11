"""
cascademodel(S::Int64, Co::Float64)

    Return a network of randomly assembled interactions according to
    the cascade model for a specified number of species and 
    connectance.

    Note this function has been ported from the EcologicalNetworks.jl
    source code and has been (minimally) modified for the purpose of 
    this project.


    #### References

    Cohen, J.E., Newman, C.M., 1985. A stochastic theory of community food webs I.
    Models and aggregated data. Proceedings of the Royal Society of London. Series
    B. Biological Sciences 224, 421–448. https://doi.org/10.1098/rspb.1985.0042

    Poisot, Timothée, Zacharie Belisle, Laura Hoebeke, Michiel Stock, and Piotr
    Szefer. 2019. “EcologicalNetworks.jl - Analysing Ecological Networks.” Ecography.
    https://doi.org/10.1111/ecog.04310.
"""
function cascademodel(S::Int64, Co::Float64)

    @assert 0.0 < Co < 0.5

    # Initiate matrix
    edges = zeros(Bool, (S, S))

    # For each species randomly asscribe rank e
    e = Random.sort(rand(S); rev=false)

    # Probability for linking two species
    p = 2*Co*S/(S - 1)

    for consumer in  1:S

        # Rank for a consumer
        rank = e[consumer]

        # Get a set of resources smaller than consumer
        potentialresources = findall(x -> x > rank, e)  # indices

        # Each gets a link with a probability p
        for resource in potentialresources

            # Take the resources and for each of them check a random number
            # if it is smaller than probability of linking two speceis create
            # a link in the A matrix

            # Random number for a potential resource
            rand() < p && (edges[consumer, resource] = true)
        end
    end

    edges = Binary(edges)
    nodes = Unipartite(edges)
    return SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
end
