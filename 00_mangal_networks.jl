using Mangal
using SpeciesInteractionNetworks

network_count = 3
mangal_networks = Mangal.networks("type" => "herbivory", "count" => network_count)
# 20 datasets for now just for initial testing

# make a nice 'dataframe' to store network data
mutable struct MangalNetworks
    id::Vector{Int8}
    richness::Vector{Int8}
    links::Vector{Int8}
    connectance::Vector{Float64}
    # we can add some 'summary stats' here as well
end

function MangalNetworks(C::Int64)
    return MangalNetworks(
        zeros(Int8, C),
        zeros(Int8, C),
        zeros(Int8, C),
        zeros(Float64, C),
    )
end

mangal_topology = MangalNetworks(network_count); 

for i in eachindex(mangal_networks)
    N = simplify(mangalnetwork(mangal_networks[i].id))

    mangal_topology.id[i] = mangal_networks[i].id
    mangal_topology.richness[i] = richness(N)
    mangal_topology.links[i] = links(N)
    mangal_topology.connectance[i] = connectance(N)
end
