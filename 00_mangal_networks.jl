using CSV
using DataFrames
using Mangal
using SpeciesInteractionNetworks

network_count = 3
mangal_networks = Mangal.networks("type" => "herbivory", "count" => network_count)
# 20 datasets for now just for initial testing

# make a nice 'dataframe' to store network data
mutable struct MangalNetworks
    id::Vector{Int64}
    richness::Vector{Int64}
    links::Vector{Int64}
    connectance::Vector{Float64}
    # we can add some 'summary stats' here as well
end

function MangalNetworks(C::Int64)
    return MangalNetworks(
        zeros(Int64, C),
        zeros(Int64, C),
        zeros(Int64, C),
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

df = DataFrame(
    id = mangal_topology.id,
    richness = mangal_topology.richness,
    links = mangal_topology.links,
    connectance = mangal_topology.connectance
)

CSV.write("data/mangal_summary.csv", df)