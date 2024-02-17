using CSV
using DataFrames
using Mangal
using SpeciesInteractionNetworks

network_count = 20
mangal_networks = Mangal.networks("type" => "herbivory", "count" => network_count)
# 20 datasets for now just for initial testing

# make a nice 'dataframe' to store network data
mangal_topology = DataFrame(
    id = zeros(Int64, network_count),
    richness = zeros(Int64, network_count),
    links = zeros(Int64, network_count),
    connectance = zeros(Float64, network_count),
    complexity = zeros(Float64, network_count)
); 

for i in eachindex(mangal_networks)
    N = simplify(mangalnetwork(mangal_networks[i].id))

    mangal_topology.id[i] = mangal_networks[i].id
    mangal_topology.richness[i] = richness(N)
    mangal_topology.links[i] = links(N)
    mangal_topology.connectance[i] = connectance(N)
    mangal_topology.complexity[i] = complexity(N)
end

CSV.write("data/mangal_summary.csv", mangal_topology)