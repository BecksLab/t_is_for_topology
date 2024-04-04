using CSV
using DataFrames
using JLD
using Mangal
using ProgressMeter
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
    complexity = zeros(Float64, network_count),
    distance = zeros(Float64, network_count),
    basal = zeros(Float64, network_count),
    top = zeros(Float64, network_count),
    S1 = zeros(Float64, network_count),
    S2 = zeros(Float64, network_count),
    S4 = zeros(Float64, network_count),
    S5 = zeros(Float64, network_count),
); 

@showprogress for i in eachindex(mangal_networks)
    N = simplify(mangalnetwork(mangal_networks[i].id))
    N = render(Binary, N) # make binary

    gen = SpeciesInteractionNetworks.generality(N)
    ind_maxgen = findmax(collect(values(gen)))[2]
    basal = findall(x -> x == 0.0, collect(values(gen))) # find species with generality of zero

    vul = SpeciesInteractionNetworks.vulnerability(N)
    top = findall(x -> x == 0.0, collect(values(vul))) # find species with vulnerability of zero

    mangal_topology.id[i] = mangal_networks[i].id
    mangal_topology.richness[i] = richness(N)
    mangal_topology.links[i] = links(N)
    mangal_topology.connectance[i] = connectance(N)
    mangal_topology.complexity[i] = complexity(N)
    mangal_topology.distance[i] = distancetobase(N, collect(keys(gen))[ind_maxgen])
    mangal_topology.basal[i] = length(basal)
    mangal_topology.top[i] = length(top)
    mangal_topology.S1[i] = length(findmotif(motifs(Unipartite, 3)[1], N))
    mangal_topology.S2[i] = length(findmotif(motifs(Unipartite, 3)[2], N))
    mangal_topology.S4[i] = length(findmotif(motifs(Unipartite, 3)[4], N))
    mangal_topology.S5[i] = length(findmotif(motifs(Unipartite, 3)[5], N))
end

CSV.write("data/mangal_summary.csv", mangal_topology)

## Write file
CSV.write("data/mangal_networks.csv", mangal_networks)