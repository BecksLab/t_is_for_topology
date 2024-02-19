using CSV
using DataFrames
using Distributions
using Graphs
#using EcologicalNetworks
using Mangal
using Random
using SpeciesInteractionNetworks

# Load the functions we need from the lib folder
include("lib/adbm.jl")
include("lib/cascade.jl")
include("lib/random.jl")
#include("lib/nestedhierarchy.jl")

# import mangal networks
mangal_topology = CSV.read("data/mangal_summary.csv", DataFrame)

# creat df for the outputs to be stored (long format)
topology = 

# calcualte complexity (for now) for relevant models
for i in 1:nrow(topology)
    
    # create networks
    N_niche = structuralmodel(NicheModel, topology.richness[i], topology.connectance[i])
    N_cascade = cascademodel(topology.richness[i], topology.connectance[i])

    # specificity
    spe_niche = specificity(N_niche)
    ind_maxspe_niche = findmax(collect(values(spe_niche)))[2]
    spe_cascade = specificity(N_cascade)
    ind_maxspe_cascade = findmax(collect(values(spe_cascade)))[2]

    topology.niche_complexity[i] = complexity(N_niche)
    topology.cascade_complexity[i] = complexity(N_cascade)
    topology.niche_distance[i] = distancetobase(N_niche, collect(keys(spe_niche))[ind_maxspe_niche])
    topology.cascade_distance[i] = distancetobase(N_cascade, collect(keys(spe_cascade))[ind_maxspe_cascade])
    
end


CSV.write("data/topology_summary.csv", topology)
