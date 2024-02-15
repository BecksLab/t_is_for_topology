#using DataFrames
using Distributions
#using EcologicalNetworks
using Mangal
using Random
using SpeciesInteractionNetworks

# Load the functions we need from the lib folder
include("lib/adbm.jl")
include("lib/cascade.jl")
#include("lib/nestedhierarchy.jl")


mangal_networks = Mangal.datasets("type" => "predation", "count" => 20)
# 20 datasets for now just for initial testing

# there is a weird glitch that low id networks aren't being called
# so I'm just truncating it for now...
mangal_networks = mangal_networks[11:19]

mangalnetwork(19)

# nvm there are just some ids that are not a vibe...

mangal_topology = fill(0, (length(mangal_networks), 4)) 
# we can make this a nice strut later...

for i in eachindex(mangal_topology)
    N = mangalnetwork(mangal_networks[i].id)

    mangal_topology[i,1] = mangal_networks[i].id
    mangal_topology[i,2] = richness(N)
    mangal_topology[i,3] = links(N)
    mangal_topology[i,4] = connectance(N)

end

species = 10
Co = 0.2

niche_network = structuralmodel(NicheModel, species, Co)
cascade_network = cascademodel(species, Co)
hierarchy_network = nestedhierarchymodel(species, 50)

adbmm = adbmmodel(species, rand(Float64, 10), rand(Float64, 10))

