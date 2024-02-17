using CSV
using DataFrames
using Distributions
#using EcologicalNetworks
using Mangal
using Random
using SpeciesInteractionNetworks

# Load the functions we need from the lib folder
include("lib/adbm.jl")
include("lib/cascade.jl")
#include("lib/nestedhierarchy.jl")

# create new columns for models
topology = CSV.read("data/mangal_summary.csv", DataFrame)
insertcols!(topology, 5, :niche_complexity => zeros(Float64, nrow(topology)))
insertcols!(topology, 6, :cascade_complexity => zeros(Float64, nrow(topology)))

# calcualte complexity (for now) for relevant models
for i in 1:nrow(topology)
    
    topology.niche_complexity[i] = complexity(structuralmodel(NicheModel, topology.richness[i], topology.connectance[i]))
    topology.cascade_complexity[i] = complexity(cascademodel(topology.richness[i], topology.connectance[i]))
    
end

CSV.write("data/topology_summary.csv", topology)
