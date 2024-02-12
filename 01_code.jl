using Distributions
using Mangal
using Random
using SpeciesInteractionNetworks

# # Niche Model

niche_network = structuralmodel(NicheModel, 10, 0.2)
cascade_network = cascademodel(10, 0.2)

adbmm = adbmmodel(10, rand(Float64, 10), rand(Float64, 10))