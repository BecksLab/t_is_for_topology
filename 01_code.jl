using Distributions
using Mangal
using SpeciesInteractionNetworks

# # Niche Model

niche_network = structuralmodel(NicheModel, 10, 0.2)
niche_network.edges

adbmm = adbmmodel(10, rand(Float64, 10), rand(Float64, 10))