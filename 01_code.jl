using Mangal
using SpeciesInteractionNetworks

# # Niche Model

niche_network = structuralmodel(NicheModel, 10, 0.2)
niche_network.edges