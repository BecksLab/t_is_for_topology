using Distributions
using Mangal
using Random
using SpeciesInteractionNetworks

# Load the functions we need from the lib folder

include("lib/adbm.jl")
include("lib/cascade.jl")
include("lib/nestedhierarchy.jl")


species = 10
Co = 0.2

niche_network = structuralmodel(NicheModel, species, Co)
cascade_network = cascademodel(species, Co)
hierarchy_network = nestedhierarchymodel(species, 50)

adbmm = adbmmodel(species, rand(Float64, 10), rand(Float64, 10))