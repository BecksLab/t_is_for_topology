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
#include("lib/nestedhierarchy.jl")
include("lib/random.jl")

# import mangal networks
mangal_topology = CSV.read("data/mangal_summary.csv", DataFrame)

# create df for the outputs to be stored (long format)

topology  = DataFrame(
    id = Int64[],
    richness = Int64[],
    links = Int64[],
    connectance_real = Float64[],
    complexity_real = Float64[],
    distance_real = Float64[],
    basal_real = Float64[],
    model = String[],
    connectance_mod = Float64[],
    complexity_mod = Float64[],
    distance_mod = Float64[],
    basal_mod = Float64[]
);

## Structural networks

model_names = ["random", "niche", "cascade"]

    for i in 1:(nrow(mangal_topology))
        
        for (j, val) in enumerate(model_names)
            if val == "random"
                N = randommodel(mangal_topology.richness[i], mangal_topology.links[i])
            elseif val == "niche"
                N = structuralmodel(NicheModel, mangal_topology.richness[i], mangal_topology.connectance[i])
            else val == "cascade"
                N = cascademodel(mangal_topology.richness[i], mangal_topology.connectance[i])
            end
            
        spe = specificity(N)
        gen = SpeciesInteractionNetworks.generality(N)
        ind_maxgen = findmax(collect(values(gen)))[2]
        basal = findall(x -> x == 0.0, collect(values(gen)))
        #findmax(collect(values(spe)))
    
        D = Dict{Symbol, Any}()
            D[:id] = mangal_topology.id[i]
            D[:richness] = mangal_topology.richness[i]
            D[:links] = mangal_topology.links[i]
            D[:connectance_real] = mangal_topology.connectance[i]
            D[:complexity_real] = mangal_topology.complexity[i]
            D[:distance_real] = mangal_topology.distance[i]
            D[:basal_real] = mangal_topology.basal[i]
            D[:model] = val
            D[:connectance_mod] = connectance(N)
            D[:complexity_mod] = complexity(N)
            D[:distance_mod] = distancetobase(N, collect(keys(gen))[ind_maxgen])
            D[:basal_mod] = length(basal)/richness(N)
            push!(topology, D)
    end  
end

## ADBM networks

CSV.write("data/topology_summary.csv", topology)
