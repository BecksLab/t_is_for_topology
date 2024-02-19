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

model_names = ["random", "niche", "cascade"]

# create df for the outputs to be stored (long format)

topology  = DataFrame(
    id = Int64[],
    richness = Int64[],
    links = Int64[],
    connectance_real = Float64[],
    complexity_real = Float64[],
    distance_real = Float64[],
    model = String[],
    connectance_mod = Float64[],
    complexity_mod = Float64[],
    distance_mod = Float64[]
);


    for i in 1:(nrow(mangal_topology))
        
        for (j, val) in enumerate(model_names)
            if val == "random"
                N = randommodel(mangal_topology.richness[i], mangal_topology.links[i])
            elseif val == "niche"
                N = structuralmodel(NicheModel, mangal_topology.richness[i], mangal_topology.connectance[i])
            else val == "cascade"
                N = cascademodel(mangal_topology.richness[i], mangal_topology.connectance[i])
            end
            
        # specificity
        spe = specificity(N)
        ind_maxspe = findmax(collect(values(spe)))[2]
    
        D = Dict{Symbol, Any}()
            D[:id] = mangal_topology.id[i]
            D[:richness] = mangal_topology.richness[i]
            D[:links] = mangal_topology.links[i]
            D[:connectance_real] = mangal_topology.connectance[i]
            D[:complexity_real] = mangal_topology.complexity[i]
            D[:distance_real] = mangal_topology.distance[i]
            D[:model] = val
            D[:connectance_mod] = connectance(N)
            D[:complexity_mod] = complexity(N)
            D[:distance_mod] = distancetobase(N, collect(keys(spe))[ind_maxspe])
            push!(topology, D)
    end  
end

CSV.write("data/topology_summary.csv", topology)
