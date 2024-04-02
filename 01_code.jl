using CSV
using DataFrames
using Distributions
using Graphs
#using EcologicalNetworks
using Mangal
using Random
using RandomBooleanMatrices
using SpeciesInteractionNetworks
using StatsBase

# Load the functions we need from the lib folder
include("lib/adbm.jl")
include("lib/cascade.jl")
include("lib/maxent.jl")
include("lib/nestedhierarchy.jl")
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
    top_real = Float64[],
    model = String[],
    connectance_mod = Float64[],
    complexity_mod = Float64[],
    distance_mod = Float64[],
    basal_mod = Float64[],
    top_mod = Float64[]
);

## Structural networks

model_names = ["random", "niche", "cascade", "hierarchy", "maxent"]
n_reps = 10 #number of reps for each model for each network

for _ in 1:n_reps
    for i in 1:(nrow(mangal_topology))
        
        for (j, val) in enumerate(model_names)
            if val == "random"
                N = randommodel(mangal_topology.richness[i], mangal_topology.links[i])
            elseif val == "niche"
                N = structuralmodel(NicheModel, mangal_topology.richness[i], mangal_topology.connectance[i])
            elseif val == "cascade"
                N = cascademodel(mangal_topology.richness[i], mangal_topology.connectance[i])
            elseif val == "hierarchy"
                N = nestedhierarchymodel(mangal_topology.richness[i], mangal_topology.links[i])
            else val == "maxent"
                N = maxentmodel(mangal_topology.richness[i], mangal_topology.connectance[i];
                # ❗ TODO
                nchains = 2,
                nsteps = 20)
            end

        N = simplify(N)
            
        gen = SpeciesInteractionNetworks.generality(N)
        ind_maxgen = findmax(collect(values(gen)))[2]
        basal = findall(x -> x == 0.0, collect(values(gen)))

        vul = SpeciesInteractionNetworks.vulnerability(N)
        ind_minvul = findmin(collect(values(vul)))[2]
    
        D = Dict{Symbol, Any}()
            D[:id] = mangal_topology.id[i]
            D[:richness] = mangal_topology.richness[i]
            D[:links] = mangal_topology.links[i]
            D[:connectance_real] = mangal_topology.connectance[i]
            D[:complexity_real] = mangal_topology.complexity[i]
            D[:distance_real] = mangal_topology.distance[i]
            D[:basal_real] = mangal_topology.basal[i]
            D[:top_real] = mangal_topology.top[i]
            D[:model] = val
            D[:connectance_mod] = connectance(N)
            D[:complexity_mod] = complexity(N)
            D[:distance_mod] = distancetobase(N, collect(keys(gen))[ind_maxgen])
            D[:basal_mod] = length(basal)
            D[:top_mod] = length(top)
            push!(topology, D)
        end  
    end
end

## ADBM networks

CSV.write("data/topology_summary.csv", topology)
