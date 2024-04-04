using CSV
using DataFrames
using Distributions
using Graphs
#using EcologicalNetworks
using Mangal
using ProgressMeter
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
mangal_networks = DataFrame(CSV.File("data/mangal_networks.csv"))

# create df for the outputs to be stored (long format)

topology  = DataFrame(
    id = Int64[],
    model = String[],
    connectance_mod = Float64[],
    complexity_mod = Float64[],
    distance_mod = Float64[],
    basal_mod = Float64[],
    top_mod = Float64[],
    S1_mod = Float64[],
    S2_mod = Float64[],
    S4_mod = Float64[],
    S5_mod = Float64[]
);

## Structural networks

model_names = ["random", "niche", "cascade", "hierarchy", "maxent"]
n_reps = 10 #number of reps for each model for each network

for _ in 1:n_reps
    @showprogress for i in 1:(nrow(mangal_topology))
        
        mangal_network = simplify(mangalnetwork(mangal_networks.id[i]))
        mangal_network = render(Binary, mangal_network)

        for (j, val) in enumerate(model_names)
            if val == "random"
                N = nullmodel(Connectance, mangal_network)
                N = randomdraws(N)
            elseif val == "niche"
                N = structuralmodel(NicheModel, mangal_topology.richness[i], mangal_topology.connectance[i])
            elseif val == "cascade"
                N = cascademodel(mangal_topology.richness[i], mangal_topology.connectance[i])
            elseif val == "hierarchy"
                N = nestedhierarchymodel(mangal_topology.richness[i], mangal_topology.links[i])
            else val == "maxent"
                N = maxentmodel(mangal_network;
                # ❗ TODO
                nchains = 2,
                nsteps = 20)
            end
        
            N = simplify(N)
            N = render(Binary, N)
        
            gen = SpeciesInteractionNetworks.generality(N)
            ind_maxgen = findmax(collect(values(gen)))[2]
            basal = findall(x -> x == 0.0, collect(values(gen)))

            vul = SpeciesInteractionNetworks.vulnerability(N)
            top = findall(x -> x == 0.0, collect(values(vul))) # find species with vulnerability of zero
     
            D = Dict{Symbol, Any}()
                D[:id] = mangal_topology.id[i]
                D[:model] = val
                D[:connectance_mod] = connectance(N)
                D[:complexity_mod] = complexity(N)
                D[:distance_mod] = distancetobase(N, collect(keys(gen))[ind_maxgen])
                D[:basal_mod] = length(basal)
                D[:top_mod] = length(top)
                D[:S1_mod] = length(findmotif(motifs(Unipartite, 3)[1], N))
                D[:S2_mod] = length(findmotif(motifs(Unipartite, 3)[2], N))
                D[:S4_mod] = length(findmotif(motifs(Unipartite, 3)[4], N))
                D[:S5_mod] = length(findmotif(motifs(Unipartite, 3)[5], N))
            push!(topology, D)
        end  
    end
end

## Neutral networks

## ADBM networks

CSV.write("data/topology_summary.csv", topology)
