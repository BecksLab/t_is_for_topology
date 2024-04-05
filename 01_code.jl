using CSV
using DataFrames
using Distributions
using Graphs
using JLD2
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
mangal_networks = load_object("data/mangal_networks.jlds")
mangal_summary = CSV.read("data/mangal_summary.csv", DataFrame)

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
n_reps = 40 #number of reps for each model for each network

@showprogress for _ in 1:n_reps
    for i in eachindex(mangal_networks)
        
        mangal_network = render(Binary, mangal_networks[i])
        mangal_rich = richness(mangal_network)
        mangal_connectance = connectance(mangal_network)
        mangal_links = links(mangal_network)

        for (j, val) in enumerate(model_names)
            if val == "random"
                N = nullmodel(Connectance, mangal_network)
                N = randomdraws(N) # from probabalistic to binary
            elseif val == "niche"
                N = structuralmodel(NicheModel, mangal_rich, mangal_connectance)
            elseif val == "cascade"
                N = cascademodel(mangal_rich, mangal_connectance)
            elseif val == "hierarchy"
                N = nestedhierarchymodel(mangal_rich, mangal_links)
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
                D[:id] = mangal_summary.id[i]
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
