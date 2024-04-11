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
include("lib/functions/adbm.jl")
include("lib/functions/cascade.jl")
include("lib/functions/maxent.jl")
include("lib/functions/nestedhierarchy.jl")
include("lib/functions/neutral.jl")
include("lib/functions/random.jl")

# import mangal networks
mangal_networks = load_object("data/raw/mangal/mangal_networks.jlds")
mangal_summary = CSV.read("data/processed/mangal_summary.csv", DataFrame)

# create df for the outputs to be stored (long format)

topology  = DataFrame(
    id = Any[],
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

model_names = ["random", "niche", "cascade", "hierarchy", "maxent", "neutral", "adbm"]
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

# import mangal networks
nz_networks = load_object("data/raw/new_zealand/nz_networks.jlds")

@showprogress for _ in 1:n_reps
    for i in eachindex(nz_networks)
    
        network = render(Binary, nz_networks[i].network)
        abun = nz_networks[i].abundance
        mass = nz_networks[i].mass
        id = nz_networks[i].id

        for (j, val) in enumerate(model_names)
            if val == "random"
                N = nullmodel(Connectance, network)
                N = randomdraws(N) # from probabalistic to binary
            elseif val == "niche"
                N = structuralmodel(NicheModel, richness(network), connectance(network))
            elseif val == "cascade"
                N = cascademodel(richness(network), connectance(network))
            elseif val == "hierarchy"
                N = nestedhierarchymodel(richness(network), links(network))
            elseif val == "maxent"
                N = maxentmodel(network;
                # ❗ TODO
                nchains = 2,
                nsteps = 20)
            elseif  val == "neutral"
                neutral_model(abun, links(network))
            else  val == "adbm"
                parameters = adbm_parameters(network, mass)
                N = adbmmodel(richness(network), parameters, abun)
            end
        
            N = simplify(N)
            N = render(Binary, N)

            if richness(N) > 0
        
            gen = SpeciesInteractionNetworks.generality(N)
            ind_maxgen = findmax(collect(values(gen)))[2]
            basal = findall(x -> x == 0.0, collect(values(gen)))

            vul = SpeciesInteractionNetworks.vulnerability(N)
            top = findall(x -> x == 0.0, collect(values(vul))) # find species with vulnerability of zero
     
            D = Dict{Symbol, Any}()
                D[:id] = id
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
end

## ADBM networks

CSV.write("data/processed/topology_summary.csv", topology)
