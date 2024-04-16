using CSV
using DataFrames
using Distributions
using JLD2
#using EcologicalNetworks
using Mangal
using ProgressMeter
using Random
using SpeciesInteractionNetworks
using StatsBase

# Load the functions we need from the lib folder
include("lib/functions/adbm.jl")
include("lib/functions/cascade.jl")
include("lib/functions/internals.jl")
include("lib/functions/maxent.jl")
include("lib/functions/nestedhierarchy.jl")
include("lib/functions/neutral.jl")
include("lib/functions/random.jl")

# import mangal networks
mangal_networks = load_object("data/raw/mangal/mangal_networks.jlds")
mangal_summary = CSV.read("data/processed/mangal_summary.csv", DataFrame)

# create df for the outputs to be stored

topology = DataFrame(
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
    S5_mod = Float64[],
);

## Mangal networks

model_names = ["random", "niche", "cascade", "maxent"]
n_reps = 40 #number of reps for each model for each network

@showprogress for _ = 1:n_reps
    for i in eachindex(mangal_networks)

        network = render(Binary, mangal_networks[i])
        abundance = [0.0, 0.0] #dummy variable
        mass = [0.0, 0.0] #dummy variable
        id = mangal_summary.id[i]

        for (j, val) in enumerate(model_names)

            D = model_summary(network, id, val; abundance, mass)
            push!(topology, D)

        end
    end
end

## New Zealand networks

# import new zealand networks
nz_networks = load_object("data/raw/new_zealand/nz_networks.jlds")
model_names = ["random", "niche", "cascade", "maxent", "neutral", "adbm"]

@showprogress for _ = 1:n_reps
    for i in eachindex(nz_networks)

        network = render(Binary, nz_networks[i].network)
        abundance = nz_networks[i].abundance
        mass = nz_networks[i].mass
        id = nz_networks[i].id

        for (j, val) in enumerate(model_names)

            D = model_summary(network, id, val; abundance, mass)
            if typeof(D) != Nothing
                push!(topology, D)
            end

        end
    end
end

CSV.write("data/processed/topology_summary.csv", topology)
