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
include("lib/functions/confusionmatrix.jl")

# Load networks
observed_networks = load_object("data/raw/mangal/mangal_networks.jlds")
nz_networks = load_object("data/raw/new_zealand/nz_networks.jlds")
predicted_networks = load_object("data/processed/predicted_networks.jlds")

# we're going to have to reformat the nz networks a bit... maybe
# but we can append to the mangal data
for i in eachindex(nz_networks)

    N_d = Dict{Symbol,Any}()
    N_d[:id] = nz_networks[i].id
    N_d[:network] = nz_networks[i].network
    push!(observed_networks, N_d)
    
end