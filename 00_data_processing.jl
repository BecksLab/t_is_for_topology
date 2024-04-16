using CSV
using DataFrames
using Glob
using JLD2
using Mangal
using ProgressMeter
using SpeciesInteractionNetworks

# Mangal networks

include("lib/data_cleaning/mangal_networks.jl")

# NZ networks

include("lib/data_cleaning/nz_networks.jl")
