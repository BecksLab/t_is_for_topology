using CSV
using DataFrames
using Glob
using JLD2
using Mangal
using ProgressMeter
using SpeciesInteractionNetworks

# helper functions
include("lib/functions/internals.jl")

# Mangal networks

include("lib/data_cleaning/mangal_networks.jl")

# NZ networks

include("lib/data_cleaning/nz_networks.jl")
