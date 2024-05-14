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
mangal_networks = load_object("data/raw/mangal/mangal_networks.jlds")
nz_networks = load_object("data/raw/new_zealand/nz_networks.jlds")

