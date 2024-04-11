# Download and prep New-Zealand food webs for analysis

mutable struct NZ_foodweb
    abundance::Vector{Float64}
    mass::Vector{Float64}
    network::SpeciesInteractionNetwork{<:Partiteness, <:Binary}
end

# Abundance data

# read taxonomic and abundance data in New Zealand (all networks)
abund_NZ = DataFrame(CSV.File(joinpath("data", "raw", "new_zealand", "taxa_dry_weight_abundance.csv")))
rename!(abund_NZ , 1 => :food_web, 4 => :no_m2)
fw_names = unique(abund_NZ[!, :food_web])

# get the name of all food webs (adjacency matrices) in New Zealand
matrix_names = readdir(joinpath("data", "raw", "new_zealand", "adjacency_matrices"))
matrix_names = replace.(matrix_names, ".csv" => "")

"""
abundance_data_NZ(fw_name::String, N_NZ) 
    fw_name: name of the food web (file name)
Returns the density (no/m2) of all species in the simplified food web with abundance data
"""
function abundance_data_NZ(fw_name)

  # filter abundance data for this food web and get species names
  abund_fw = abund_NZ[abund_NZ[!, :food_web] .== fw_name, :]
  abund_taxa = abund_fw[!, :taxa]
  # read food web data and get species names
  N_df = DataFrame(CSV.File.(joinpath("data", "raw", "new_zealand", "adjacency_matrices", "$fw_name.csv"),
                    drop=[1]))
                    
  N_NZ = Matrix{Bool}(N_df)
  edges = Binary(N_NZ)
  nodes =  Unipartite(names(N_df))
  N = SpeciesInteractionNetwork(nodes, edges)

  # simplify network and get species names
  simplify(N)
  taxa_N = species(N)
  # get abundance data for all taxa in network
  abund = abund_fw[in(taxa_N).(abund_taxa), :no_m2]
  mass = abund_fw[in(taxa_N).(abund_taxa), :dw]
  return NZ_foodweb(abund, mass, N)
end

# get the abundance data of all food webs in New Zealand
NZ_webs = abundance_data_NZ.(fw_names)

## Write networks as object
save_object("data/raw/new_zealand/nz_networks.jlds", NZ_webs)
## Write files
#TODO