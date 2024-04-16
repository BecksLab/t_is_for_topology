# Download and prep New-Zealand food webs for analysis

# read taxonomic and abundance data in New Zealand (all networks)
abund_NZ = DataFrame(
    CSV.File(joinpath("data", "raw", "new_zealand", "taxa_dry_weight_abundance.csv")),
)
rename!(abund_NZ, 1 => :food_web, 4 => :no_m2)
fw_names = unique(abund_NZ[!, :food_web])

# get the name of all food webs (adjacency matrices) in New Zealand
matrix_names = readdir(joinpath("data", "raw", "new_zealand", "adjacency_matrices"))
matrix_names = replace.(matrix_names, ".csv" => "")

# struct for outputs
mutable struct NZ_foodweb
    abundance::Vector{Float64}
    mass::Vector{Float64}
    network::SpeciesInteractionNetwork{<:Partiteness,<:Binary}
    id::Any
end

"""
abundance_data_NZ(fw_name::String) 

  Returns the density (no/m2), biomass, and SpeciesInteractionNetwork of all species in
  the simplified food web with abundance data for the New Zealand dataset

  Note this function has been ported from the Banville (2023)
  source code and has been (minimally) modified for the purpose of 
  this project.
  
  #### References
  
  Banville, Francis, Dominique Gravel, and Timothée Poisot. 2023. 
  “What Constrains Food Webs? A Maximum Entropy Framework for Predicting 
  Their Structure with Minimal Biases.” PLOS Computational Biology 19 (9): 
  e1011458. https://doi.org/10.1371/journal.pcbi.1011458. 
"""
function abundance_data_NZ(fw_name::String)

    # filter abundance data for this food web and get species names
    abund_fw = abund_NZ[abund_NZ[!, :food_web].==fw_name, :]
    abund_taxa = abund_fw[!, :taxa]
    # read food web data for the specifc dataset
    N_df = DataFrame(
        CSV.File.(
            joinpath("data", "raw", "new_zealand", "adjacency_matrices", "$fw_name.csv"),
            drop = [1],
        ),
    )

    # 'fiter' networks using abundance/taxonomy combos
    N_NZ = Matrix{Bool}(N_df)
    taxa_N = names(N_df)
    N_NZ = N_NZ[in(abund_taxa).(taxa_N), in(abund_taxa).(taxa_N)]
    sp_names = taxa_N[in(abund_taxa).(taxa_N)]

    # turn into network
    edges = Binary(N_NZ)
    nodes = Unipartite(sp_names)
    N = SpeciesInteractionNetwork(nodes, edges)

    # simplify network and get species names
    simplify(N)
    taxa_N = species(N)
    # get abundance data for all taxa in network
    abund = abund_fw[in(taxa_N).(abund_taxa), :no_m2]
    # get biomass data for all taxa in network
    mass = abund_fw[in(taxa_N).(abund_taxa), :dw]

    return NZ_foodweb(abund, mass, N, fw_name)
end

# get the abundance data of all food webs in New Zealand
NZ_webs = abundance_data_NZ.(fw_names)

## Write networks as object
save_object("data/raw/new_zealand/nz_networks.jlds", NZ_webs)

# Create summary of networks

# make a nice 'dataframe' to store network data
nz_topology = DataFrame(
    id = Any[],
    richness = Int64[],
    links = Int64[],
    connectance = Float64[],
    complexity = Float64[],
    distance = Float64[],
    basal = Float64[],
    top = Float64[],
    S1 = Float64[],
    S2 = Float64[],
    S4 = Float64[],
    S5 = Float64[],
);

@showprogress for i in eachindex(NZ_webs)

    N = simplify(NZ_webs[i].network)
    N = render(Binary, N) # make binary

    d = _network_summary(N)

    D = Dict{Symbol,Any}()
    D[:id] = NZ_webs[i].id
    D[:richness] = d[:richness]
    D[:links] = d[:links]
    D[:connectance] = d[:connectance]
    D[:complexity] = d[:complexity]
    D[:distance] = d[:distance]
    D[:basal] = d[:basal]
    D[:top] = d[:top]
    D[:S1] = d[:S1]
    D[:S2] = d[:S2]
    D[:S4] = d[:S4]
    D[:S5] = d[:S5]
    push!(nz_topology, D)
end


## Write file
CSV.write("data/processed/nz_summary.csv", nz_topology)
