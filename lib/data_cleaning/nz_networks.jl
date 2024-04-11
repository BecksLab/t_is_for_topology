# Download and prep New-Zealand food webs for analysis

mutable struct NZ_foodweb
    abundance::Vector{Float64}
    mass::Vector{Float64}
    network::SpeciesInteractionNetwork{<:Partiteness, <:Binary}
    id::Any
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
  taxa_N = names(N_df)
  N_NZ = N_NZ[in(abund_taxa).(taxa_N), in(abund_taxa).(taxa_N)]
  sp_names = taxa_N[in(abund_taxa).(taxa_N)]

  edges = Binary(N_NZ)
  nodes =  Unipartite(sp_names)
  N = SpeciesInteractionNetwork(nodes, edges)

  # simplify network and get species names
  simplify(N)
  taxa_N = species(N)
  # get abundance data for all taxa in network
  abund = abund_fw[in(taxa_N).(abund_taxa), :no_m2]
  mass = abund_fw[in(taxa_N).(abund_taxa), :dw]
  return NZ_foodweb(abund, mass, N, fw_name)
end

# get the abundance data of all food webs in New Zealand
NZ_webs = abundance_data_NZ.(fw_names)

## Write networks as object
save_object("data/raw/new_zealand/nz_networks.jlds", NZ_webs)

# create summary of networks

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
  
  # now we calcualte relevant 'summary statistics'
  gen = SpeciesInteractionNetworks.generality(N)
  ind_maxgen = findmax(collect(values(gen)))[2] # find a species with maximal generality
  basal = findall(x -> x == 0.0, collect(values(gen))) # find species with generality of zero
  
  vul = SpeciesInteractionNetworks.vulnerability(N)
  top = findall(x -> x == 0.0, collect(values(vul))) # find species with vulnerability of zero
  
  D = Dict{Symbol, Any}()
    D[:id] = NZ_webs[i].id
    D[:richness] = richness(N)
    D[:links] = links(N)
    D[:connectance] = connectance(N)
    D[:complexity] = complexity(N)
    D[:distance] = distancetobase(N, collect(keys(gen))[ind_maxgen])
    D[:basal] = length(basal)
    D[:top] = length(top)
    D[:S1] = length(findmotif(motifs(Unipartite, 3)[1], N))
    D[:S2] = length(findmotif(motifs(Unipartite, 3)[2], N))
    D[:S4] = length(findmotif(motifs(Unipartite, 3)[4], N))
    D[:S5] = length(findmotif(motifs(Unipartite, 3)[5], N))
    push!(nz_topology, D)
end


## Write file
CSV.write("data/processed/nz_summary.csv", nz_topology)