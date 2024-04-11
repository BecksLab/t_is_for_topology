
network_count = 20
mangal_networks = Mangal.networks("type" => "herbivory", "count" => network_count)
# 20 datasets for now just for initial testing

# make a nice 'dataframe' to store network data
mangal_topology = DataFrame(
    id = Int64[],
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

# make object to store each network so we can import it later
networks = Any[]

@showprogress for i in eachindex(mangal_networks)
    N = simplify(mangalnetwork(mangal_networks[i].id))
    
    if richness(N) <= 200 # remove large networks because they make things slow
        
        push!(networks, N) # push netwrork 'as is'

        N = render(Binary, N) # make binary

        # now we calcualte relevant 'summary statistics'
        gen = SpeciesInteractionNetworks.generality(N)
        ind_maxgen = findmax(collect(values(gen)))[2] # find a species with maximal generality
        basal = findall(x -> x == 0.0, collect(values(gen))) # find species with generality of zero
    
        vul = SpeciesInteractionNetworks.vulnerability(N)
        top = findall(x -> x == 0.0, collect(values(vul))) # find species with vulnerability of zero
    
        D = Dict{Symbol, Any}()
            D[:id] = mangal_networks[i].id
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
        push!(mangal_topology, D)
    end
end

## Write networks as object
save_object("data/raw/mangal/mangal_networks.jlds", networks)

## Write file
CSV.write("data/raw/mangal/mangal_networks.csv", mangal_networks)
CSV.write("data/processed/mangal_summary.csv", mangal_topology)