
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

        push!(networks, N) # push network 'as is'

        N = render(Binary, N) # make binary

        d = _network_summary(N)

        D = Dict{Symbol,Any}()
        D[:id] = mangal_networks[i].id
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
        push!(mangal_topology, D)
    end
end

## Write networks as object
save_object("data/raw/mangal/mangal_networks.jlds", networks)

## Write file
CSV.write("data/raw/mangal/mangal_networks.csv", mangal_networks)
CSV.write("data/processed/mangal_summary.csv", mangal_topology)
