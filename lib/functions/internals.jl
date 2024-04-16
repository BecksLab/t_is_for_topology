# General sundry internal functions

"""
_get_matrix(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    Internal function to return a matrix of interactions from a
    SpeciesInteractionNetwork
"""
function _get_matrix(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    species = richness(N)
    n = zeros(Int64, (species, species))
    for i in axes(n, 1)
        for j in axes(n, 2)
            if N.edges[i, j] == true
                n[i, j] = 1
            end
        end
    end

    return n
end

"""
_network_summary(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    returns the 'summary statistics' for a network
"""
function _network_summary(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    A = _get_matrix(N)

    gen = SpeciesInteractionNetworks.generality(N)
    ind_maxgen = findmax(collect(values(gen)))[2]

    D = Dict{Symbol,Any}(
        :richness => richness(N),
        :links => links(N),
        :connectance => connectance(N),
        :complexity => complexity(N),
        :distance => distancetobase(N, collect(keys(gen))[ind_maxgen]),
        :basal => sum(vec(sum(A, dims = 2) .== 0)),
        :top => sum(vec(sum(A, dims = 1) .== 0)),
        :S1 => length(findmotif(motifs(Unipartite, 3)[1], N)),
        :S2 => length(findmotif(motifs(Unipartite, 3)[2], N)),
        :S4 => length(findmotif(motifs(Unipartite, 3)[4], N)),
        :S5 => length(findmotif(motifs(Unipartite, 3)[5], N)),
    )

    return D
end

"""
model_summary(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    returns the 'summary statistics' for a network
"""
function model_summary(
    network::SpeciesInteractionNetwork{<:Partiteness,<:Binary},
    network_id::Any,
    model_name::String;
    abundance::Vector{Float64} = [0.0, 0.0],
    mass::Vector{Float64} = [0.0, 0.0],
)

    # data checks
    if model_name ∉ ["random", "niche", "cascade", "hierarchy", "maxent", "neutral", "adbm"]
        error(
            "Invalid value for model_name -- must be one of random, niche, cascade, hierarchy, maxent, neutral, adbm",
        )
    end
    if model_name == "adbm" && length(abundance) != richness(network)
        error("Invalid length for abundance -- must be length $(richness(network))")
    end
    if model_name == "neutral" && length(abundance) != richness(network)
        error("Invalid length for abundance -- must be length $(richness(network))")
    end
    if model_name == "adbm" && length(mass) != richness(network)
        error("Invalid length for mass -- must be length $(richness(network))")
    end

    # generate network based on specified model
    if model_name == "random"
        N = nullmodel(Connectance, network)
        N = randomdraws(N) # from probabalistic to binary
    elseif model_name == "niche"
        N = structuralmodel(NicheModel, richness(network), connectance(network))
    elseif model_name == "cascade"
        N = cascademodel(richness(network), connectance(network))
    elseif model_name == "hierarchy"
        N = nestedhierarchymodel(richness(network), links(network))
    elseif model_name == "maxent"
        N = maxentmodel(
            network;
            # ❗ TODO
            nchains = 2,
            nsteps = 20,
        )
    elseif model_name == "neutral"
        N = neutralmodel(abundance, links(network), species(network))
    else
        model_name == "adbm"
        parameters = adbm_parameters(network, mass)
        N = adbmmodel(species(network), parameters, abundance)
    end

    N = simplify(N)
    N = render(Binary, N)

    if richness(N) > 0
        d = _network_summary(N)
        D = Dict{Symbol,Any}()
        D[:id] = network_id
        D[:model] = model_name
        D[:connectance_mod] = d[:connectance]
        D[:complexity_mod] = d[:complexity]
        D[:distance_mod] = d[:distance]
        D[:basal_mod] = d[:basal]
        D[:top_mod] = d[:top]
        D[:S1_mod] = d[:S1]
        D[:S2_mod] = d[:S2]
        D[:S4_mod] = d[:S4]
        D[:S5_mod] = d[:S5]
    else
        D = nothing
    end
    return D
end
