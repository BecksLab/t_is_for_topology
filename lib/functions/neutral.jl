"""
neutral_model(abundance::Vector, L::Int64)
    
    Return a network where links are assigned based on the abundance
    of species.

    Note this function has been ported from the Banville (2023)
    source code and has been (minimally) modified for the purpose of 
    this project.

    #### References

    Banville, Francis, Dominique Gravel, and Timothée Poisot. 2023. 
    “What Constrains Food Webs? A Maximum Entropy Framework for Predicting 
    Their Structure with Minimal Biases.” PLOS Computational Biology 19 (9): 
    e1011458. https://doi.org/10.1371/journal.pcbi.1011458.
"""
function neutral_model(abundance::Vector, L::Int64)
    rel_abun = abundance ./ sum(abundance) # vector of relative abundances
    abun_mat = rel_abun * rel_abun' # neutral abundance matrix (proportional to adjacency matrix)

    S = length(abundance) # species richness
    
    # draw L interactions using the neutral abundance matrix as weights
    # repeat multiple times and keep the L higher values
    N_sim = zeros(Float64, S, S)
    N_neutral =  zeros(Bool, S, S)

    for i in 1:100
        N_samp = sample(eachindex(abun_mat), Weights(vec(abun_mat)), L, replace=false)
        N_sim[N_samp] .= N_sim[N_samp] .+ 1.0
    end

    thres = sort(vec(N_sim), rev=true)[L] # Lth maximum value
    N_neutral[N_sim .< thres] .= 0.0
    N_neutral[N_sim .> thres] .= 1.0

    i = sample(findall(N_sim .== thres), Int64(L-sum(N_neutral)), replace=false)
    N_neutral[i] .= 1

    # convert to unipartite network 
    edges = Binary(N_neutral)
    nodes = Unipartite(edges)
    return SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
end