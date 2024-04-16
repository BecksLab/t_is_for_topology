using RandomBooleanMatrices

include("internals.jl")

"""
maxentmodel(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}; 
        nchains::Int64 = 4, nsteps::Int64 = 2000)

    Returns the adjacency matrix of maximum SVD-entropy constrained 
    by the joint degree sequence of a networks using a simulating 
    annealing algorithm

    Note this function has been ported from the Banville (2023)
    source code and has been (minimally) modified for the purpose of 
    this project.

    #### References

    Banville, Francis, Dominique Gravel, and Timothée Poisot. 2023. 
    “What Constrains Food Webs? A Maximum Entropy Framework for Predicting 
    Their Structure with Minimal Biases.” PLOS Computational Biology 19 (9): 
    e1011458. https://doi.org/10.1371/journal.pcbi.1011458. 
"""
function maxentmodel(
    N::SpeciesInteractionNetwork{<:Partiteness,<:Binary};
    nchains::Int64 = 4,
    nsteps::Int64 = 2000,
)

    # matrix generator object
    n = _get_matrix(N)

    rmg = matrixrandomizer(n)
    # initial vector for SVD-entropies and object for best matrices
    entropies = zeros(Float64, nsteps, nchains)
    A = []

    for j = 1:nchains
        # generate a new random matrix with the same row and column sums (joint degree sequence) as N
        edges = Binary(rand(rmg))
        nodes = Unipartite(edges)
        A0 = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
        A0 = simplify(A0)
        best = complexity(A0)
        entropies[1, j] = best

        # initial temperature of the simulating annealing algorithm
        T = 0.2

        # simulating annealing algorithm
        for i = 2:nsteps
            # propose a new constrained random matrix and compute the difference in SVD-entropy 
            A1 = swap!(A0)
            candidate = complexity(A1)
            delta = candidate - best
            # accept if the difference is positive or with a probability p if it's negative
            if delta > 0
                A0 = A1
                best = candidate
            else
                p = exp(delta / T)
                P = rand(Uniform(0, 1))
                if P < p
                    A0 = A1
                    best = candidate
                end
            end
            entropies[i, j] = best
            # update the temperature
            T = T * 0.99
        end
        push!(A, A0)
    end
    # find the network with maximum entropy among all chains
    imax = findmax(entropies[nsteps, :])[2]
    Amax = (A = A[imax], entropies = entropies)
    return simplify(Amax.A)
end
