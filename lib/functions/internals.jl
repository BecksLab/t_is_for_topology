# General sundry internal functions

"""
_get_matrix(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    Internal function to return a matrix of interactions from a
    SpeciesInteractionNetwork
"""
function _get_matrix(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    species = richness(N)
    n = zeros(Bool, (species, species))
    for i in axes(n, 1)
        for j in axes(n, 2)
            if N.edges[i,j] == true
                n[i,j] = 1
            end
        end
    end

    return n
end