"""
    nestedhierarchymodel(species::Int64, connectance::Int64)

    Return a network where resources are assigned to consumers according
    to the nested hierarchy model for the specified number of species and
    links.

    Note this function has been ported from the EcologicalNetworks.jl
    source code and has been (minimally) modified for the purpose of 
    this project.
    
    #### References
    
    Cattin, M.-F., Bersier, L.-F., Banašek-Richter, C., Baltensperger, R., Gabriel,
    J.-P., 2004. Phylogenetic constraints and adaptation explain food-web structure.
    Nature 427, 835–839. https://doi.org/10.1038/nature02327
    
    Poisot, Timothée, Zacharie Belisle, Laura Hoebeke, Michiel Stock, and Piotr
    Szefer. 2019. “EcologicalNetworks.jl - Analysing Ecological Networks.” Ecography.
    https://doi.org/10.1111/ecog.04310.
"""
function nestedhierarchymodel(S::Int64, L::Int64)

    # Initiate matrix
    edges = Binary(zeros(Bool, (S, S)))
    nodes = Unipartite(Binary(edgs))
    A = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)

    Co = (L/(S*S))

    # Assign values to all species and sort in an ascending order
    e = sort(rand(S))

    # Beta parameter (after Allesina et al. 2008)
    # Classic niche: β = 1.0/(2.0*C)-1.0
    β = (S - 1.0)/(2.0*Co*S) - 1.0

    # Random variable from the Beta distribution
    X = rand(Beta(1.0, β), 1)

    # 'r' values are assigned to all species
    r = e .* X

    # Number of prey for predator i, li rounded to the nearest integer
    expectedlinks = ((r ./ sum(r)) .* L)
    observedlinks = map(x -> round(Int, x), expectedlinks)

    totallinks = sum(observedlinks)

    # User message: number of links to be assigned
    # println("Nested hierarchy model with ", totallinks, " links to assign")

    # Links which may happen to exceed S are set to S-1
    for link in 1:S

        if observedlinks[link] >= S

            observedlinks[link] = S-1

        end

    end

    # Ensure that the first species has no links
    observedlinks[1] = 0

    # Species pool
    species_pool = species(A)

    # Create dictionaries with all the species indices and names

    # Strings into indices
    sp_pos = Dict{Symbol, Int64}()

    for (n, f) in enumerate(species_pool)

        sp_pos[f] = n

    end

    # Indices into strings
    sp_str = Dict{Int64, Symbol}()

    for (n, f) in enumerate(species_pool)

        sp_str[n] = f

    end

    # Species already predated
    predated = []

    # Assign links to all consumers except for the first one (basal)
    for consumer in 2:S

        # STAGE 1: Assign species with smaller niche value

        # Assign species randomly unless one has other predators
        eligible = StatsBase.sample(1:(consumer-1), (consumer-1), replace = false)

        linkstoassign = observedlinks[consumer]

        # Hold the last resource assagned in at this stage
        last_resource = 0

        # Loop through a random resources
        for resource in eligible

            # Update last resource holder
            last_resource = resource

            # Break the loop in case of reaching desired number of links
            linkstoassign == 0 && break

            # Assign link
            A[consumer,resource] = true
            linkstoassign = linkstoassign - 1

            # Update the predatored species array
            push!(predated, sp_str[resource])

            # If there are other predators break this
            (sum(A[:,resource]) > 1) && break

        end

        # Go to another consumer in case there are no more links
        linkstoassign == 0 && continue

        # STAGE 2: if the previos resource had a consumer assign resource from th eunion

        # Make a union of all consumers feeding on last_resource
        shared_consumers = predecessors(A, species_pool[last_resource])
        resource_union = Set{Symbol}()
        [union!(resource_union, predecessors(A, cons)) for cons in shared_consumers]

        # Remove species already utilized by current consumer
        resource_union = setdiff(resource_union, predecessors(A, species_pool[consumer]))

        # Shuffle the resource_union
        resource_union = StatsBase.sample(collect(resource_union), length(resource_union), replace = false)

        for resource in resource_union

            linkstoassign == 0 && break

            A[consumer, sp_pos[resource]] = true
            linkstoassign = linkstoassign - 1

            # Update the no_predator array
            push!(predated, resource)

        end

        # Go to another consumer in case there are no more links
        linkstoassign == 0 && continue

        # STAGE 3: randomly assign resources without predators

        # Remove already utilized species
        not_predated = setdiff(species_pool, predated)

        # Shuffle the current set
        not_predated = StatsBase.sample(collect(not_predated), length(not_predated), replace = false)

        # Assign links
        for resource in not_predated

            linkstoassign == 0 && break

            A[consumer, sp_pos[resource]] = true
            linkstoassign = linkstoassign - 1

        end

        # Go to another consumer in case there are no more links
        linkstoassign == 0 && continue

        # STAGE 4: Randomly assign species from the rest of the resources

        # Rest of the unassigned species!
        diet = collect(A.edges[consumer])

        # Get rid of specis already assigned to the current consumer
        not_consumed = setdiff(species_pool, diet)

        # Shuffle the set
        not_consumed = StatsBase.sample(collect(not_consumed), length(not_consumed), replace = false)

        # Assign links
        for resource in not_consumed

            linkstoassign == 0 && break

            A.edges[consumer, sp_pos[resource]] = true
            linkstoassign = linkstoassign - 1

        end

    end

    return A

end