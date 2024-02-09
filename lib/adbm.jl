"""
    _set_attack_constant

internal for adbm function, sets attack constant for species
"""
function _set_attack_constant(
    attack_constant::Vector{Float64}
)
    a = rand(-1.0:0.25:1.0, length(attack_constant))
    for i in axes(attack_constant)
        attack_constant[i] = a[i]
    end
    return attack_constant
end 

"""
        _set_energy_content

internal for adbm function, sets energy content for species
"""
function _set_energy_content(
    energy_content::Vector{Float64},
    mass::Vector{Float64}
)
    for i in axes(energy_content)
        energy_content[i] = mass[i]
    end
    return energy_content
end   


"""
    adbmmodel

function to generate netwrok using ADBM, Petchy and friends 2008
"""
function adbmmodel(
    species_richness::Int64,
    mass::Vector{Float64},
    resource_density::Vector{Float64};
    breadth::AbstractRange = -8.0:1.0:1.0,
    resource_constant::Float64 = -0.75,
    energy_constant::Float64 = 1.0,
    handling_constant::Float64 = 1.0,
)
    # create matrices
    attack_constant = zeros(Float64, species_richness)
    adbmmatrix = zeros(Float64, (species_richness, species_richness))
    energy_content = zeros(Float64, species_richness)

    # create attack constant
    attack_constant = _set_attack_constant(attack_constant)

    # energy content based on mass
    energy_content = _set_energy_content(energy_content, mass)

    b = rand(breadth)
    b = 2^b


    for j in axes(adbmmatrix, 1)
        for i in axes(adbmmatrix, 2)
            if mass[i]/mass[j] ≥ b #immediately set to zero
                adbmmatrix[j,i] = 0.0  
            end
            _handling_time = handling_constant/(b - (mass[i]/mass[j]))
            _resource_density = resource_density[i]^resource_constant
            _energy_content = energy_content[i]*energy_constant
            _attack_rate = (mass[i]^attack_constant[i])*(mass[j]^attack_constant[j])
            _encounter_rate = _resource_density*_attack_rate
            
            adbmmatrix[j,i] = (_encounter_rate*_energy_content)/(_resource_density*_handling_time)
        end
    end
    return adbmmatrix
end
