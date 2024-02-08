    
"""
    adbm

function to generate netwrok using ADBM, Petchy and friends 2008
"""
function adbm(
    adbmmatrix::Matrix{Float64},
    mass::Vector{UInt8},
    energy_content::Vector{UInt8};
    resource_constant::Float64 = -0.75,
    energy_constant::Float64 = 1.0,
    attack_constant::Float64 = 1.0,
    handling_constant::Float64 = 1.0,
    b::Float64 = 1.0,
)
    a = attack_constant
    for j in axes(adbmmatrix, 1)
        for i in axes(adbmmatrix, 2)
            if mass[i]/mass[j] ≥ b #immediately set to zero
                adbmmatrix[j,i] = 0.0  
            end
            _handling_time = handling_constant/(b - (mass[i]/mass[j]))
            _resource_density = resource_density[i]*resource_constant
            _energy_content = energy_content[i]*energy_constant
            _attack_rate = a*(mass[i]^a)*(mass[j]^a)
            _encounter_rate = _resource_density*_attack_rate
            
            adbmmatrix[j,i] = (_encounter_rate*_energy_content)/(_resource_density*_handling_time)

    end
    return adbmmatrix
end