S = 10

parameters = Dict{Symbol, Any}()
parameters[:is_producer] = rand(Bool, S)
parameters[:e] = 1.0
parameters[:b] = 1.0
parameters[:bodymass] = rand(Float64, S)
parameters[:Nmethod] = :biomass
parameters[:a_adbm] = 1.0
parameters[:aj] = 1.0
parameters[:ai] = 1.0
parameters[:Hmethod] = :ratio
parameters[:h_adbm] = 1.0
parameters[:hj] = 1.0
parameters[:n] = 1.0
parameters[:ni] = 1.0

adbmmodel(S, parameters, rand(Float64, S))

nestedhierarchymodel(S, 10)