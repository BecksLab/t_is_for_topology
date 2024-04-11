"""
adbm_parameters
"""
function adbm_parameters(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary},
                        bodymass::Vector{Float64}; 
                        e::Float64 = 1.0, 
                        a_adbm::Float64 = 0.0189, 
                        ai::Float64 = -0.491, 
                        aj::Float64 = -0.465, 
                        b::Float64 = 0.401, 
                        h_adbm::Float64 = 1.0, 
                        hi::Float64 = 1.0, 
                        hj::Float64 = 1.0, 
                        n::Float64 = 1.0, 
                        ni::Float64 = 0.75, 
                        Hmethod::Symbol = :ratio, 
                        Nmethod::Symbol = :original)

                        
    A = zeros(Int64, (richness(N),richness(N)))
    for i in axes(A, 1)
      for j in axes(A, 2)
        if N.edges[i,j] == true
          A[i,j] = 1
        end
      end
    end

    parameters = Dict{Symbol,Any}(
      :e          => e,
      :a_adbm     => a_adbm,
      :ai         => ai,
      :aj         => aj,
      :b          => b,
      :h_adbm     => h_adbm,
      :hi         => hi,
      :hj         => hj,
      :n          => n,
      :ni         => ni,
      :A          => A
      )
    
    #check Hmethod
    if Hmethod ∈ [:ratio, :power]
      parameters[:Hmethod] = Hmethod
    else
      error("Invalid value for Hmethod -- must be :ratio or :power")
    end
    # check Nmethod
    if Nmethod ∈ [:original, :biomass]
      parameters[:Nmethod] = Nmethod
    else
      error("Invalid value for Nmethod -- must be :original or :biomass")
    end
    # add empty cost matrix
    S = size(parameters[:A],2)
    parameters[:costMat] = ones(Float64,(S,S))

    # Identify producers
    is_producer = vec(sum(A, dims = 2) .== 0)
    parameters[:is_producer] = is_producer 

    # add bodymass
    parameters[:bodymass] = bodymass

    return parameters
end


"""
_get_adbm_terms

This function takes the parameters for the ADBM model and returns
the final terms used to determine feeding patterns. It is used internally by  adbmmodel().
"""
function _get_adbm_terms(S::Int64, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})
  E = parameters[:e] .* parameters[:bodymass]
  if parameters[:Nmethod] == :original
    N = parameters[:n] .* (parameters[:bodymass] .^ parameters[:ni])
  elseif parameters[:Nmethod] == :biomass
    N = biomass
  end
  A_adbm = parameters[:a_adbm] * (parameters[:bodymass].^parameters[:aj]) * (parameters[:bodymass].^parameters[:ai])' # a * pred * prey
  for i = 1:S #for each prey
    A_adbm[:,i] = A_adbm[:,i] .* N[i]
  end
  λ = A_adbm
  if parameters[:Hmethod] == :ratio
    H = zeros(Float64,(S,S))
    ratios = (parameters[:bodymass] ./ parameters[:bodymass]')' #PREDS IN ROWS : PREY IN COLS
    for i = 1:S , j = 1:S
      if ratios[j,i] < parameters[:b]
      H[j,i] =  parameters[:h_adbm] / (parameters[:b] - ratios[j,i])
      else
      H[j,i] = Inf
      end
    end
  elseif parameters[:Hmethod] == :power
    H = parameters[:h_adbm] * (parameters[:bodymass].^parameters[:hj]) * (parameters[:bodymass].^parameters[:hi])' # h * pred * prey
  end

  adbmTerms = Dict{Symbol,Any}(
  :E => E,
  :λ => λ,
  :H => H)

  return(adbmTerms)
end

"""
_get_feeding_links

This function takes the terms calculated by _get_adbm_terms() and uses them to determine the feeding
links of species j. Used internally by adbmmodel().
"""
function _get_feeding_links(S::Int64,E::Vector{Float64}, λ::Array{Float64},
   H::Array{Float64},biomass::Vector{Float64},j)

  profit = E ./ H[j,:]
  # Setting profit of species with zero biomass  to -1.0
  # This prevents them being included in the profitSort
  profit[vec(biomass .== 0.0)] .= -1.0

  profs = sortperm(profit,rev = true)

  λSort = λ[j,profs]
  HSort = H[j,profs]
  ESort = E[profs]

  λH = cumsum(λSort .* HSort)
  Eλ = cumsum(ESort .* λSort)

  λH[isnan.(λH)] .= Inf
  Eλ[isnan.(Eλ)] .= Inf

  cumulativeProfit = Eλ ./ (1 .+ λH)

  if all(0 .== cumulativeProfit)
  feeding = []
  else
  feeding = profs[1:maximum(findall(cumulativeProfit .== maximum(cumulativeProfit)))]
  end

  #cumulativeProfit[end] = NaN
  #feeding = profs[(append!([true],cumulativeProfit[1:end-1] .< profitSort[2:end]))]
  return(feeding)
end


"""
adbmmodel

This function returns the food web based on the ADBM model of Petchey et al. 2008.
The function takes the paramteres created by rewire_parameters() and uses 
getADBM_Terms() and getFeedingLinks() to detemine the web structure. This function
is called using the callback to include rewiring into biomass simulations.

#### References

Petchey, Owen L., Andrew P. Beckerman, Jens O. Riede, and Philip H. Warren.
2008. “Size, Foraging, and Food Web Structure.” Proceedings of the National
Academy of Sciences 105 (11): 4191–96. https://doi.org/10.1073/pnas.0710672105.

"""
function adbmmodel(S::Int64, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})
  adbmMAT = zeros(Bool,(S,S))
  adbmTerms = _get_adbm_terms(S,parameters,biomass)
  E = adbmTerms[:E]
  λ = adbmTerms[:λ]
  H = adbmTerms[:H]
  for j = 1:S
    if !parameters[:is_producer][j]
      if biomass[j] > 0.0
        feeding = _get_feeding_links(S,E,λ,H,biomass,j)
        adbmMAT[j,feeding] .= 1
      end
    end
  end
  edges = Binary(adbmMAT)
  nodes = Unipartite(edges)
  return SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
end