# Modified from Poisot, Timothée. 2023. “Guidelines for the Prediction of
# Species Interactions Through Binary Classification.” Methods in Ecology and
# Evolution 14 (5): 1333–45. https://doi.org/10.1111/2041-210X.14071.

include("internals.jl")

struct ConfusionMatrix
    tp::Any
    tn::Any
    fp::Any
    fn::Any
    function ConfusionMatrix(a, b, c, d)
        s = a + b + c + d
        return new(a / s, b / s, c / s, d / s)
    end
end

tpr(M::ConfusionMatrix) = M.tp / (M.tp + M.fn)
tnr(M::ConfusionMatrix) = M.tn / (M.tn + M.fp)
ppv(M::ConfusionMatrix) = M.tp / (M.tp + M.fp)
npv(M::ConfusionMatrix) = M.tn / (M.tn + M.fn)
fnr(M::ConfusionMatrix) = M.fn / (M.fn + M.tp)
fpr(M::ConfusionMatrix) = M.fp / (M.fp + M.tn)
informedness(M::ConfusionMatrix) = tpr(M) + tnr(M) - 1.0
mcc(M::ConfusionMatrix) =
    (M.tp * M.tn - M.fp * M.fn) /
    sqrt((M.tp + M.fp) * (M.tp + M.fn) * (M.tn + M.fp) * (M.tn + M.fn))

function ∫(x::Array{T}, y::Array{T}) where {T<:Number}
    S = zero(Float64)
    for i = 2:length(x)
        S += (x[i] - x[i-1]) * (y[i] + y[i-1]) * 0.5
    end
    return .-S
end

function benchmark(obs, pred; levels = 500)
    thresholds = LinRange(minimum(pred), maximum(pred), levels)
    obs = Bool.(obs) # make Bool
    M = Vector{ConfusionMatrix}(undef, length(thresholds))
    for (i, τ) in enumerate(thresholds)
        binpred = pred .>= τ
        tp = sum(obs .& binpred)
        tn = sum(.!obs .& .!binpred)
        fp = sum(.!obs .& binpred)
        fn = sum(obs .& .!binpred)
        M[i] = ConfusionMatrix(tp, tn, fp, fn)
    end
    AUPRC = ∫(tpr.(M), ppv.(M))
    _mcc = mcc.(M)
    𝑚𝑐𝑐 = findmax(_mcc[.!isnan.(_mcc)])[1] # need to avoid NaN bias
    τ = thresholds[last(findmax(informedness.(M)))]
    return AUPRC, 𝑚𝑐𝑐, τ
end
