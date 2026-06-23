Zygote.@adjoint function ScatteredOut(output, k=1)
    function ∇scattered(δ)
        return (δ.output, nothing)
    end
    ScatteredOut(output, k), ∇scattered
end
Zygote.@adjoint function ScatteredOut(m, k, output)
    function ∇scattered(δ)
        return (nothing, nothing, δ.output)
    end
    ScatteredOut{eltype(output),length(output)}(m, k, output), ∇scattered
end

Zygote.@adjoint function getindex(F::T, i::Integer) where {T<:Scattered}
    function getInd_rrule(Ȳ)
        # zeroNonRefed = map(ii -> ii - 1 == i ? Ȳ : zeros(eltype(F.output[ii]), size(F.output[ii])...), (1:length(F.output)...,))
        zeroNonRefed = map(ii -> ii - 1 == i ? Ȳ : zero(F.output[ii]), (1:length(F.output)...,))
        ∂F = T(F.m, F.k, zeroNonRefed)
        return ∂F, nothing
    end
    return getindex(F, i), getInd_rrule
end

Zygote.@adjoint function getindex(F::T, inds::AbstractArray) where {T<:Scattered}
    function getInd_rrule(Ȳ)
        # zeroNonRefed = map(ii -> ii - 1 in inds ? Ȳ[indexin(ii - 1, inds)[1]] : zeros(eltype(F.output[ii]), size(F.output[ii])...), (1:length(F.output)...,))
        zeroNonRefed = map(ii -> ii - 1 in inds ? Ȳ[indexin(ii - 1, inds)[1]] : zero(F.output[ii]), (1:length(F.output)...,))
        ∂F = T(F.m, F.k, zeroNonRefed)
        return ∂F, nothing
    end
    return getindex(F, inds), getInd_rrule
end

Zygote.@adjoint function getindex(x::T, p::pathLocs) where {T<:Scattered}
    function getInd_rrule(Δ)
        # zeroNonRefed = map(ii -> zeros(eltype(x.output[ii]), size(x.output[ii])...), (1:length(x.output)...,))
        zeroNonRefed = map(ii -> zero(x.output[ii]), (1:length(x.output)...,))
        ∂x = T(x.m, x.k, zeroNonRefed)
        ∂x[p] = Δ
        return ∂x, nothing
    end
    return getindex(x, p), getInd_rrule
end

function rrule(::typeof(flatten), scatRes)
    function ∇flatten(Δarray)
        return (NoTangent(), roll(Δarray, scatRes),)
    end
    return flatten(scatRes), ∇flatten
end
function rrule(::typeof(roll), toRoll, stOutput)
    function ∇roll(Δ)
        return NoTangent(), flatten(Δ), NoTangent()
    end
    return roll(toRoll, stOutput), ∇roll
end


function ChainRulesCore.rrule(::typeof(normalize), x, Nd)
    n = ndims(x)
    totalThisLayer = prod(size(x)[(Nd+1):(n-1)])
    sumSqDims = 1:(n-1)
    normSq = sum(abs.(x) .^ 2, dims=sumSqDims)
    scale = totalThisLayer ./ sqrt.(normSq)
    scale = ifelse.(isnan.(scale) .| isinf.(scale), one(eltype(scale)), scale)
    y = x .* scale

    function normalize_pullback(Δy)
        # All operations stay on the same device as x
        ∂x = Δy .* scale
        return NoTangent(), ∂x, NoTangent()
    end
    return y, normalize_pullback
end

Zygote.@adjoint function normalize(x, Nd)
    n = ndims(x)
    totalThisLayer = prod(size(x)[(Nd+1):(n-1)])
    sumSqDims = 1:(n-1)
    normSq = sum(abs.(x) .^ 2, dims=sumSqDims)
    scale = totalThisLayer ./ sqrt.(normSq)
    scale = ifelse.(isnan.(scale) .| isinf.(scale), one(eltype(scale)), scale)
    y = x .* scale
    function normalize_pullback(Δy)
        scale_adapted = adapt(typeof(Δy), scale)
        ∂x = Δy .* scale_adapted
        return ∂x, nothing
    end
    return y, normalize_pullback
end