using DataStructures

"""
Classic Direct method for exponential transitions
"""
struct MarkovDirect
end


function Next(rm::MarkovDirect, process, rng)
    total = 0.0
    cumulative = zeros(Float64, 0)
    keys = Array{Any,1}()
    Hazards(process, rng) do clock, now, enabled, rng2
        total += Parameters(clock.intensity.distribution)[1]
        push!(cumulative, total)
        push!(keys, clock)
    end

    if total > eps(Float64)
        chosen = searchsortedfirst(cumulative,rand(rng) * total)
        assert(chosen < length(cumulative) + 1)
        return (Time(process) - log(rand(rng)) / total, keys[chosen])
    else
        return (Inf, nothing)
    end
end

Observer(fr::DirectMethod) = (hazard, time, updated, rng) -> nothing
