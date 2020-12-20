import LightGraphs: SimpleGraph, add_vertices!, neighbors
import Distributions: Uniform, Exponential, Weibull, Gamma, ContinuousUnivariateDistribution, params
import Random: MersenneTwister, rand
import Base: iterate, eltype
using Traceur
using Profile
using DataFrames
using BenchmarkTools


const CTDist = ContinuousUnivariateDistribution

### Version 1: Visitor pattern on an array

function visit(process, when::Float64, rng, Key::Type{T}) where {T}
    total::Float64 = 0.0
    cumulative = zeros(Float64, 0)
    keys = Array{T,1}()
    hazards(process) do clock_key::T, distribution::Exponential
        total += params(distribution)[1]::Float64
        push!(cumulative, total)
        push!(keys, clock_key)
    end

    if total > eps(Float64)
        chosen = searchsortedfirst(cumulative, rand(rng, Uniform(0, total)))
        @assert chosen < length(cumulative) + 1
        return (when - log(rand(rng)::Float64) / total, keys[chosen])
    else
        return (Inf, zero(T))
    end
end


function visit(process, when, rng)
    visit(process, when, rng, keytype(process))
end


function hazards(visitor::Function, arr::Array{T, 1}) where T <: CTDist
    for (i, d) in enumerate(arr)
        visitor(i, d)
    end
end


keytype(::Array) = Int

const same_seed = 9237842
const cnt = 10
rng = MersenneTwister(same_seed)
distribution_array = CTDist[Exponential(1) for i in 1:cnt]
visit(distribution_array, 0.0, rng, Int64)
@code_warntype visit(distribution_array, 0.0, rng, Int64)

rng = MersenneTwister(same_seed)
exponential_array = [Exponential(1) for i in 1:10]
visit(exponential_array, 0.0, rng)


### Version 2: Iterator on an array


function iterate_over(hazards, when, rng, Key::Type{T}) where T
    total = 0.0
    cumulative = zeros(Float64, 0)
    keys = Array{T, 1}()
    for (clock_key::T, distribution::Exponential{Float64}) in hazards
        total += params(distribution)[1]::Float64
        push!(cumulative, total)
        push!(keys, clock_key)
    end

    if total > eps(Float64)
        chosen = searchsortedfirst(cumulative, rand(rng, Uniform(0, total)))
        @assert chosen < length(cumulative) + 1
        return (when - log(rand(rng)) / total, keys[chosen])
    else
        return (Inf, zero(T))
    end
end


function iterate_over(hazards, when, rng)
    iterate_over(hazards, when, rng, eltype(hazards).parameters[1])
end

exponential_array = fill(Exponential(), 10)
iterate_over(enumerate(exponential_array), 0.0, rng)
@code_warntype iterate_over(enumerate(exponential_array), 0.0, rng, Int64)

### Version 3: Visit a SimpleGraph

struct GraphTransitions{T <: CTDist}
    g::SimpleGraph{Int64}
    distributions::Array{T, 1}
end


mutable struct GraphProcess{T <: CTDist}
    transitions::GraphTransitions{T}
    fired_transition::Int64
end


function hazards(visitor::Function, process::GraphProcess{T}) where T <: CTDist
    for n in neighbors(process.transitions.g, process.fired_transition)
        visitor(n, process.transitions.distributions[n])
    end
end

keytype(::GraphProcess) = Int

const graph_seed = 903472
vertex_cnt = 10
edge_cnt = 20
g = SimpleGraph{Int64}(vertex_cnt, edge_cnt; seed = graph_seed)
distribution_array = CTDist[Exponential(1) for i in 1:vertex_cnt]
process = GraphProcess(GraphTransitions(g, distribution_array), 1)
visit(process, 0.0, rng)
@code_warntype visit(process, 0.0, rng)


### Version 4: Iterate over a SimpleGraph

function iterate(process::GraphProcess{T}) where T <: CTDist
    inner_iterator = neighbors(process.transitions.g, process.fired_transition)
    iter_pair = iterate(inner_iterator)
    if !isnothing(iter_pair)
        i, inner_state = iter_pair
        ((i, process.transitions.distributions[i]), (inner_iterator, inner_state))
    else
        nothing
    end
end


function iterate(process::GraphProcess{T}, state) where T <: CTDist
    inner_iterator, inner_state = state
    iter_pair = iterate(inner_iterator, inner_state)
    if !isnothing(iter_pair)
        i, state = iter_pair
        ((i, process.transitions.distributions[i]), (inner_iterator, state))
    else
        nothing
    end
end

eltype(::GraphProcess{T}) where T = Tuple{Int64, T}
iterate(process)
@code_warntype iterate(process)

iterate_over(process, 0.0, rng)
@code_warntype iterate_over(process, 0.0, rng, Int64)


### Now compare

# We have 8 cases and want to compare visitor pattern to iterator pattern,
# so that's 4 comparisons.
#
# 1. Array
#    a. All exponentials
#    b. Continuous-time distributions
# 2. Graph
#    a. All exponentials
#    b. Continuous-time distributions


function visit_array(DistType, seed, dist_cnt, draw_cnt)
    rng = MersenneTwister(seed)
    distribution_array = DistType[Exponential(1) for i in 1:dist_cnt]
    when = 0.0
    for draw_idx in draw_cnt
        when, which = visit(distribution_array, when, rng)
    end
end

function iterate_array(DistType, seed, dist_cnt, draw_cnt)
    rng = MersenneTwister(seed)
    distribution_array = DistType[Exponential(1) for i in 1:dist_cnt]
    when = 0.0
    for draw_idx in draw_cnt
        when, which = iterate_over(enumerate(exponential_array), when, rng)
    end
end


function visit_graph(DistType, seed, dist_cnt, edge_cnt, draw_cnt)
    rng = MersenneTwister(seed)
    g = SimpleGraph{Int64}(dist_cnt, edge_cnt; seed = graph_seed)
    distribution_array = DistType[Exponential(1) for i in 1:dist_cnt]
    process = GraphProcess(GraphTransitions(g, distribution_array), 1)

    when = 0.0
    for draw_idx in draw_cnt
        when, which = visit(process, when, rng)
        process.fired_transition = ifelse(isfinite(which), which, 1)
    end
end



function iterate_graph(DistType, seed, dist_cnt, edge_cnt, draw_cnt)
    rng = MersenneTwister(seed)
    g = SimpleGraph{Int64}(dist_cnt, edge_cnt; seed = graph_seed)
    distribution_array = DistType[Exponential(1) for i in 1:dist_cnt]
    process = GraphProcess(GraphTransitions(g, distribution_array), 1)

    when = 0.0
    for draw_idx in draw_cnt
        when, which = iterate_over(process, when, rng)
        process.fired_transition = ifelse(isfinite(which), which, 1)
    end
end


const same_seed = 9237842
const dist_cnt = 1000
const draw_cnt = 1000
# 68 μs, 84 kb
res = @benchmark visit_array(Exponential{Float64}, same_seed, dist_cnt, draw_cnt)
median(res).time
# 98 μs, 138 kb
@benchmark visit_array(CTDist, same_seed, dist_cnt, draw_cnt)
# 49 μs, 36 kb
@benchmark iterate_array(Exponential{Float64}, same_seed, dist_cnt, draw_cnt)
# 52 μs, 52 kb
@benchmark iterate_array(CTDist, same_seed, dist_cnt, draw_cnt)

const dist_cnt = 1000
const edge_cnt = 2000
const draw_cnt = 1000
# 91 μs, 143 kb
@benchmark visit_graph(Exponential{Float64}, same_seed, dist_cnt, edge_cnt, draw_cnt)
# 95 μs, 158 kb
@benchmark visit_graph(CTDist, same_seed, dist_cnt, edge_cnt, draw_cnt)
# 90 μs, 143 kb
@benchmark iterate_graph(Exponential{Float64}, same_seed, dist_cnt, edge_cnt, draw_cnt)
# 96 μs, 158 kb
@benchmark iterate_graph(CTDist, same_seed, dist_cnt, edge_cnt, draw_cnt)
