# A Vector Addition System

"""
A `VectorAdditionSystem` is an older form of simulation that is a lot
like Petri nets. The state is a vector.
The system is a list of transitions. Each transition
is an array of values. Negative numbers mean the transition needs to
take this many tokens from the state, meaning the state at those indices
must be an integer at least that large. Positive numbers mean the
transition places tokens into the state. Unlike chemical simulations,
the rate doesn't depend on the number of combinations of species present.
"""
struct VectorAdditionSystem
    transitions::Array{Int, 2}  # states x transitions
    rates::Array{Float64, 1}  # length is transitions
end


function zero_state(vas::VectorAdditionSystem)
    zeros(Int, size(vas.transitions, 1))
end


function hazards!(visitor::Function, vas::VectorAdditionSystem, state)
    for rate_idx in eachindex(vas.rates)
        if all(vas.transitions[:, rate_idx] .+ state >= 0)
            cnt += 1
            visitor(rate_idx, vas.rates[rate_idx], :Enabled, rng)
        end
    end
end


function fire(vas::VectorAdditionSystem, state, transition_idx)
    state .+= vas.transitions[:, transition_idx]
end


function hazards!(visitor::Function, vas::VectorAdditionSystem, state, transition_idx)
    delta_state = vas.transitions[:, transition_idx]
    for rate_idx in eachindex(vas.rates)
        summed = vas.transitions[:, rate_idx] .+ state
        was_enabled = all(summed .>= 0)
        now_enabled = all(summed .+ delta_state .>= 0)
        if was_enabled && !now_enabled
            visitor(rate_idx, vas.rates[rate_idx], :Disabled, rng)
        elseif !was_enabled && now_enabled
            visitor(rate_idx, vas.rates[rate_idx], :Enabled, rng)
        end
    end
end
