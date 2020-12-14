using Logging
@Logging.configure(level=DEBUG)
using CTDE

import Base: getindex, setindex!
import CTDE: Update!, Reset!


struct PhysicalState
    v::Array{Int,1}
    N::Int
    PhysicalState(N)=new(zeros(Int, N), N)
end

function Reset(state::PhysicalState)
    state.v=zeros(Int, state.N)
    state[1, 1]=true
end

"""
The accessor functions treat S, I, and R as different
substates, even though it's stored in an array.
The values of the state are true or false.
"""
function getindex(state::PhysicalState, who, disease)
    state.v[who]==disease
end

function setindex!(state::PhysicalState, v, who, disease)
    if v
        state.v[who]=disease
    elseif state.v[who]==disease
        state.v[who]=-1
    else
        # Nothing to do if it already isn't in this state.
    end
end



function Recover(state, infectious, recovered)
    state[infectious[1], recovered[2]]=true
    [infectious, recovered]
end

function Infect(state, susceptible, infectious)
    state[susceptible[1], infectious[2]]=true
    [susceptible, infectious]
end


struct RecoverIntensity <: Intensity
    distribution::TransitionExponential
    enabled::Bool
    RecoverIntensity(dist)=new(dist, false)
end


function Update!(intensity::RecoverIntensity, time, state, who)
    modified=:Undefined
    enabled=state[who[1], who[2]]
    if enabled!=intensity.enabled
        if enabled
            intensity.distribution.enabling_time=time
            modified=:Enabled
        else
            modified=:Disabled
        end
        intensity.enabled=enabled
    else
        modified=:Unmodified
    end
    modified
end

struct InfectIntensity <: Intensity
    distribution::TransitionExponential
    enabled::Bool
    InfectIntensity(dist)=new(dist, false)
end


function Update!(intensity::InfectIntensity, time, state, who, whom)
    modified=:Undefined
    enabled=(state[who[1], who[2]] && state[whom[1], whom[2]])
    if enabled!=intensity.enabled
        if enabled
            intensity.distribution.enabling_time=time
            modified=:Enabled
        else
            modified=:Disabled
        end
        intensity.enabled=enabled
    else
        modified=:Unmodified
    end
    modified
end


function MakeProcess(N, parameters, rng)
    state=PhysicalState(N)
    # individual, disease state (0, 1, 2)
    state[1, 1]=true
    process=PartialProcess(state)

    for midx = 1:N
        hazard=RecoverIntensity(TransitionExponential(parameters['r'], 0))
        AddTransition!(process,
            hazard, [(midx, 1)],
            Recover, [(midx, 1), (midx, 2)],
            "r$midx")

        for sidx=1:N
            if sidx!=midx
                infect=InfectIntensity(
                        TransitionExponential(parameters['i'], 0))
                AddTransition!(process,
                    infect, [(midx, 1), (sidx, 0)],
                    Infect, [(sidx, 0), (sidx, 1)],
                    "i$midx$sidx")
            end
        end
    end
    (process, state)
end


############ The Observer

function SeeNothing(state::PhysicalState, affected, clock_name, time::Float64)
    true
end

############ Runs the simulation
function Run()
    if length(ARGS)>3
        individual_cnt=Int(ARGS[1])
        run_cnt=Int(ARGS[2])
        beta=Float64(ARGS[3])
        gamma=Float64(ARGS[4])
        seed=Int(ARGS[5])
    else
        individual_cnt=10
        run_cnt=100
        beta=2
        gamma=1.0
        seed=333333
    end
    parameters=Dict(
        'i'=>beta/individual_cnt, # rate of infection of neighbor
        'r'=>gamma, # latent to infectious
    )
    rng=MersenneTwister(seed)
    process, state=MakeProcess(individual_cnt, parameters, rng)

    final_size=zeros(Int, individual_cnt+1)
    for run_idx=1:run_cnt
        Reset(state)
        sampler=NextReactionHazards()
        RunSimulation(process, sampler, SeeNothing, rng)
        final_size[individual_cnt-countnz(state.v)+1]+=1
    end

    final_size/=sum(final_size)
    for p_idx=1:(individual_cnt+1)
        print("$(p_idx-1)\t$(final_size[p_idx])\n")
    end
end

Run()

