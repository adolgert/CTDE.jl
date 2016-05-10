# This is a way to check that samplers are correct.
# It runs a well-mixed Susceptible-Infected-Susceptible
# and measures not the average residence in each state
# but the probability of being in a particular state at a
# particular time, so it looks at the master equation.
using CTDE
using Gadfly

import Base: getindex, setindex!
import CTDE: Update!, Reset!

function Recover(state, who)
    state[who]=0
    [who]
end

function RecoverParameters(time, state, who, others...)
    enabled=(state[who]==1)
    # Forbid recovery if this is the only one infectious.
    found_nonzero=false
    for nz_idx = 1:length(others)
        if state[others[nz_idx]]>0
            found_nonzero=true
            break
        end
    end
    (enabled && found_nonzero, [1.0, 2.0])
end

function Infect(state, who)
    state[who]=1
    [who]
end

function InfectParameters(time, state, who, whom)
    (state[who]==1 && state[whom]==0, [1.0])
end


function MakeProcess(N, parameters, rng)
    # N=10 # number of individuals
    state=zeros(Int, N)
    state[1]=1

    process=PartialProcess(state)
    transition_idx=1

    for midx = 1:N
        hazard=MemorylessIntensity(RecoverParameters,
                TransitionWeibull(1.0, 2.0))
                # TransitionExponential(parameters[:Gamma]))
        depends=[midx]
        for dep_idx=1:N
            if dep_idx!=midx
                push!(depends, dep_idx)
            end
        end
        # We add the index=transition_idx only for Direct methods.
        AddIntegratedTransition!(process,
            hazard, depends,
            Recover, [midx],
            "r$midx", index=transition_idx )
        transition_idx+=1

        for sidx=1:N
            if sidx!=midx
                infect=MemorylessIntensity(InfectParameters,
                        TransitionExponential(parameters[:Beta]))
                        #TransitionWeibull(parameters[:Beta], 1.5, 0))
                AddIntegratedTransition!(process,
                    infect, [midx, sidx],
                    Infect, [sidx],
                    "i$midx$sidx", index=transition_idx)
                transition_idx+=1
            end
        end
    end

    (process, state)
end

type StateHistory
    enter::Array{Float64,1}
    leave::Array{Float64,1}
    StateHistory()=new(Array(Float64,0), Array(Float64,0))
end

function AsXY(s::StateHistory, normalization)
    assert(length(s.enter)>0 && length(s.leave)>0)
    sort!(s.enter)
    sort!(s.leave)

    enters=0
    leaves=0
    x=zeros(Float64, length(s.enter)+length(s.leave))
    y=zeros(Float64, length(s.enter)+length(s.leave))
    running=0
    while enters<length(s.enter) || leaves<length(s.leave)
        when=0.0
        if length(s.enter)-enters==0
            leaves+=1
            when=s.leave[leaves]
        elseif length(s.leave)-leaves==0
            enters+=1
            when=s.enter[enters]
        else
            if s.enter[enters+1]<s.leave[leaves+1]
                enters+=1
                when=s.enter[enters]
            else
                leaves+=1
                when=s.leave[leaves]
            end
        end
        x[enters+leaves]=when
        y[enters+leaves]=(enters-leaves)
    end
    y/=normalization
    (x, y)
end

type SamplingObserver
    measure_state::Array{StateHistory, 1}
    measure_time
    measure_summary
    measure_idx
    measure_cnt
    regeneration_time
    regeneration_cnt
    SamplingObserver(state_cnt, measure_cnt)=new(
        [StateHistory() for i in 1:state_cnt],
        0.0, 1, 0, measure_cnt, 0.0, 0)
end

function ObserveState(so::SamplingObserver, state::Array{Int,1},
        affected, clock_name, time::Float64)
    start=so.measure_time-so.regeneration_time
    finish=time-so.regeneration_time
    push!(so.measure_state[so.measure_summary].enter, start)
    push!(so.measure_state[so.measure_summary].leave, finish)
    so.measure_time=time
    so.measure_idx+=1
    so.measure_summary=countnz(state)
    if so.measure_summary==1
        so.regeneration_time=so.measure_time
        so.regeneration_cnt+=1
    end
    so.measure_idx < so.measure_cnt
end

function Observer(so::SamplingObserver)
    function sobserve(state::Array{Int,1}, affected, clock_name, time::Float64)
        ObserveState(so, state, affected, clock_name, time)
    end
end

function MakePlots(so::SamplingObserver)
    for plot_idx in 1:length(so.measure_state)
        fraction=length(so.measure_state[plot_idx].enter)/so.regeneration_time
        print("Plot $plot_idx fraction $fraction\n")
        x1, y1=AsXY(so.measure_state[plot_idx], so.regeneration_cnt)
        levelplot=plot(x=x1, y=y1, Geom.line)
        draw(PDF("level$(plot_idx).pdf", 4inch, 3inch), levelplot)
        # for print_idx = 1:length(x1)
        #     print(x1[print_idx], "\t", y1[print_idx], "\n")
        # end
    end
end

function ShowFractions(so::SamplingObserver)
    for plot_idx in 1:length(so.measure_state)
        fraction=length(so.measure_state[plot_idx].enter)/so.regeneration_time
        print("Plot $plot_idx fraction $fraction\n")
    end
end

function Run()
    rng=MersenneTwister(333333)
    N=3
    parameters=Dict(:Gamma =>1.0, :Beta => 1.0)
    process, state=MakeProcess(N, parameters, rng)
    observer=SamplingObserver(N, 1000000)
    sampler=FixedDirect(TransitionCount(process))
    # sampler=NaiveSampler()
    # sampler=NextReactionHazards()
    # sampler=FirstReaction()

    RunSimulation(process, sampler, Observer(observer), rng)

    MakePlots(observer)
    # ShowFractions(observer)
end

Run()
