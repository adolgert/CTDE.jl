
using CTDE

function Initiate(state, who)
    state[who]=1
    [who]
end


"""
This function builds an initiation function so that it can
enclose the rate parameter.
"""
function BuildInitiateRate(params)
    function IntiateRate(time, state, who, neighbors...)
        if state[who]==1
            for empty_idx = 1:length(neighbors)
                if state[empty_idx]==1
                    return (false, [])
                end
            end
            return (true, [params[:alpha]])
        else
            return (false, [])
        end
    end
end


function Terminate(state, who)
    state[who]=0
    [who]
end


function BuildTerminateRate(params)
    function TerminateRate(time, state, who)
        (state[who]==1, [params[:beta]])
    end
end


function Elongate(state, start, finish)
    state[start]=0
    state[finish]=0
    [start, finish]
end


function BuildElongateRate(params)
    function ElongateRate(time, state, who, next)
        # The k and n used in the article are a shape and a rate,
        # not a shape and a scale.
        (state[who]==1 && state[next]==0, [params[:n], params[:k]])
    end
end


function BuildElongateEndRate(params)
    function ElongateEndRate(time, state, who)
        (state[who]==1, [params[:n], params[:k]])
    end
end


function Decay(state, who)
    state[who]=0
    [who]
end

function BuildDecayRate(params)
    function DecayRate(time, state, who)
        (state[who]==1, [params[:lambda]])
    end
end



function MakeProcess(parameters, rng)
    # The +1 is because we put into the last place whether the mRNA has decayed.
    L=parameters[:L]
    DecayPlace=L+1
    state=zeros(Int, DecayPlace)
    state[length(state)]=1 # The mRNA starts whole.
    # Skip is the width of the mRNA. An mRNA at 1 can't move
    # if there is another mRNA at 1+skip.
    skip=parameters[:l]

    process=PartialProcess(state)

    initiate_rate=MemoryIntensity(BuildInitiateRate(parameters),
            TransitionExponential(1.0))
    AddTransition!(process,
        initiate_rate, Array{Int,1}(1:skip),
        Initiate, [1],
        "initiate")

    terminate_rate=MemoryIntensity(BuildTerminateRate(parameters),
            TransitionExponential(1.0))
    AddTransition!(process,
        terminate_rate, [L],
        Terminate, [L],
        "terminate")

    decay_rate=MemoryIntensity(BuildDecayRate(parameters),
            TransitionExponential(1.0))
    AddTransition!(process,
        decay_rate, [DecayPlace],
        Decay, [DecayPlace],
        "decay")

    for codon_idx = 1:(L-skip)
        elongate_rate=MemoryIntensity(BuildElongateRate(parameters),
                TransitionWeibull(1.0, 2.0))
        AddTransition!(process,
            elongate_rate, [codon_idx, codon_idx+skip],
            Elongate, [codon_idx, codon_idx+1],
            "elongate$(codon_idx)")
    end

    # The last few codons can't have an mRNA in front of them.
    # For instance, say the mRNA width, skip=3. Label the locations
    # L-4, L-3, L-2, L-1, L, L+1, L+2
    #        X,   X,   X, O,
    # Then L is where the mRNA terminates. L+2 is the end of the mRNA
    # when it terminates. The last mRNA location to have interference
    # is at L-3, because it can be blocked by the mRNA at L.
    for codon_idx = (L-skip+1):(L-1)
        elongate_rate=MemoryIntensity(BuildElongateEndRate(parameters),
                TransitionWeibull(1.0, 2.0))
        AddTransition!(process,
            elongate_rate, [codon_idx],
            Elongate, [codon_idx, codon_idx+1],
            "elongate$(codon_idx)")
    end

    (process, state)
end

type Observations
    enter::Array{Float64,1}
    leave::Array{Float64,1}
    Observations()=new(Array{Float64,1}(), Array{Float64,1}())
end


function Observer(store::Observations)
    function Observe(state::Array{Int,1}, affected, clock_name, time::Float64)
        if clock_name=="initiate"
            push!(store.enter, time)
        elseif clock_name=="terminate"
            push!(store.leave, time)
        else
            0 # ignore other stuff
        end
        state[length(state)]==1
    end
end


function Run()
    rng=MersenneTwister(333333)
    run_cnt=1000
    parameters=Dict(
        :n => 0.5,
        :k => 0.5,
        :alpha => 1/45,
        :beta => 2/15,
        :lambda => 1/1350,
        :L => 815,
        :l => 12
        )
    copy_count=Dict{Int,Int}()
    process, state=MakeProcess(parameters, rng)
    for run_idx = 1:run_cnt
        obs=Observations()
        sampler=NextReactionHazards()

        RunSimulation(process, sampler, Observer(obs), rng)
        if haskey(copy_count, length(obs.leave))
            copy_count[length(obs.leave)]+=1
        else
            copy_count[length(obs.leave)]=1
        end
    end
    print(copy_count)
end

Run()
