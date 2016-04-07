
using CTDE

function Initiate(time, state, who)
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
            return (true, [params["lambda"]])
        else
            return (false, [])
        end
    end
end


function Terminate(time, state, who)
    state[who]=0
    [who]
end


function BuildTerminateRate(params)
    function TerminateRate(time, state, who)
        (state[who]==1, [params["beta"]])
    end
end


function Elongate(time, state, start, finish)
    state[start]=0
    state[finish]=0
    [start, finish]
end

function BuildElongateRate(params)
    function ElongateRate(time, state, who, next)
        (state[who]==1 && state[next]==0, [params["k"], params["lambda"]])
    end
end


function BuildElongateEndRate(params)
    function ElongateEndRate(time, state, who)
        (state[who]==1, [params["k"], params["lambda"]])
    end
end


function Decay(time, state, who)
    state[who]=0
end

function BuildDecayRate(params)
    function DecayRate(time, state, who, next)
        (state[who]==1, [params["gamma"]])
    end
end



function MakeProcess(parameters, rng)
    # The +1 is because we put into the last place whether the mRNA has decayed.
    L=parameters[:L]
    state=zeros(Int, L+1)
    state[length(state)]=1 # The mRNA starts whole.
    skip=parameters[:skip]

    process=PartialProcess(state)

    initiate_rate=MemoryIntensity(BuildTerminateRate(parameters),
            TransitionExponential(1.0))
    AddTransition!(process,
        initiate_rate, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        Initiate, [1],
        "initiate")

    terminate_rate=MemoryIntensity(BuildTerminateRate(parameters),
            TransitionExponential(1.0))
    AddTransition!(process,
        terminate, [L],
        Terminate, [L],
        "terminate")

    decay_rate=MemoryIntensity(BuildDecayRate(params),
            TransitionExponential(1.0))
    AddTransition!(process,
        decay_rate, [length(state)],
        Decay, [length(state)],
        "decay")

    for codon_idx = 1:(L-skip)
        elongate_rate=MemoryIntensity(BuildElongateRate(params),
                TransitionWeibull(1.0, 2.0))
        AddTransition!(process,
            elongate_rate, [codon_idx, codon_idx+skip],
            Elongate, [codon_idx, codon_idx+1],
            "elongate$(codon_idx)")
    end

    # The last few codons can't have an mRNA in front of them.
    for codon_idx = (L-skip):(L-1)
        elongate_rate=MemoryIntensity(BuildElongateRate(params),
                TransitionWeibull(1.0, 2.0))
        AddTransition!(process,
            elongate_rate, [codon_idx],
            Elongate, [codon_idx, codon_idx+1],
            "elongate$(codon_idx)")
    end

    (process, state)
end

function Observe(state::Array{Int,1}, affected, clock_name, time::Float64)
    state[length(state)]==1
end


function Run()
    rng=MersenneTwister(333333)
    N=3
    parameters=Dict(
        :alpha => 1/45,
        :beta => 2/15,
        :lambda => 1/1350,
        :L => 815,
        :skip => 12
        )
    process, state=MakeProcess(N, parameters, rng)
    sampler=NextReactionHazards()

    RunSimulation(process, sampler, Observe, rng)

    # MakePlots(observer)
    # ShowFractions(observer)
end

Run()
