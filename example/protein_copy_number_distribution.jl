# This is an implementation of a model from Gorissen
# and Vanderzande [1].
#
# [1] M. Gorissen and C. Vanderzande, “Ribosome Dwell Times and
# the Protein Copy Number Distribution,” J Stat Phys, vol. 148,
# pp. 628–636, 2012.

using CTDE
using UNURAN
import Base: start, next, done
using Logging
@Logging.configure(level=DEBUG)

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
        if state[who]==0
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
    state[finish]=1
    [start, finish]
end


function BuildElongateRate(params)
    function ElongateRate(time, state, who, next)
        # The k and n used in the article are a shape and a rate,
        # not a shape and a scale.
        (state[who]==1 && state[next]==0, [params[:n], 1/params[:k]])
    end
end


function BuildElongateEndRate(params)
    function ElongateEndRate(time, state, who)
        (state[who]==1, [params[:n], 1/params[:k]])
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



function MakeProcess(parameters, rng, decay=true)
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
            TransitionExponential(rng, 1.0))
    AddIntegratedTransition!(process,
        initiate_rate, Array{Int,1}(1:skip),
        Initiate, [1],
        "initiate")

    terminate_rate=MemoryIntensity(BuildTerminateRate(parameters),
            TransitionExponential(rng, 1.0))
    AddIntegratedTransition!(process,
        terminate_rate, [L],
        Terminate, [L],
        "terminate")

    if decay
        decay_rate=MemoryIntensity(BuildDecayRate(parameters),
                TransitionExponential(rng, 1.0))
        AddIntegratedTransition!(process,
            decay_rate, [DecayPlace],
            Decay, [DecayPlace],
            "decay")
    end

    for codon_idx = 1:(L-skip)
        elongate_rate=MemoryIntensity(BuildElongateRate(parameters),
                CTDE.GammaUnur(rng, 1.0, 2.0))
        AddIntegratedTransition!(process,
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
                GammaUnur(rng, 1.0, 2.0))
        AddIntegratedTransition!(process,
            elongate_rate, [codon_idx],
            Elongate, [codon_idx, codon_idx+1],
            "elongate$(codon_idx)")
    end

    (process, state)
end


type Observations
    row::Array{Int,1}
    leave::Array{Float64,1}
    Observations()=new(Array{Int,1}(), Array{Float64,1}())
end

function Write(obs::Observations, name)
    writedlm("$(name)row.txt", obs.row)
    writedlm("$(name)leave.txt", obs.leave)
end

function Clear!(obs::Observations)
    obs.row=Array{Int,1}()
    obs.leave=Array{Float64,1}()
end

function FromFile(name::AbstractString)
    obs=Observations()
    obs.row=readdlm("$(name)row.txt")
    obs.leave=readdlm("$(name)leave.txt")
end

function FromFiles(name::AbstractString)
    obs=Observations()
    all_files=readdir()
    for f in all_files
        if startswith(f, name) && endswith(f, ".txt")
            row=readdlm(f)
            leave=readdlm(f)
            row .+= length(obs.row)
            obs.row=[obs.row ; row]
            obs.leave=[obs.leave ; leave]
        end
    end
    obs
end


function Observer(store::Observations)
    function Observe(state::Array{Int,1}, affected, clock_name, time::Float64)
        println("observe clock $clock_name $time")
        if clock_name=="terminate"
            # print("TERMINATE $time\n")
            push!(store.leave, time)
        else
            0 # ignore other stuff
        end
        state[length(state)]==1
    end
end

function TimedObserver(store::Observations)
    function Observe(state::Array{Int,1}, affected, clock_name, time::Float64)
        if clock_name=="terminate"
            push!(store.leave, time)
        else
            0 # ignore other stuff
        end
        time<50000
    end
end

# Make an iterator over the exit times for each run.
# It returns a list of [time, time, time] for exits.
function start(store::Observations)
    1
end

function next(obs::Observations, row_idx::Int)
    row_begin=obs.row[row_idx]
    if row_idx<length(obs.row)
        row_end=obs.row[row_idx+1] - 1
    else
        row_end=length(obs.row)
    end
    (obs.leave[row_begin:row_end], row_idx+1)
end

function done(store::Observations, row_idx::Int)
    row_idx>length(store.row)
end


"""
What is the distribution of how many proteins it produces
before decay?
"""
function DistributionOfProteinsProduced(obs::Observations)
    histogram=Dict{Int,Int}()
    for exits in obs
        exit_cnt=length(exits)
        if haskey(histogram, exit_cnt)
            histogram[exit_cnt]+=1
        else
            histogram[exit_cnt]=1
        end
    end

    m=maximum([k for k::Int in keys(histogram)])
    last_nonzero=-1
    total=0
    for i=0:m
        if haskey(histogram, i)
            total+=histogram[i]
            print("$i\t$(histogram[i])\n")
        else
            if last_nonzero<0
                last_nonzero=i-1
            end
        end
    end
    reduced_histogram=zeros(Float64, last_nonzero+1)
    for nz=0:last_nonzero
        reduced_histogram[nz+1]=histogram[nz]/total
    end
    histogram, reduced_histogram
end


"""
How much time does it take the second mRNA to leave after
the first one? So on up to 6. For the first one, with what distribution
does it leave?
"""
function InterDepartureWaitingTime(obs::Observations)
    histogram=Array{EmpiricalDistribution,1}()
    for gen_idx = 1:7
        push!(histogram, EmpiricalDistribution())
    end

    for exits in obs
        if length(exits)>0
            push!(histogram[1], exits[1])
        end
        for l_idx = 2:min(length(exits), 7)
            push!(histogram[l_idx], exits[l_idx]-exits[l_idx-1])
        end
    end
    histogram
end


"""
Assuming decay doesn't happen, how quickly do mRNA leave,
depending on the time since start? This constructs an
empirical distribution on which to perform an estimate of
the hazard rate for leaving.
"""
function LeavingTimeDistribution(obs::Observations)
    ed=EmpiricalDistribution()
    cnt=0
    for exits in obs
        for e in exits
            push!(ed, e)
        end
        cnt+=1
    end
    (ed, cnt)
end


function Reset!(state)
    state[1:length(state)]=0
    state[length(state)]=1
end


function Run()
    #rng=MersenneTwister(333333)
    seed=UInt64[111, 222, 333, 444, 555, 666]
    RngStream_SetPackageSeed(seed)
    rng=unur_urng_rngstream_new("urng-1")

    run_cnt=50
    parameters=Dict(
        :n => 0.5,
        :k => 0.5,
        :alpha => 1/45,
        :beta => 2/15,
        :lambda => 1/1350,
        :L => 815,
        :l => 12
        )
    process, state=MakeProcess(parameters, rng)
    obs=Observations()
    for run_idx = 1:run_cnt
        Reset!(state)
        sampler=NextReactionHazards()

        push!(obs.row, length(obs.leave)+1)
        RunSimulation(process, sampler, Observer(obs), rng)
    end

    DistributionOfProteinsProduced(obs)
    InterDepartureWaitingTime(obs)
    LeavingTimeDistribution(obs)
end
