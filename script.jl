include("semimarkov.jl")

using SemiMarkov

include("samplemodels.jl")
cnt=10
seed=32
if length(ARGS)>0
    cnt=int(ARGS[1])
end
if length(ARGS)>1
    seed=int(ARGS[2])
end

rng=MersenneTwister(seed)
#model=sir_explicit(cnt, [0.1, 0.4])
model=sirs_birth_death(cnt)


type OutputObserver
    observations::Array{(Int,Int,Int,Float64),1}
    cnt::Int
    filename::String
    OutputObserver(name)=new(Array((Int,Int,Int,Float64), 10000), 1, name)
end

function save(obs::OutputObserver)
    if obs.cnt>1
        open(obs.filename, "a") do f
            writedlm(f, obs.observations[1:(obs.cnt-1)])
        end
        obs.cnt=1
    end
end

function observe(obs::OutputObserver, state)
    setindex!(obs.observations, (
        length(state.marking["s"]),
        length(state.marking["i"]),
        length(state.marking["r"]),
        state.current_time
        ),
        obs.cnt)
    obs.cnt+=1
    if obs.cnt>length(obs.observations)
        save(obs)
    end
end

function run_until(model, sampling, observer, end_time)
    #sampling=SampleSemiMarkov.FirstReaction()

    running=true
    trans=NRTransition(model.state.last_fired, current_time(model))
    steps=0
    while running && trans.time<end_time
    	trans=choose(sampling, model, rng)
    	#println(id_time)
    	if trans.time!=Inf
            #println("fire ", id_time)
    		fire(model, trans)
            #println("marking ", model.state.marking)
    	else
    		running=false
    	end
        observe(observer, model.state)
        steps+=1
    end
    println("steps ", steps, " end time ", trans.time)
end

observer=OutputObserver("sirs.txt")
sampling=NextReactionHazards()

@time(run_until(model, sampling, observer, 30.0))
save(observer)

