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


function run_until(model, sampling, end_time)
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
        steps+=1
    end
    println("steps ", steps, " end time ", trans.time)
end

sampling=NextReactionHazards()

@time(run_until(model, sampling, 30.0))

