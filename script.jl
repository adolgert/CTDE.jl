include("semimarkov.jl")

using SemiMarkov

include("samplemodels.jl")
cnt=10
seed=32
beta=2.0
gamma=1.0
if length(ARGS)>0
    cnt=int(ARGS[1])
end
if length(ARGS)>1
    seed=int(ARGS[2])
end
if length(ARGS)>2
    beta=float(ARGS[3])
end
if length(ARGS)>3
    gamma=float(ARGS[4])
end

rng=MersenneTwister(seed)
model=sir_explicit(cnt, [beta, gamma])
R0=beta/gamma

null_observer(state)=nothing


function print_observer(state)
    m=state.marking
    #println("(",length(m,"s"),", ",length(m,"i"),", ",length(m,"r"),")")
end

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

function run_until(model, sampling, report, end_time)
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
        report(model.state)
        steps+=1
    end
    #println("steps ", steps, " end time ", current_time(model))
end

observer=OutputObserver("sirs.txt")
sampling=NextReactionHazards()

run_cnt=10000
final_size=Array(Int, run_cnt)
final_time=Array(Float64, run_cnt)
for i in 1:run_cnt
    model.state=TokenState(int_marking())
    add_tokens(model.state.marking, "s", cnt-1)
    add_tokens(model.state.marking, "i", 1)
    # state->observe(print_observer, state)
    run_until(model, sampling, print_observer, Inf)
    final_size[i]=length(model.state.marking, "r")
    final_time[i]=current_time(model)
end
println("R0 ", R0)
println("average size ", mean(final_size))
println("average final time ", mean(final_time))
println("figures from Linda Allen's paper")
println("http://eaton.math.rpi.edu/csums/papers/epidemic/allenstochasticepidemic.pdf")
println("R0\tN=20\tN=100")
vals=[0.5 1.0 2.0 5.0 10.0;
     1.76 3.34 8.12 15.66 17.98;
     1.93 6.10 38.34 79.28 89.98]
for i in 1:3
    println(vals[1,i],'\t', vals[2,i], '\t', vals[3,i])
end
#save(observer)

