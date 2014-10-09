include("herd.jl")
include("herd_plot.jl")
include("herd_observer.jl")

cnt=20
run_cnt=100_000
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

disease_exponential={
    'l'=>3.59, # latent to infectious
    's'=>2.04, # not clinical to clinical
    'i'=>4.39, # infectious to recovered
    'e'=>10, # clinical back to subclinical (unknown and less important)
    'g'=>0.5 # rate of infection of neighbor
}

disease_nonexponential={
    'l'=>(1.782, 3.974),
    's'=>(1.222, 1.672),
    'i'=>(3.969, 1.107),
    'e'=>1/10,
    'g'=>1/0.5
}

# Convert scale in days to a hazard per day.
de_params=Dict{Char,Float64}()
for (k, v) in disease_exponential
    de_params[k]=1/v
end

exponential_transition_distributions={
    "infectious"=>(lm, when, p)->TransitionExponential(p['l'], when)
}

function herd_model(params, rng)
    params['c']=cnt
    model=individual_exponential_graph(de_params)
    sampling=FirstReaction()
    observer=HerdDiseaseObserver(params['c'])
    run_steps(model, sampling, s->observe(observer, s), rng)
    show(observer)
end

herd_model(de_params, rng)
