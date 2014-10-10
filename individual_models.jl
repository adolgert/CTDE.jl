include("herd.jl")
include("herd_plot.jl")
include("herd_observer.jl")

run_cnt=100_000
seed=32
if length(ARGS)>0
    run_cnt=int(ARGS[1])
end
if length(ARGS)>1
    seed=int(ARGS[2])
end

rng=MersenneTwister(seed)

function individual_graph_set(model, title, run_cnt)
    sampling=FirstReaction()
    observer=IndividualDiseaseObserver(run_cnt)

    run_steps(model, sampling, s->observe(observer, s), rng)
    println("Ran model")
    plot_observer(observer, title)
end


disease_exponential={
    'l'=>3.59, # latent to infectious
    's'=>2.04, # not clinical to clinical
    'i'=>4.39, # infectious to recovered
    'e'=>10, # clinical back to subclinical (unknown and less important)
    'g'=>0.5 # rate of infection of neighbor
}

disease_nonexponential={
    'l'=>(3.974, 1.782),    # 'l'=>(1.782, 3.974),
    's'=>(1.222, 1.672),
    'i'=>(3.969, 1.107),
    'u'=>(5.3, 4.02),
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

function run_individual_model()
    # model=individual_exponential(de_params, exponential_transition_distributions)
    # individual_graph_set(model, "Simulation from Exponential Fits to Data", run_cnt)
    # nonexp_model=individual_nonexponential(disease_nonexponential)
    # individual_graph_set(nonexp_model, "Simulation with Dependent Clinical Period", run_cnt)
    indep_model=individual_independent(disease_nonexponential)
    individual_graph_set(indep_model,
            "Simulation with Independent Clinical Period", run_cnt)
end

run_individual_model()
