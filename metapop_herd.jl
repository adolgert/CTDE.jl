include("herd.jl")
include("herd_plot.jl")
include("herd_observer.jl")

using SemiMarkov.SmallGraphs

cnt=20
run_cnt=10
seed=32
if length(ARGS)>0
    cnt=int(ARGS[1])
end
if length(ARGS)>1
    seed=int(ARGS[2])
end

rng=MersenneTwister(seed)

disease_exponential={
    'l'=>3.59, # latent to infectious
    's'=>2.04, # not clinical to clinical
    'i'=>4.39, # infectious to recovered
    'e'=>5, # clinical back to subclinical (unknown and less important)
    'g'=>0.5, # rate of infection of neighbor
}

disease_nonexponential={
    'l'=>(1.782, 3.974),
    's'=>(1.222, 1.672),
    'i'=>(3.969, 1.107),
    'e'=>1/5,
    'g'=>1/0.5,
    'f'=>1/4, # rate of infection to next pen
    'w'=>1/20 # rate of infection to any other animal in domain.
}

function pen_graph(block_cnt, block_length)
    g=UndirectedGraph()
    block_base=1
    for block_idx in 1:block_cnt
        add_edge(g, block_base, block_base+1)
        for row in 2:block_length
            top=block_base+2*(row-1)
            bottom=block_base+2*(row-1)+1
            add_edge(g, top, bottom)
            add_edge(g, top, top-2)
            add_edge(g, bottom, bottom-2)
        end
        block_base+=2*block_length
    end
    g
end

# Convert scale in days to a hazard per day.
de_params=Dict{Char,Float64}()
for (k, v) in disease_exponential
    de_params[k]=1/v
end

function herd_model(params, individuals_per_pen, block_cnt, block_length, rng)
    pen_contact=pen_graph(block_cnt, block_length)
    model=explicit_metapopulation(params, pen_contact)
    total=length(pen_contact)*individuals_per_pen
    sampling=NextReactionHazards()
    observer=HerdDiseaseObserver(total)
    run_steps(model, sampling, s->observe(observer, s), rng)
    show(observer)
    save_observer(observer, "metapop.h5")
    plot_observer(observer, "Metapopulation Trajectory")
end

de_params['c']=cnt
disease_nonexponential['c']=cnt

block_cnt=2
block_length=2
herd_model(disease_nonexponential, cnt, block_cnt, block_length, rng)
