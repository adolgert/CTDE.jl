include("semimarkov.jl")
using DataFrames
using Gadfly
using SmoothingKernels
using SemiMarkov
import SemiMarkov: enabled_transitions, current_time, current_time!
import SemiMarkov: fire, init


disease_exponential={
	'l'=>3.59, # latent to infectious
	's'=>2.04, # not clinical to clinical
	'i'=>4.39, # infectious to recovered
	'e'=>10, # clinical back to subclinical (unknown and less important)
}


function individual_exponential(params)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    structure=model.structure
    add_place(structure, 'l') # latent
    add_place(structure, 'i') # infectious
    add_place(structure, 'r') # recovered
    add_place(structure, 'n') # not clinical
    add_place(structure, 'c') # clinical

    # The basic SIR
    infectious=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionExponential(params['l'], when), Int[])
        end)
    recover=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionExponential(params['i'], when), Int[])
        end)
    add_transition(structure, 'f', infectious,
        [('l',-1),('i',1)],
        [])
    add_transition(structure, 'd', recover,
        [('i', -1), ('r', 1)],
        [])

    # The clinical state.
    clinical=ConstExplicitTransition( (lm, when)->begin
    	(TransitionExponential(params['s'], when), Int[])
    	end)
    endclinical=ConstExplicitTransition( (lm, when)->begin
    	(TransitionExponential(params['e'], when), Int[])
    	end)
    add_transition(structure, 'a', clinical,
    	[('n', -1), ('i', -1), ('c', 1), ('i', 1)],
    	[])
    add_transition(structure, 'e', endclinical,
    	[('c', -1), ('n', 1)],
    	[])

    reset=ConstExplicitTransition( (lm, when)->begin
    	(TransitionExponential(1.0, when), Int[])
	end)
    add_transition(structure, 'z', endclinical,
    	[('n', -1), ('r', -1), ('n', 1), ('l', 1)],
    	[])


    add_tokens(model.state.marking, 'l', 1)
    add_tokens(model.state.marking, 'n', 1)

    model
end

type IndividualDiseaseObserver
	t::Array{Float64,2}
	latent::Float64 # When latent was started.
	previous_time::Float64
	idx
	IndividualDiseaseObserver()=new(zeros(Float64, 10_000, 4), 0., 0., 1)
end

function observe(eo::IndividualDiseaseObserver, state)
	last_fired=state.last_fired
	delta=state.current_time-eo.previous_time
	const ind_fired={ 'f'=>1, 'a'=>2, 'd'=>3, 'e'=>4, 'z'=>5 }
	idx=ind_fired[state.last_fired]
	# nonzero=Char[]
	# for p in Char['l', 'i', 'r', 'n', 'c']
	# 	if length(state.marking, p)>0
	# 		push!(nonzero, p)
	# 	end
	# end
    # println(join(nonzero), " ", eo.idx, " ", idx, " ", state.current_time, " ", eo.latent)
	if idx<5
		eo.t[eo.idx, idx]=state.current_time-eo.latent
	elseif idx==5
		# Back to the beginning.
		eo.latent=state.current_time
    	eo.idx+=1
	end

	eo.previous_time=state.current_time
	eo.idx<=size(eo.t)[1]
end


function plot_density(v, name)
	plot_cnt=100
	cumulative=zeros(Float64, plot_cnt)
	times=zeros(Float64, plot_cnt)
	bandwidth=2
	lambda=1/bandwidth
	vmax=maximum(v)
	println("density ", length(v), " ", minimum(v), " ", maximum(v), " ", mean(v))
	for i in 1:plot_cnt
		dx=(i-1)*1.1*vmax/plot_cnt
		s=lambda*SmoothingKernels.kernels[:epanechnikov](lambda*(v-dx))
		cumulative[i]=sum(s)
		times[i]=dx
	end
	df=DataFrame(PDF=cumulative, Times=times)
	#myplot=plot(df, x="Times", y="Original", Geom.line)
	myplot=plot(df, x="Times", y="PDF", Geom.line)
	draw(PDF("$(name).pdf", 4inch, 3inch), myplot)
end

function plot_observer(eo::IndividualDiseaseObserver)
	d=vec(eo.t[1,:])
	sort!(d)
	plot_density(d, "time to infectious")
end

function run_steps(model, sampling, report, rng)
    #sampling=SampleSemiMarkov.FirstReaction()

    running=true
    init(sampling, model, rng)
    trans=NRTransition(model.state.last_fired, current_time(model))
    steps=0
    while running
    	trans=next(sampling, model, rng)
    	#println(id_time)
    	if trans.time!=Inf
            #println("fire ", id_time)
    		fire(sampling, model, trans, rng)
            #println("marking ", model.state.marking)
            running=report(model.state)
    	else
    		running=false
    	end
        steps+=1
    end
    #println("steps ", steps, " end time ", current_time(model))
end

cnt=20
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

# Convert scale in days to a hazard per day.
de_params=Dict{Char,Float64}()
for (k, v) in disease_exponential
	de_params[k]=1/v
end
model=individual_exponential(de_params)
sampling=FirstReaction()
observer=IndividualDiseaseObserver()

run_steps(model, sampling, s->observe(observer, s), rng)
println("Ran model")
plot_observer(observer)
