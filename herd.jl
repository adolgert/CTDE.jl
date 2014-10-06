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
    'g'=>0.5 # rate of infection of neighbor
}

disease_nonexponential={
    'l'=>(1.782, 3.974),
    's'=>(1.222, 1.672),
    'i'=>(3.969, 1.107),
    'e'=>1/10,
    'g'=>1/0.5
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


function individual_nonexponential(params)
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
            (TransitionWeibull(params['l'][1], params['l'][2], when), Int[])
        end)
    recover=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionGamma(params['i'][1], params['i'][2], when), Int[])
        end)
    add_transition(structure, 'f', infectious,
        [('l',-1),('i',1)],
        [])
    add_transition(structure, 'd', recover,
        [('i', -1), ('r', 1)],
        [])

    # The clinical state.
    clinical=ConstExplicitTransition( (lm, when)->begin
        (TransitionGamma(params['s'][1], params['s'][2], when), Int[])
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
	bandwidth=0.2
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
	myplot=plot(df, x="Times", y="PDF", Geom.line,
        Guide.xlabel("time interval since infection [days]"),
        Guide.ylabel("probability distribution of firing"),
        Guide.title(name))
	draw(SVG("$(name).svg", 20cm, 15cm), myplot)
end

function plot_one_observer(d, name)
    sort!(d)
    plot_density(d, name)
end

function plot_observer(eo::IndividualDiseaseObserver)
	plot_one_observer(eo.t[:,1], "time to infectious")
    plot_one_observer(eo.t[:,2], "time to clinical")
    plot_one_observer(eo.t[:,3], "time to recovered")
    plot_one_observer(eo.t[:,4], "time to subclinical")
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


function individual_exponential_graph(params)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    cnt=params['c']
    structure=model.structure
    for i in 1:cnt
        add_place(structure, (i,'s')) # susceptible
        add_place(structure, (i,'l')) # latent
        add_place(structure, (i,'i')) # infectious
        add_place(structure, (i,'r')) # recovered
        add_place(structure, (i,'n')) # not clinical
        add_place(structure, (i,'c')) # clinical
    end

    # The basic SIR
    for i in 1:cnt
        infectious=ConstExplicitTransition(
            (lm, when)->begin
                (TransitionExponential(params['l'], when), Int[])
            end)
        recover=ConstExplicitTransition(
            (lm, when)->begin
                (TransitionExponential(params['i'], when), Int[])
            end)
        add_transition(structure, (i,'f'), infectious,
            [((i,'l'),-1),((i,'i'),1)],
            [])
        add_transition(structure, (i,'d'), recover,
            [((i,'i'), -1), ((i,'r'), 1)],
            [])

        # The clinical state.
        clinical=ConstExplicitTransition( (lm, when)->begin
            (TransitionExponential(params['s'], when), Int[])
            end)
        endclinical=ConstExplicitTransition( (lm, when)->begin
            (TransitionExponential(params['e'], when), Int[])
            end)
        add_transition(structure, (i,'a'), clinical,
            [((i,'n'), -1), ((i,'i'), -1), ((i,'c'), 1), ((i,'i'), 1)],
            [])
        add_transition(structure, (i,'e'), endclinical,
            [((i,'c'), -1), ((i,'n'), 1)],
            [])
    end

    for i in 1:cnt
        for j in 1:cnt
            infect=ConstExplicitTransition((lm, when)->begin
                (TransitionExponential(params['g']/cnt, when), Int[])
                end)
            add_transition(structure, (i,j,'g'), infect,
            [((i,'i'),-1), ((j,'s'), -1), ((i,'i'),1), ((j,'l'),1)],
            [])
        end
    end

    for i in 1:(cnt-1)
        add_tokens(model.state.marking, (i,'s'), 1)
        add_tokens(model.state.marking, (i,'n'), 1)
    end
    add_tokens(model.state.marking, (cnt,'i'), 1)
    add_tokens(model.state.marking, (cnt,'n'), 1)

    model
end

typealias Time Float64
typealias TrajectoryEntry (Int64,Int64,Int64,Float64)
type TrajectoryStruct
  s::Int64
  i::Int64
  r::Int64
  t::Time
  TrajectoryStruct()=new(0,0,0,0.)
  TrajectoryStruct(s_,i_,r_,t_)=new(s_,i_,r_,t_)
end

function convert(::Type{TrajectoryEntry}, x::TrajectoryStruct)
    (x.s, x.i, x.r, x.t)
end

type HerdDiseaseObserver
    t::Array{TrajectoryEntry,1}
    cnt::Int
    sir::TrajectoryStruct
    previous_time::Float64
    HerdDiseaseObserver(cnt)=new(Array(TrajectoryEntry, 10_000),
            1, TrajectoryStruct(int(cnt-1), 1, 0, 0.), 0.)
end

function observe(eo::HerdDiseaseObserver, state)
    last_fired=state.last_fired
    delta=state.current_time-eo.previous_time
    const ind_fired={ 'f'=>1, 'a'=>2, 'd'=>3, 'e'=>4, 'z'=>5, 'g'=>6 }
    shiftstore= xo->begin
        xo.sir.t=state.current_time
        xo.t[xo.cnt]=convert(TrajectoryEntry, xo.sir)
        xo.cnt+=1
    end
    if length(last_fired)==2
        if last_fired[2]=='d'
            eo.sir.i-=1
            eo.sir.r+=1
            shiftstore(eo)
        end
    else
        eo.sir.s-=1
        eo.sir.i+=1
        shiftstore(eo)
    end
    if eo.cnt>length(eo.t)
        new_len=2*eo.cnt
        new_t=Array(TrajectoryEntry, new_len)
        new_t[1:length(eo.t)]=eo.t
        eo.t=new_t
    end
    eo.previous_time=state.current_time
    true # run until you can't run any more.
end

function show(eo::HerdDiseaseObserver)
  for i in 1:(eo.cnt-1)
      println(eo.t[i][1], " ", eo.t[i][2], " ", eo.t[i][3], " ", eo.t[i][4])
  end
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

function individual_graph_set(model)
    sampling=FirstReaction()
    observer=IndividualDiseaseObserver()

    run_steps(model, sampling, s->observe(observer, s), rng)
    println("Ran model")
    plot_observer(observer)
end

# Convert scale in days to a hazard per day.
de_params=Dict{Char,Float64}()
for (k, v) in disease_exponential
    de_params[k]=1/v
end

function run_individual_model()
    #model=individual_exponential(de_params)
    #individual_graph_set(exp_model)
    nonexp_model=individual_nonexponential(disease_nonexponential)
    individual_graph_set(nonexp_model)
end

function herd_model(params, rng)
    params['c']=cnt
    model=individual_exponential_graph(de_params)
    sampling=FirstReaction()
    observer=HerdDiseaseObserver(params['c'])
    run_steps(model, sampling, s->observe(observer, s), rng)
    show(observer)
end

#herd_model(de_params, rng)
run_individual_model()
