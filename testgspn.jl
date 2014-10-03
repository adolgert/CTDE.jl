include("semimarkov.jl")
using DataFrames
using Gadfly
using SemiMarkov
import SemiMarkov: enabled_transitions, current_time, current_time!
import SemiMarkov: fire, init


function kernel_model_exp()
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    structure=model.structure
    add_place(structure, "a")

    transition1=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionExponential(2.0, when), Int[])
        end)
    transition2=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionExponential(0.3, when), Int[])
        end)
    add_transition(structure, 1, transition1,
        [("a",-1),("a",1)],
        [])
    add_transition(structure, 2, transition2,
        [("a", -1), ("a", 1)],
        [])

    add_tokens(model.state.marking, "a", 1)

    model
end


function kernel_model_weib()
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    structure=model.structure
    add_place(structure, "a")

    transition1=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionWeibull(1.0, 1.7, when), Int[])
        end)
    transition2=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionWeibull(0.8, 0.4, when), Int[])
        end)
    add_transition(structure, 1, transition1,
        [("a",-1),("a",1)],
        [])
    add_transition(structure, 2, transition2,
        [("a", -1), ("a", 1)],
        [])

    add_tokens(model.state.marking, "a", 1)

    model
end


function independent_model_exp()
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    structure=model.structure
    add_place(structure, "a")
    add_place(structure, "b")

    transition1=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionExponential(2.0, when), Int[])
        end)
    transition2=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionExponential(0.3, when), Int[])
        end)
    add_transition(structure, 1, transition1,
        [("a",-1),("a",1)],
        [])
    add_transition(structure, 2, transition2,
        [("b", -1), ("b", 1)],
        [])

    add_tokens(model.state.marking, "a", 1)
    add_tokens(model.state.marking, "b", 1)

    model
end


function independent_model_weib()
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    structure=model.structure
    add_place(structure, "a")
    add_place(structure, "b")

    transition1=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionWeibull(1.0, 1.7, when), Int[])
        end)
    transition2=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionWeibull(0.8, 0.4, when), Int[])
        end)
    add_transition(structure, 1, transition1,
        [("a",-1),("a",1)],
        [])
    add_transition(structure, 2, transition2,
        [("b", -1), ("b", 1)],
        [])

    add_tokens(model.state.marking, "a", 1)
    add_tokens(model.state.marking, "b", 1)

    model
end


type EmpiricalObserver
	empirical::Array{EmpiricalDistribution,1}
	previous_time::Float64
	function EmpiricalObserver(cnt)
		ea=Array(EmpiricalDistribution,cnt)
		for x in 1:cnt
			ea[x]=EmpiricalDistribution()
		end
		new(ea, 0.0)
	end
end

function show(eo::EmpiricalObserver)
	total=sum([length(x) for x in eo.empirical])
    println("dist    percent")
	for i in 1:length(eo.empirical)
		println(i, " ", length(eo.empirical[i]), " ", length(eo.empirical[i])/total)
	end
end

function observe(eo::EmpiricalObserver, state)
	last_fired=state.last_fired
	delta=state.current_time-eo.previous_time
	push!(eo.empirical[last_fired], delta)
	eo.previous_time=state.current_time
end

function run_steps(model, sampling, report, end_step)
    #sampling=SampleSemiMarkov.FirstReaction()

    running=true
    init(sampling, model, rng)
    trans=NRTransition(model.state.last_fired, current_time(model))
    steps=0
    while running && steps<=end_step
    	trans=next(sampling, model, rng)
    	#println(id_time)
    	if trans.time!=Inf
            #println("fire ", id_time)
    		fire(sampling, model, trans, rng)
            #println("marking ", model.state.marking)
            report(model.state)
    	else
    		running=false
    	end
        steps+=1
    end
    #println("steps ", steps, " end time ", current_time(model))
end

function try_one(model, sampling)
    #sampling=NextReactionHazards()
    empirical=EmpiricalObserver(2)
    step_cnt=100000
    run_steps(model, sampling, x->observe(empirical, x), step_cnt)
    show(empirical)
end

seed=20
rng=MersenneTwister(seed)
model=kernel_model_exp()
sampling=FirstReaction()
try_one(model, sampling)

model=kernel_model_exp()
sampling=NextReactionHazards()
try_one(model, sampling)

model=kernel_model_weib()
sampling=FirstReaction()
try_one(model, sampling)

model=kernel_model_weib()
sampling=NextReactionHazards()
try_one(model, sampling)

model=independent_model_exp()
sampling=FirstReaction()
try_one(model, sampling)

model=independent_model_exp()
sampling=NextReactionHazards()
try_one(model, sampling)

model=independent_model_weib()
sampling=FirstReaction()
try_one(model, sampling)

model=independent_model_weib()
sampling=NextReactionHazards()
try_one(model, sampling)