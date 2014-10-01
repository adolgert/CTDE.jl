# This file tests whether propagators are working by creating
# very simple models, one which resets with each step and one
# where the firing of transitions is independent.
# It also, along the way, shows what a propagator expects
# of the models.
include("semimarkov.jl")
using DataFrames
using Gadfly
using SemiMarkov
import SemiMarkov: all_transitions, modified_transitions, current_time, current_time!

# This model always has the same two distributions enabled with
# the last transition as the enabling time. It disables both
# and enables both at each step.
# It is the simplest test of the samplers.
# Transition keys are 1 and 2.
type AlwaysEnabledModel
	current_time::Float64
	distribution_cnt::Int
	enabling_time::Array{Float64,1}
	make_distribution::Function
	AlwaysEnabledModel(distributions)=new(0.0, 2, Float64[-1, -1], distributions)
end

function distribution(model::AlwaysEnabledModel, which::Int)
	model.make_distribution(which, model.enabling_time[which])
end

function current_time(m::AlwaysEnabledModel)
	m.current_time
end

function current_time!(m::AlwaysEnabledModel, t::Float64)
	m.current_time=t
end

function exponential_distribution(which::Int, when::Float64)
	if which==1
		return TransitionExponential(2.0, when)
	elseif which==2
		return TransitionExponential(0.3, when)
	end
end

function always_enabled_exponential()
	AlwaysEnabledModel(exponential_distribution)
end

function weibull_distribution(which::Int, when::Float64)
	if which==1
		return TransitionWeibull(1.0, 1.7, when)
	elseif which==2
		return TransitionWeibull(0.8, 0.4, when)
	end
end

function always_enabled_weibull()
	AlwaysEnabledModel(weibull_distribution)
end

function all_transitions(enable::Function, disable::Function,
        model::AlwaysEnabledModel, rng::MersenneTwister)
	for j in 1:model.distribution_cnt
		if model.enabling_time[j]>=0
			disable(j, model.current_time)
		end
	end
	for i in 1:model.distribution_cnt
		model.enabling_time[i]=model.current_time
		d=distribution(model, i)
		enable(i, d, model.current_time, rng)
	end	
end


function modified_transitions(model::AlwaysEnabledModel, enable::Function,
		disable::Function, rng::MersenneTwister)
	for di in 1:model.distribution_cnt
		if model.enabling_time[di]>=0
			disable(di, model.current_time)
		end
	end
	for ei in 1:model.distribution_cnt
		model.enabling_time[ei]=model.current_time
		d=distribution(model, ei)
		enable(ei, d, model.current_time, rng)
	end
end


function fire(model::AlwaysEnabledModel, id_time::NRTransition)
	model.current_time=id_time.time
	model.enabling_time[id_time.key]=-1
end


function plot_comparison(nelson::NelsonAalenDistribution,
		comparison, name::String)
	cumulative=Array(Float64, length(nelson.integrated_hazard))
	times=Array(Float64, length(nelson.integrated_hazard))
	fixed=Array(Float64, length(nelson.integrated_hazard))
	for i in 1:length(cumulative)
		when, how_much=cdf(nelson, i)
		times[i]=when
		cumulative[i]=how_much
		fixed[i]=cdf(comparison, when, 0.0)
	end
	df=DataFrame(NelsonAalen=cumulative, Times=times, Original=fixed)
	#myplot=plot(df, x="Times", y="Original", Geom.line)
	 myplot=plot(df, layer(x="Times", y="Original", Geom.line),
	 	layer(x="Times", y="NelsonAalen", Geom.line))
	draw(PDF("$(name).pdf", 4inch, 3inch), myplot)
end


function always_test(propagator, model, cnt)
	rng=MersenneTwister(1)
	empirical=Array(EmpiricalDistribution, model.distribution_cnt)
	for create_idx=1:model.distribution_cnt
		empirical[create_idx]=EmpiricalDistribution()
	end
	for i in 1:cnt
	    choice=choose(propagator, model, rng)
	    push!(empirical[choice.key], choice.time-model.current_time)
	    fire(model, choice)
	end
	total=sum([length(x) for x in empirical])
    println("dist\tpercent")
	for p_idx in 1:model.distribution_cnt
		println(p_idx, "\t", length(empirical[p_idx])/total)
	end
	nelsonaalen=multiple_measures(empirical)
end

# This is the actual test. It makes a set of PDFs of the
# distributions coming in and going out.
function test_always_enable()
	na=always_test(FirstReaction(), always_enabled_exponential(), 10000)
	plot_comparison(na[1], exponential_distribution(1, 0.0), "fexp1")
	plot_comparison(na[2], exponential_distribution(2, 0.0), "fexp2")
	na=always_test(NextReactionHazards(), always_enabled_exponential(), 10000)
	plot_comparison(na[1], exponential_distribution(1, 0.0), "nexp1")
	plot_comparison(na[2], exponential_distribution(2, 0.0), "nexp2")
	na=always_test(FirstReaction(), always_enabled_weibull(), 10000)
	plot_comparison(na[1], weibull_distribution(1, 0.0), "fweib1")
	plot_comparison(na[2], weibull_distribution(2, 0.0), "fweib2")
	na=always_test(NextReactionHazards(), always_enabled_weibull(), 10000)
	plot_comparison(na[1], weibull_distribution(1, 0.0), "nweib1")
	plot_comparison(na[2], weibull_distribution(2, 0.0), "nweib2")
end


# This is independently-firing transitions.
# Because they are independent, no firing will disable another
# transition.
type BlinkingLightsModel
	current_time::Float64
	distribution_cnt::Int
	enabling_time::Array{Float64,1}
	make_distribution::Function
	last_transition::NRTransition
	function BlinkingLightsModel(distributions)
		last=NRTransition(nothing, -1)
		enabling=Float64[0,0]
		new(0.0, 2, enabling, distributions, last)
	end
end


function distribution(model::BlinkingLightsModel, which::Int)
	model.make_distribution(which, model.enabling_time[which])
end

function current_time(m::BlinkingLightsModel)
	m.current_time
end

function current_time!(m::BlinkingLightsModel, t::Float64)
	m.current_time=t
end

function blinking_exponential()
	BlinkingLightsModel(exponential_distribution)
end

function blinking_weibull()
	BlinkingLightsModel(weibull_distribution)
end

function all_transitions(enable::Function, disable::Function,
        model::BlinkingLightsModel, rng::MersenneTwister)
	if model.last_transition.key!=nothing
		#disable(model.last_transition.key, model.current_time)
		model.enabling_time[model.last_transition.key]=current_time(model)
	end
	for i in 1:model.distribution_cnt
		d=distribution(model, i)
		enable(i, d, model.current_time, rng)
	end
end


function modified_transitions(model::BlinkingLightsModel, enable::Function,
		disable::Function, rng::MersenneTwister)
	which_fired=model.last_transition.key
	if which_fired!=nothing
		# The one that fired is already disabled.
		#disable(which_fired, model.current_time)
		model.enabling_time[which_fired]=current_time(model)
		d=distribution(model, which_fired)
		enable(which_fired, d, model.current_time, rng)
	else
		# Have to enable everything on the first time through.
		for i in 1:model.distribution_cnt
			d=distribution(model, i)
			enable(i, d, model.current_time, rng)
		end
	end
end


function fire(model::BlinkingLightsModel, id_time::NRTransition)
	model.current_time=id_time.time
	model.enabling_time[id_time.key]=-1
	model.last_transition=id_time
end


function plot_comparison(empir::EmpiricalDistribution,
		comparison, name::String)
	cumulative=Array(Float64, length(empir.samples))
	times=Array(Float64, length(empir.samples))
	fixed=Array(Float64, length(empir.samples))
	for i in 1:length(cumulative)
		when, how_much=cdf(empir, i)
		times[i]=when
		cumulative[i]=how_much
		fixed[i]=cdf(comparison, when, 0.0)
	end
	df=DataFrame(NelsonAalen=cumulative, Times=times, Original=fixed)
	#myplot=plot(df, x="Times", y="Original", Geom.line)
	 myplot=plot(df, layer(x="Times", y="Original", Geom.line),
	 	layer(x="Times", y="NelsonAalen", Geom.line))
	draw(PDF("$(name).pdf", 4inch, 3inch), myplot)
end


function blink_test(propagator, model, cnt)
	rng=MersenneTwister(1)
	empirical=Array(EmpiricalDistribution, model.distribution_cnt)
	for create_idx=1:model.distribution_cnt
		empirical[create_idx]=EmpiricalDistribution()
	end
	for i in 1:cnt
	    choice=choose(propagator, model, rng)
	    enabling=model.enabling_time[choice.key]
	    push!(empirical[choice.key], choice.time-enabling)
	    fire(model, choice)
	end
	for bidx in 1:model.distribution_cnt
		build!(empirical[bidx])
	end
	total=sum([length(x) for x in empirical])
    println("dist\tpercent")
	for p_idx in 1:model.distribution_cnt
		println(p_idx, "\t", length(empirical[p_idx])/total)
	end
	empirical
end


function test_blinking()
	na=blink_test(FirstReaction(), blinking_exponential(), 10000)
	plot_comparison(na[1], exponential_distribution(1, 0.0), "bfexp1")
	plot_comparison(na[2], exponential_distribution(2, 0.0), "bfexp2")
	na=blink_test(NextReactionHazards(), blinking_exponential(), 10000)
	plot_comparison(na[1], exponential_distribution(1, 0.0), "bnexp1")
	plot_comparison(na[2], exponential_distribution(2, 0.0), "bnexp2")
	na=blink_test(FirstReaction(), blinking_weibull(), 10000)
	plot_comparison(na[1], weibull_distribution(1, 0.0), "bfweib1")
	plot_comparison(na[2], weibull_distribution(2, 0.0), "bfweib2")
	na=blink_test(NextReactionHazards(), blinking_weibull(), 10000)
	plot_comparison(na[1], weibull_distribution(1, 0.0), "bnweib1")
	plot_comparison(na[2], weibull_distribution(2, 0.0), "bnweib2")
end


test_always_enable()
test_blinking()
