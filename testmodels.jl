include("semimarkov.jl")
using SemiMarkov
import SemiMarkov: all_transitions, modified_transitions, current_time, current_time!

# This model always has the same two distributions enabled with
# the last transition as the enabling time.
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
	model.make_distribution(which, model.current_time)
end

function current_time(m::AlwaysEnabledModel)
	m.current_time
end

function current_time!(m::AlwaysEnabledModel, t::Float64)
	m.current_time=t
end

function exponential_distribution(which::Int, when::Float64)
	if which==1
		return TransitionExponential(1.0, when)
	elseif which==2
		return TransitionExponential(0.8, when)
	end
end

function always_enabled_exponential()
	AlwaysEnabledModel(exponential_distribution)
end

function weibull_distribution(which::Int, when::Float64)
	if which==1
		return TransitionWeibull(1.0, 1.2, when)
	elseif which==2
		return TransitionWeibull(0.8, 1.1, when)
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
		d=distribution(model, i)
		model.enabling_time[i]=model.current_time
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

always_test(FirstReaction(), always_enabled_exponential(), 10000)
always_test(NextReactionHazards(), always_enabled_exponential(), 10000)
always_test(FirstReaction(), always_enabled_weibull(), 10000)
always_test(NextReactionHazards(), always_enabled_weibull(), 10000)
