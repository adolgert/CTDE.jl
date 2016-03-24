using Logging
@Logging.configure(level=DEBUG)
include("samplesemimarkov.jl")
include("transitiondistributions.jl")
include("partialprocess.jl")
# SIR Example
# This is a simple one.

"""
We make a SpeciesCount object so that the count of species
is a type that is passed by reference.
"""
type SpeciesCount
	v::Int
	SpeciesCount(n)=new(n)
end

type StoichiometricOperator
	stoichiometry::Vector{Int}
end

function Fire(so::StoichiometricOperator, args...)
	for arg_idx = 1:length(args)
		args[arg_idx].v+=so.stoichiometry[arg_idx]
	end
	@debug("StoichiometricOperator.Fire $args")
	# No reason to have a stoichiometry of 0, so report all
	# arguments were modified.
	args
end

function BuildStoichiometricFiring(stoichiometry)
	so=StoichiometricOperator(stoichiometry)
	(args...)->Fire(so, args...)
end


type SpeciesIntensity <: Intensity
	distribution::TransitionExponential
	enabled::Bool
	SpeciesIntensity(distribution)=new(distribution, false)
end

Enabled(si::SpeciesIntensity)=si.enabled

function Update!(si::SpeciesIntensity, time, state...)
	modified=:Undefined
	enabled_now=true
	for s in state
		if s.v<1
			enabled_now=false
			break
		end
	end

	if enabled_now != si.enabled
		if enabled_now
			si.distribution.enabling_time=time
			modified=:Enabled
		else
			modified=:Disabled
		end
		si.enabled=enabled_now
	else
		modified=:Unmodified
	end
	modified
end


type SIR
	s::SpeciesCount
	i::SpeciesCount
	r::SpeciesCount
	SIR()=new(SpeciesCount(1), SpeciesCount(0), SpeciesCount(0))
end


function MakeSIR(N=10)
	state=Vector{SIR}(N)
	for i = 1:N
		state[i]=SIR()
	end
	state[1].s.v=0
	state[1].i.v=1

	process=PartialProcess()

	# Add recoveries.
	for ridx = 1:N
		hazard=SpeciesIntensity(TransitionExponential(1.0, 0.0))
		recovery=BuildStoichiometricFiring([-1, 1])
		AddTransition!(process, hazard, (state[ridx].i,),
				recovery, (state[ridx].i, state[ridx].r), "r$ridx")
	end

	# Add infections.
	for iidx = 1:N
		for sidx=1:N
			if sidx!=iidx
				hazard=SpeciesIntensity(TransitionExponential(1.5, 0.0))
				infection=BuildStoichiometricFiring([-1, 1])
				AddTransition!(process, hazard,
						(state[iidx].i, state[sidx].s), infection,
						(state[sidx].s, state[sidx].i), "i$iidx-$sidx")
			end
		end
	end
	(process, state)
end

function PrintState(state::Vector{SIR})
	for idx = 1:length(state)
		assert(-1<state[idx].s.v<2)
		if state[idx].s.v>0
			print("s")
		end
		assert(-1<state[idx].i.v<2)
		if state[idx].i.v>0
			print("i")
		end
		assert(-1<state[idx].r.v<2)
		if state[idx].r.v>0
			print("r")
		end
	end
	print("\n")
end

function TestHash()
	rng=MersenneTwister(333333)
	process, state=MakeSIR(3)
	Init(process)
	#sampler=FirstReaction()
	sampler=NextReactionHazards()
	time, clock=Dynamics(process, sampler, rng)
	while !isinf(time)
		print("$time $clock ")
		PrintState(state)
		time, clock=Dynamics(process, sampler, rng)
	end
	print("Successful completion\n")
end

TestHash()
