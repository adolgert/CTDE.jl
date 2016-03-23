using Graphs

include("partialprocess_private.jl")

"""
The PartialProcess is responsible for providing enough information
to find the next state and time of the process.
"""
type PartialProcess
    time::Float64
    dependency_graph::DependencyGraph
    clocks::Vector{Clocks}
    PartialProcess()=new(0.0, DependencyGraph(), Vector{Clocks}())
end

Time(pp::PartialProcess)=pp.time


function AddTransition!(pp::PartialProcess, intensity, int_deps,
		firing, fire_deps, name)
	clock=Clock(intensity, firing, name)
	push!(pp.clocks, clock)
	AddIntensity!(pp.dependency_graph, intensity, int_deps)
	AddFiring!(pp.dependency_graph, firing, fire_deps)
end


function Init(pp::PartialProcess)
	for clock in pp.clocks
		FireIntensity!(clock, pp.time,
			IntensityProject(pp.dependency_graph, clock))
	end
end

function Hazards(f::Function, pp::PartialProcess)
	for clock in pp.clocks
		if Enabled(clock)
			f(clock)
		end
	end
end


function Fire!(pp::PartialProcess, time, clock, rng, intensity_observer)
	affected_clocks=FiringProject!(pp.dependency_graph, clock)
	FireIntensity!(clock, time, IntensityProject(pp.dependency_graph, clock))
	intensity_observer(clock, time, rng)
	for affected in affected_clocks
		updated=UpdateIntensity!(affected, time,
				IntensityProject(pp.dependency_graph, clock))
		if updated
			intensity_observer(affected, time, rng)
		end
	end
	pp.time=time
end


# SIR Example
# This is a simple one.

type SpeciesCount
	v::Int
	SpeciesCount(n)=new(n)
end

type SIR
	s::SpeciesCount
	i::SpeciesCount
	r::SpeciesCount
	SIR()=new(SpeciesCount(1), SpeciesCount(0), SpeciesCount(0))
end

function MoveToken(a, b)
	b.v=a.v
	a.v=0
end

type SpeciesIntensity
	distribution
	enabled::Bool
	SpeciesIntensity(distribution)=new(distribution, false)
end

Enabled(si::SpeciesIntensity)=si.enabled

function Update!(si::SpeciesIntensity, time, state...)
	enabled_now=true
	for s in state
		if s.v<1
			enabled_now=false
		end
	end

	if enabled_now != si.enabled
		if enabled_now
			distparams.enabling_time=time
		else
			0
		end
		si.enabled=enabled_now
		return true
	else
		return false
	end
end

function MakeSIR(N=10)
	state=Vector{SIR}(N)
	for i = 1:N
		state[i]=SIR()
	end
	state[1].s=0
	state[1].i=1

	process=PartialProcess()
end

type A
  v::Int
end

immutable B
  v::Int
end


function TestHash()
	a1=A(3)
	a2=A(4)
	b1=B(1)
	b2=B(2)
	pp=PartialProcess(a1)
	dg=DependencyGraph()
	AddIntensity!(dg, a1, [a2, b2, b1])
	AddFiring!(dg, a1, [a2, b2, b1])
end

TestHash()
