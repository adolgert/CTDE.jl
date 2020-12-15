export Intensity, Enabled, Update!, Reset!, Sample, Putative, HazardIntegral
export PartialProcess, Time, AddTransition!, Init, Hazards, Fire!
export TransitionCount, MemoryIntensity, MemorylessIntensity

include("intensity.jl")

include("partialprocess_private.jl")


"""
The PartialProcess is responsible for providing enough information
to find the next state and time of the process.
"""
struct PartialProcess
	state
    time::Float64
    dependency_graph::DependencyGraph
    clocks::Vector{Clock}
    PartialProcess(state)=new(state, 0.0, DependencyGraph(), Vector{Clock}())
end

Time(pp::PartialProcess)=pp.time

TransitionCount(pp::PartialProcess)=length(pp.clocks)

"""
The dependencies, int_deps and fire_deps, have to be arrays of keys
into the state.
"""
function AddTransition!(pp::PartialProcess, intensity::Intensity,
		int_deps, firing::Function, fire_deps, name; sampler_args...)
	clock=Clock(IntegratedIntensity(intensity), firing, name, sampler_args)
	push!(pp.clocks, clock)
	AddIntensity!(pp.dependency_graph, clock, int_deps)
	AddFiring!(pp.dependency_graph, clock, fire_deps)
end


function Init(pp::PartialProcess)
	pp.time=0.0
	for clock in pp.clocks
		FireIntensity!(clock.intensity, pp.time, pp.state,
			IntensityProject(pp.dependency_graph, clock)...)
	end
end

"""
This iterates over transitions. It passes the whole "clock", which
is the transition, so that we have visibility into which transition is
firing, if we want. What's needed from the clock is its intensity, which
combines the distribution with a history of when it is on or off.
"""
function Hazards(f::Function, pp::PartialProcess, rng)
	for clock in pp.clocks
		if Enabled(clock.intensity)
			f(clock, pp.time, :Enabled, rng)
		end
	end
end


function Fire!(pp::PartialProcess, time, clock, rng, intensity_observer,
		state_observer)
	affected_clocks, affected_places=FiringProject!(
			pp.dependency_graph, clock, pp.state, clock.firing)
	fireupdate=FireIntensity!(clock.intensity, time, pp.state,
			IntensityProject(pp.dependency_graph, clock)...)
	intensity_observer(clock, time, :Fired, rng)
	if Enabled(clock.intensity)
		intensity_observer(clock, time, :Enabled, rng)
	end
	for affected in affected_clocks
		updated=Update!(affected.intensity, time, pp.state,
				IntensityProject(pp.dependency_graph, affected))
		if updated!=:Unmodified
			intensity_observer(affected, time, updated, rng)
		end
	end
	pp.time=time
	state_observer(pp.state, affected_places, clock.name, time)
end
