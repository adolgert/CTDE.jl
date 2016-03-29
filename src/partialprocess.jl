export Intensity, Enabled, Update!, Reset!, Sample, Putative, HazardIntegral
export PartialProcess, Time, AddTransition!, Init, Hazards, Fire!
export TransitionCount
include("partialprocess_private.jl")

"""
An `Intensity` is an abstract class for intensities, or hazard rates.
This defines some default behavior. The intensity is responsible for
knowing, given the past states since this transition last fired,
what is the current distribution going forward.
"""
abstract Intensity

"""
If the intensity defines an enabled member, this will
take care of the Enabled() function.
"""
Enabled(intensity::Intensity)=intensity.enabled

Update!(intensity::Intensity, time, state, keys...)=nothing

"""
Most intensities don't need to reset when the transition fires.
"""
function Reset!(intensity::Intensity, time, state, keys...)
	Update!(intensity, time, state, keys...)
end

"""
If the intensity defines a distribution member, this will do sampling.
"""
Sample(intensity::Intensity, when, rng)=rand(intensity.distribution, when, rng)

"""
If the intensity defines a distribution member, this will do sampling.
"""
function Putative(intensity::Intensity, when, exponential_interval)
	implicit_hazard_integral(intensity.distribution, exponential_interval, when)
end

"""
If the intensity defines a distribution member, this will do the
hazard integral.
"""
function HazardIntegral(intensity::Intensity, time0, time1)
	hazard_integral(intensity.distribution, time0, time1)
end


"""
The PartialProcess is responsible for providing enough information
to find the next state and time of the process.
"""
type PartialProcess
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
	clock=Clock(intensity, firing, name, sampler_args)
	push!(pp.clocks, clock)
	AddIntensity!(pp.dependency_graph, clock, int_deps)
	AddFiring!(pp.dependency_graph, clock, fire_deps)
end


function Init(pp::PartialProcess)
	for clock in pp.clocks
		FireIntensity!(clock, pp.time, pp.state,
			IntensityProject(pp.dependency_graph, clock)...)
	end
end

function Hazards(f::Function, pp::PartialProcess, rng)
	for clock in pp.clocks
		if Enabled(clock)
			f(clock, pp.time, :Enabled, rng)
		end
	end
end


function Fire!(pp::PartialProcess, time, clock, rng, intensity_observer,
		state_observer)
	affected_clocks, affected_places=FiringProject!(
			pp.dependency_graph, clock, pp.state, clock.firing)
	fireupdate=FireIntensity!(clock, time, pp.state,
			IntensityProject(pp.dependency_graph, clock)...)
	intensity_observer(clock, time, :Fired, rng)
	if Enabled(clock)
		intensity_observer(clock, time, :Enabled, rng)
	end
	for affected in affected_clocks
		updated=UpdateIntensity!(affected, time, pp.state,
				IntensityProject(pp.dependency_graph, affected))
		if updated!=:Unmodified
			intensity_observer(affected, time, updated, rng)
		end
	end
	pp.time=time
	state_observer(pp.state, affected_places, clock.name, time)
end
