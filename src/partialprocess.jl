include("partialprocess_private.jl")

"""
An `Intensity` is an abstract class for intensities, or hazard rates.
This defines some default behavior. The intensity is responsible for
knowing, given the past states since this transition last fired,
what is the current distribution going forward.
"""
abstract Intensity

Enabled(intensity::Intensity)=intensity.enabled

Sample(intensity::Intensity, when, rng)=rand(intensity.distribution, when, rng)

function Putative(intensity::Intensity, when, exponential_interval)
	implicit_hazard_integral(intensity.distribution, exponential_interval, when)
end

function HazardIntegral(intensity::Intensity, time0, time1)
	hazard_integral(intensity.distribution, time0, time1)
end


"""
The PartialProcess is responsible for providing enough information
to find the next state and time of the process.
"""
type PartialProcess
    time::Float64
    dependency_graph::DependencyGraph
    clocks::Vector{Clock}
    PartialProcess()=new(0.0, DependencyGraph(), Vector{Clock}())
end

Time(pp::PartialProcess)=pp.time


function AddTransition!(pp::PartialProcess, intensity::Intensity,
		int_deps::Tuple, firing::Function, fire_deps::Tuple, name)
	clock=Clock(intensity, firing, name)
	push!(pp.clocks, clock)
	AddIntensity!(pp.dependency_graph, clock, int_deps)
	AddFiring!(pp.dependency_graph, clock, fire_deps)
end


function Init(pp::PartialProcess)
	for clock in pp.clocks
		FireIntensity!(clock, pp.time,
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


function Fire!(pp::PartialProcess, time, clock, rng, intensity_observer)
	affected_clocks=FiringProject!(pp.dependency_graph, clock, clock.firing)
	fireupdate=FireIntensity!(clock, time,
			IntensityProject(pp.dependency_graph, clock)...)
	intensity_observer(clock, time, :Disabled, rng)
	if Enabled(clock)
		intensity_observer(clock, time, :Enabled, rng)
	end
	for affected in affected_clocks
		updated=UpdateIntensity!(affected, time,
				IntensityProject(pp.dependency_graph, affected))
		if updated!=:Unmodified
			intensity_observer(affected, time, updated, rng)
		end
	end
	pp.time=time
end


function Dynamics(pp::PartialProcess, sampler, rng)
	time, clock=Next(sampler, pp, rng)
	if !isinf(time)
		Fire!(pp, time, clock, rng, Observer(sampler))
	end
	(time, clock)
end
