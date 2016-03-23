
# Dependency Graph for causality in the process

type ClockAdjacency
    hazard::Vector{Any}
    firing::Vector{Any}
    ClockAdjacency()=new()
end

type DependencyGraph
	place::Dict{Any,Vector{Any}}
	clock::Dict{Any,ClockAdjacency}

	DependencyGraph()=new(Dict{Any,Vector{Any}}(), Dict{Any,ClockAdjacency}())
end


function AddIntensity!(dg::DependencyGraph, clock, places)
	if haskey(dg.clock, clock)
		dg.clock[clock].hazard=places
	else
		ca=ClockAdjacency()
		ca.hazard=places
		dg.clock[clock]=ca
	end

	for place_idx = 1:length(places)
		if haskey(dg.place, places[place_idx])
			push!(dg.place[places[place_idx]], clock)
		else
			places[place_idx]=[clock]
		end
	end
end

function AddFiring!(dg::DependencyGraph, clock, places)
	if haskey(dg.clock, clock)
		dg.clock[clock].firing=places
	else
		ca=ClockAdjacency()
		ca.firing=places
		dg.clock[clock]=ca
	end
end

function IntensityProject(dg::DependencyGraph, clock)
	dg.clock[clock].hazard
end

function FiringProject!(dg::DependencyGraph, clock, operator)
	affected_places=operator(dg.clock[clock].firing)
	affected_clocks=Set()
	for place in affected_places
		union!(affected_clocks, dg.place[place])
	end
	setdiff!(affected_clocks, clock)
	affected_clocks
end



# A Clock is a multiparameter clock that fires at intervals.
type Clock
	intensity
	firing
	last_modification_time::Float64
	integrated_hazard::Float64
	name
	Clock(intensity, firing, name)=new(intensity, firing,
		0.0, 0.0, name)
end

Enabled(c::Clock)=Enabled(c.intensity)

function FireIntensity!(c::Clock, time, state)
	c.last_modification_time=time
	c.integrated_hazard=0
	c.enabled=false
	Update!(c.intensity, time, state...)
end

function UpdateIntensity!(c::Clock, time, state)
	if Enabled(c.intensity)
		c.integrated_hazard+=HazardIntegral(c.intensity,
				c.last_modification_time, time)
	end
	Update!(c.intensity, time, state...)
end


function Sample(c::Clock, when, rng)
	Sample(c.intensity, when, rng)
end

function Putative(c::clock, when, exponential_interval)
	ImplicitHazardIntegral(exponential_interval-c.integrated_hazard, when)
end

