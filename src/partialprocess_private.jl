import Base: show

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

show(io::IO, c::Clock)=show(io, c.name)

Enabled(c::Clock)=Enabled(c.intensity)

function FireIntensity!(c::Clock, time, state, keys...)
	# @debug("FireIntensity state $state")
	c.last_modification_time=time
	c.integrated_hazard=0
	Update!(c.intensity, time, state, keys...)
end

function UpdateIntensity!(c::Clock, time, state, keys)
	if Enabled(c.intensity)
		c.integrated_hazard+=HazardIntegral(c.intensity,
				c.last_modification_time, time)
	end
	Update!(c.intensity, time, state, keys...)
end


function Sample(c::Clock, when, rng)
	Sample(c.intensity, when, rng)
end

function Putative(c::Clock, when, exponential_interval)
	Putative(c.intensity, when, exponential_interval-c.integrated_hazard)
end


# Dependency Graph for causality in the process
type ClockAdjacency
    hazard::Tuple
    firing::Tuple
    ClockAdjacency()=new()
end

type DependencyGraph
	place::Dict{Any,Vector{Clock}}
	clock::Dict{Clock,ClockAdjacency}

	DependencyGraph()=new(Dict{Any,Vector{Clock}}(), Dict{Clock,ClockAdjacency}())
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
			dg.place[places[place_idx]]=[clock]
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

	# Even if nothing depends on places affected by firing,
	# they should still be vertices in the graph.
	for place in places
		if !haskey(dg.place, place)
			dg.place[place]=[]
		end
	end
end

function IntensityProject(dg::DependencyGraph, clock)
	dg.clock[clock].hazard
end

function FiringProject!(dg::DependencyGraph, clock, state, operator)
	affected_places=operator(state, dg.clock[clock].firing...)
	affected_clocks=Set()
	for place in affected_places
		union!(affected_clocks, dg.place[place])
	end
	setdiff!(affected_clocks, [clock])
	affected_clocks
end

