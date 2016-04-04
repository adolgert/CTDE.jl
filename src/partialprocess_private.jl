import Base: show

# A Clock is a multiparameter clock that fires at intervals.
type Clock
	intensity::Intensity
	firing
	name # A name assigned by the client API.
	kind # classification of this clock for the sampler.
	Clock(intensity, firing, name, sampler_args)=new(intensity, firing,
		name, sampler_args)
end

show(io::IO, c::Clock)=show(io, c.name)


# Dependency Graph for causality in the process
type ClockAdjacency
    hazard
    firing
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
	(affected_clocks, affected_places)
end

