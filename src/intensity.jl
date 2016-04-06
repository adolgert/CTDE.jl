
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

Distribution(intensity::Intensity)=intensity.distribution

"""
If the intensity defines a distribution member, this will do sampling.
"""
Sample(intensity::Intensity, when, rng)=Sample(intensity.distribution, when, rng)

"""
If the intensity defines a distribution member, this will do sampling.
"""
function Putative(intensity::Intensity, when, exponential_interval, used)
	Putative(intensity.distribution, when, exponential_interval, used)
end

"""
If the intensity defines a distribution member, this will do the
hazard integral.
"""
function HazardIntegral(intensity::Intensity, time0, time1)
	HazardIntegral(intensity.distribution, time0, time1)
end


"""
The IntegratedIntensity wraps an intensity in order to
ensure that its integrated hazard is always calculated.
"""
type IntegratedIntensity <: Intensity
	last_modification_time::Float64
	integrated_hazard::Float64
	intensity::Intensity
	IntegratedIntensity(wrapped)=new(0.0, 0.0, wrapped)
end

Enabled(ii::IntegratedIntensity)=Enabled(ii.intensity)

function FireIntensity!(ii::IntegratedIntensity, time, state, keys...)
	@debug("Reset modification time for $(c.name)")
	ii.last_modification_time=time
	ii.integrated_hazard=-1
	Reset!(ii.intensity, time, state, keys...)
end

function Update!(ii::IntegratedIntensity, time, state, keys)
	if Enabled(ii.intensity)
		ii.integrated_hazard=ConsumeSample(Distribution(ii.intensity),
				ii.integrated_hazard, ii.last_modification_time, time)
		@debug("Added $added to integrated hazard of $(c.name)")
	end
	ii.last_modification_time=time
	Update!(ii.intensity, time, state, keys...)
end

function Reset!(ii::IntegratedIntensity, time, state, keys...)
	Reset!(ii.intensity, time, state, keys...)
end

Distribution(ii::IntegratedIntensity)=Distribution(ii.intensity)

function Sample(ii::IntegratedIntensity, when, rng)
	Sample(ii.intensity, when, rng)
end

function MeasuredSample(ii::IntegratedIntensity, when, rng)
	MeasuredSample(Distribution(ii.intensity), when, rng)
end

function Putative(ii::IntegratedIntensity, when, exponential_interval)
	Putative(ii.intensity, when, exponential_interval, ii.integrated_hazard)
end


"""
The MemoryIntensity follows a set pattern. There are three parts,
the invariant to determine when it is enabled, the parameters
function, which sets parameters from state, and the distribution
whose parameters are set.
"""
type MemoryIntensity <: Intensity
  invariant::Function
  distribution::TransitionDistribution
  enabled::Bool
  MemoryIntensity(invariant, distribution)=new(invariant, distribution, false)
end

Enabled(mi::MemoryIntensity)=mi.enabled

function Reset!(mi::MemoryIntensity, time, state, keys...)
	(mi.enabled, params)=mi.invariant(time, state, keys...)
	if mi.enabled
		Parameters!(mi.distribution, params...)
		EnablingTime!(mi.distribution, time)
	end
end

function Update!(mi::MemoryIntensity, time, state, keys...)
	modified=:Undefined
	(enabled, params)=mi.invariant(time, state, keys...)
	if mi.enabled!=enabled
		if enabled
			Parameters!(mi.distribution, params...)
			EnablingTime!(mi.distribution, time)
			modified=:Enabled
		else
			modified=:Disabled
		end
		mi.enabled=enabled
	else
		if params!=Parameters(mi.distribution)
			Parameters!(mi.distribution, params...)
			modified=:Modified
		else
			modified=:Unmodified
		end
	end
	modified
end



"""
The MemorylessIntensity follows a set pattern. There are three parts,
the invariant to determine when it is enabled, the parameters
function, which sets parameters from state, and the distribution
whose parameters are set.

This differs from MemoryIntensity in one line. If an enabled
intensity sees a change of state which changes its parameters,
the MemorylessIntensity changes the enabling time to the current
time, while the MemoryIntensity does not.
"""
type MemorylessIntensity <: Intensity
  invariant::Function
  distribution::TransitionDistribution
  enabled::Bool
  MemorylessIntensity(invariant, distribution)=new(
  		invariant, distribution, false)
end

Enabled(mi::MemorylessIntensity)=mi.enabled

function Reset!(mi::MemorylessIntensity, time, state, keys...)
	(mi.enabled, params)=mi.invariant(time, state, keys...)
	if mi.enabled
		Parameters!(mi.distribution, params...)
		EnablingTime!(mi.distribution, time)
	end
end

function Update!(mi::MemorylessIntensity, time, state, keys...)
	modified=:Undefined
	(enabled, params)=mi.invariant(time, state, keys...)
	if mi.enabled!=enabled
		if enabled
			Parameters!(mi.distribution, params...)
			EnablingTime!(mi.distribution, time)
			modified=:Enabled
		else
			modified=:Disabled
		end
		mi.enabled=enabled
	else
		if params!=Parameters(mi.distribution)
			Parameters!(mi.distribution, params...)
			# This is the only line different from MemoryIntensity.
			EnablingTime!(mi.distribution, time)
			modified=:Modified
		else
			modified=:Unmodified
		end
	end
	modified
end
