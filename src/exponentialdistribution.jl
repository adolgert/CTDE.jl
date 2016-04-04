
"""
An exponential distribution with an enabling time. No,
it doesn't use the enabling time, but it has one
for consistency.
"""
type TransitionExponential <: TransitionDistribution
    hazard::Float64
    enabling_time::Float64
    TransitionExponential(rate::Real)=new(rate, 0.0)
end


Parameters(d::TransitionExponential)=[d.hazard]


function Parameters!(d::TransitionExponential, rate::Real)
    d.hazard=rate
end


function EnablingTime!(d::TransitionExponential, t::Float64)
    d.enabling_time=t
end


function Sample(distribution::TransitionExponential, now::Float64, rng)
    now-randexp(rng)/distribution.hazard
end

function MeasuredSample(d::TransitionExponential, now::Float64, rng)
    u=randexp(rng)
    (now+u/d.hazard, u)
end

function HazardIntegral(dist::TransitionExponential, start, finish)
    @assert(finish>=start)
    (finish-start)*dist.hazard
end


function ConsumeSample(dist::TransitionExponential, xa, start, finish)
    xa = xa<0 ? 0 : xa
    xa+HazardIntegral(dist, start, finish)
end

function CumulativeDistribution(dist::TransitionExponential, when, now)
    1-exp(-hazard_integral(dist, now, when))
end


function ImplicitHazardIntegral(dist::TransitionExponential,
        cumulative_hazard, current_time)
    @assert(cumulative_hazard>=0)
    current_time+cumulative_hazard/dist.hazard
end


function Putative(dist::TransitionExponential, when,
        interval, consumed_interval)
    ImplicitHazardIntegral(dist, interval-consumed_interval, when)
end


function test(TransitionExponential)
    rng=MersenneTwister()
    rate=2.0
    dist=TransitionExponential(rate, 0.0)
    ed=EmpiricalDistribution()
    for i in 1:100000
        push!(ed, Sample(dist, 0.0, rng))
    end
    lambda_estimator=1/Base.mean(ed.samples)
    too_low=(rate<lambda_estimator*(1-1.96/sqrt(length(ed))))
    too_high=(rate>lambda_estimator*(1+1.96/sqrt(length(ed))))
    @debug("TransitionExponential low ", too_low, " high ", too_high)
end

