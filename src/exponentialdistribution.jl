
"""
An exponential distribution with an enabling time. No,
it doesn't use the enabling time, but it has one
for consistency.
"""
type TransitionExponential <: TransitionDistribution
    relative_distribution::Distributions.Exponential # relative time
    enabling_time::Float64
end


function TransitionExponential(rate::Real, enabling_time::Real)
    dist=Distributions.Exponential(1.0/rate)
    TransitionExponential(dist, enabling_time)
end


function TransitionExponential(rate::Real)
    TransitionExponential(rate, 0.0)
end


Parameters(d::TransitionExponential)=[1.0/scale(d.relative_distribution),
        d.enabling_time]


function Parameters!(d::TransitionExponential, rate::Real, enabling_time::Real)
    d.dist=Distributions.Exponential(1.0/rate)
    d.enabling_time=enabling_time
end


function EnablingTime!(d::TransitionExponential, t::Float64)
    d.enabling_time=t
end


function Sample(distribution::TransitionExponential, now::Float64, rng)
    # We store the distribution for this call. Doing the inverse with
    # a log() is very slow compared to the Ziggurat method, which should
    # be available here.
    now+quantile(distribution.relative_distribution,
            rand(rng))
end

function MeasuredSample(d::TransitionExponential, now::Float64, rng)
    u=randexp(rng)
    value=now+u*scale(d.relative_distribution)
    (value, u)
end

function HazardIntegral(dist::TransitionExponential, start, finish)
    @assert(finish>=start)
    (finish-start)/scale(dist.relative_distribution)
end


function ConsumeSample(dist::TransitionExponential, xa, start, finish)
    if xa<0
        xa=0
    end
    xa+HazardIntegral(dist, start, finish)
end

function CumulativeDistribution(dist::TransitionExponential, when, now)
    1-exp(-hazard_integral(dist, now, when))
end


function ImplicitHazardIntegral(dist::TransitionExponential,
        cumulative_hazard, current_time)
    @assert(cumulative_hazard>=0)
    current_time+cumulative_hazard*scale(dist.relative_distribution)
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

