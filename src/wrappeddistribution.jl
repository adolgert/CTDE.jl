
"""
A `WrappedDistribution` combines a `ContinuousUnivariateDistribution`
with an enabling time in order to produce a distribution which
has an enabling time and can be sampled with a left shift,
as required for continuous-time simulation.
"""
type WrappedDistribution <: TransitionDistribution
    # Relative to the enabling time.
    relative_distribution::Distributions.ContinuousUnivariateDistribution
    enabling_time::Float64
    WrappedDistribution(d)=new(d, 0.0)
end


function Parameters(distribution::WrappedDistribution)
    params(distribution.relative_distribution)
end


function Parameters!(distribution::WrappedDistribution, params_list...)
    # Use the type of the wrapped distribution, which is immutable,
    # to create a new copy using the modified parameters
    dist_type=typeof(distribution.relative_distribution)
    distribution.relative_distribution=dist_type(params_list...)
end


function EnablingTime!(distribution::WrappedDistribution, Te::Float64)
    distribution.enabling_time=Te
end


"""
Given a distribution F
P(x<=t | x>t0)=P(x<=t && x>t0) / P(x>t0).
This is the left shift, as described by Gibson and Bruck.
"""
function Sample(distribution::WrappedDistribution, now::Float64,
        rng::MersenneTwister)
    quantile(distribution, now, rand(rng))
end


"""
Given a distribution F
P(x<=t | x>t0)=P(x<=t && x>t0) / P(x>t0)
U is a uniform variable between 0<U<=1
"""
function Quantile(distribution::WrappedDistribution, t0::Float64,
        U::Float64)
    te=distribution.enabling_time
    te+quantile(distribution.relative_distribution,
        U+(1-U)*cdf(distribution.relative_distribution, t0-te))
end


function MeasuredSample(d::WrappedDistribution, t0::Float64, rng)
    U=rand(rng)
    te=d.enabling_time
    value=te+quantile(d.relative_distribution,
        U+(1-U)*cdf(d.relative_distribution, t0-te))
    (value, -log(1-U))
end


"""
Cumulative distribution function.
The current time of the system is "now".
The cdf is being evaluated for a future time, "when".
"""
function CumulativeDistribution(dist::WrappedDistribution, when::Float64, now::Float64)
    t0te=cdf(dist.relative_distribution, now-dist.enabling_time)
    tte=cdf(dist.relative_distribution, when-dist.enabling_time)
    (tte-t0te)/(1-t0te)
end


function ConsumeSample(dist::WrappedDistribution, xa, start, finish)
    if xa<0
        xa=0
    end
    xa+HazardIntegral(dist, start, finish)
end


"""
Given a hazard, the integral of that hazard over a time.
int_{t0}^{t1} hazard(s, te) ds.
The hazard is f(t)/(1-F(t)), where F is the cumulative distribution
function and f is its derivative.
"""
function HazardIntegral(dist::WrappedDistribution, t1, t2)
    # logccdf is log(1-cdf(d, x))
    rel=dist.relative_distribution
    te=dist.enabling_time
    logccdf(rel, t1-te)-logccdf(rel, t2-te)
end


function Putative(dist::WrappedDistribution, when,
        interval, consumed_interval)
    ImplicitHazardIntegral(dist, interval-consumed_interval, when)
end

"""
This is the inverse of the hazard integral. At what time
would the cumulative hazard equal `xa`?
xa = int_{t0}^{t} hazard(s, te) ds. Solve for t.
"""
function ImplicitHazardIntegral(dist::WrappedDistribution, xa, t0)
    rel=dist.relative_distribution
    te=dist.enabling_time
    t=te+invlogccdf(rel, -xa+logccdf(rel, t0-te))
    @assert(t>=t0)
    t
end
