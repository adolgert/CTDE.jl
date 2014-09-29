
using Distributions
import Distributions: quantile, rand, cdf
import Base: rand

export TransitionDistribution, WrappedDistribution, TransitionExponential
export TransitionWeibull
export rand, test, hazard_integral, implicit_hazard_integral, cdf
export parameters, quantile

# These are the distributions of stochastic processes in absolute time.
# So it's exponential, starting at some enabling time, or
# Weibull, but with an enabling time. They can be sampled at any
# time with the assumption that they have not fired up until
# some time called "now."
abstract TransitionDistribution

type WrappedDistribution <: TransitionDistribution
    # Relative to the enabling time.
    relative_distribution::Distributions.ContinuousUnivariateDistribution
    enabling_time::Float64
    WrappedDistribution(d, e)=new(d, e)
end

# Given a distribution F
# P(x<=t | x>t0)=P(x<=t && x>t0) / P(x>t0)
function rand(distribution::WrappedDistribution, now::Float64,
        rng::MersenneTwister)
    quantile(distribution, now, rand(rng))
end

# Given a distribution F
# P(x<=t | x>t0)=P(x<=t && x>t0) / P(x>t0)
# U is a uniform variable between 0<U<=1
function quantile(distribution::WrappedDistribution, now::Float64,
        U::Float64)
    time_delta=now-distribution.enabling_time
    now+quantile(distribution.relative_distribution,
        1-(1-U)*ccdf(distribution.relative_distribution, time_delta))-time_delta
end

# The current time of the system is "now".
# The cdf is being evaluated for a future time, "when".
function cdf(dist::WrappedDistribution, when::Float64, now::Float64)
    1 - ((1-cdf(dist.relative_distribution, when-dist.enabling_time))/
            (1-cdf(dist.relative_distribution, now-dist.enabling_time)))
end

# Given a hazard, the integral of that hazard over a time.
# int_{t0}^{t1} hazard(s, te) ds
function hazard_integral(dist::WrappedDistribution, last, now)
    # logccdf is log(1-cdf(d, x))
    rel=dist.relative_distribution
    logccdf(rel, last-dist.enabling_time)-logccdf(rel, now-dist.enabling_time)
end

# xa = int_{t0}^{t} hazard(s, te) ds. Solve for t.
function implicit_hazard_integral(dist::WrappedDistribution, xa, t1)
    rel=dist.relative_distribution
    -invlogccdf(rel, xa-logccdf(t1-dist.enabling_time))
end


######### Exponential
type TransitionExponential <: TransitionDistribution
    relative_distribution::Distributions.Exponential # relative time
    enabling_time::Float64
end

function TransitionExponential(rate::Real, enabling_time::Real)
    dist=Distributions.Exponential(1.0/rate)
    TransitionExponential(dist, enabling_time)
end

parameters(d::TransitionExponential)=[1.0/scale(d.dist), d.enabling_time]

function rand(distribution::TransitionExponential, now::Float64, rng)
    # We store the distribution for this call. Doing the inverse with
    # a log() is very slow compared to the Ziggurat method, which should
    # be available here.
    now+quantile(distribution.relative_distribution,
            rand(rng))
end

function hazard_integral(dist::TransitionExponential, start, finish)
    (finish-start)/scale(dist.relative_distribution)
end

function cdf(dist::TransitionExponential, when, now)
    1-exp(-hazard_integral(dist, now, when))
end

function implicit_hazard_integral(dist::TransitionDistribution,
        cumulative_hazard, current_time)
    current_time+cumulative_hazard*scale(dist.relative_distribution)
end

function test(TransitionExponential)
    rng=MersenneTwister()
    rate=2.0
    dist=TransitionExponential(rate, 0.0)
    ed=EmpiricalDistribution()
    for i in 1:100000
        push!(ed, rand(dist, 0.0, rng))
    end
    lambda_estimator=1/Base.mean(ed.samples)
    too_low=(rate<lambda_estimator*(1-1.96/sqrt(length(ed))))
    too_high=(rate>lambda_estimator*(1+1.96/sqrt(length(ed))))
    @trace("TransitionExponential low ", too_low, " high ", too_high)
end

########### Weibull
# F(T)=1-exp(-((T-Te)/lambda)^k)
type TransitionWeibull <: TransitionDistribution
    parameters::Array{Float64,1}
end
function TransitionWeibull(lambda, k, enabling_time)
    TransitionWeibull([lambda, k, enabling_time])
end
parameters(tw::TransitionWeibull)=tw.parameters

function rand(dist::TransitionWeibull, now::Float64, rng::MersenneTwister)
    (λ, k, tₑ)=dist.parameters
    d=now-tₑ
    value=0
    U=Base.rand(rng)
    if d>0
        value=λ*(-log(1-U)+(d/λ)^k)^(1/k)-d
    else
        value=-d+λ*(-log(1-U))^(1/k)
    end
    now+value
end

function hazard_integral(dist::TransitionWeibull, last, now)
    (λ, k, tₑ)=dist.parameters
    ((now-tₑ)/λ)^k - ((last-tₑ)/λ)^k
end

function cdf(dist::TransitionWeibull, when, now)
    1-exp(-hazard_integral(dist, now, when))
end

function implicit_hazard_integral(dist::TransitionWeibull,
        cumulative_hazard, now)
    (λ, k, tₑ)=dist.parameters
    tₑ + λ*(cumulative_hazard + ((now-tₑ)/λ)^k)^(1.0/k)
end

function test(dist::TransitionWeibull)
    rng=MersenneTwister()
    (λ, k, tₑ)=dist.parameters
    ed=EmpiricalDistribution()
    for i in 1:10000
        push!(ed, rand(dist, 0.0, rng))
    end
    expected_mean=λ*gamma(1+1/k)
    actual_mean=mean(ed)
    @trace("mean expected ", expected_mean, " actual ", actual_mean,
        " diff ", abs(expected_mean-actual_mean))

    expected_variance=λ^2*(gamma(1+2/k)-gamma(1+1/k)^2)
    obs_var=variance(ed)
    @trace("variance expected ", expected_variance, " actual ", obs_var,
        " diff ", abs(expected_variance-obs_var))

    min_value=min(ed)
    mink=min_value^k
    total=0.0
    for i in 1:length(ed)
        total+=ed.samples[i]^k
    end
    λ_estimator=(total/length(ed)-mink)^(1/k)
    @trace("λ expected ", λ, " actual ", λ_estimator,
        " diff ", abs(λ-λ_estimator))

    numerator=0.0
    denominator=0.0
    logsum=0.0
    for s in ed.samples
        numerator+=s^k*log(s) - mink*log(min_value)
        denominator+=s^k - mink
        logsum+=log(s)
    end
    k_est_inv=numerator/denominator - logsum/length(ed)
    k_est=1.0/k_est_inv
    @trace("k expected ", k, " actual ", k_est,
        " diff ", abs(k-k_est))
end


#################################
type EmpiricalDistribution
    samples::Array{Float64,1}
    built::Bool
    EmpiricalDistribution()=new(Array(Float64,0), false)
end

function push!(ed::EmpiricalDistribution, value)
    push!(ed.samples, value)
end

length(ed::EmpiricalDistribution)=length(ed.samples)

function mean(ed::EmpiricalDistribution)
    Base.mean(ed.samples)
end

function min(ed::EmpiricalDistribution)
    Base.minimum(ed.samples)
end

function variance(ed::EmpiricalDistribution)
    m=Base.mean(ed.samples)
    total=0.0
    for d in ed.samples
        total+=(m+d)^2
    end
    total/length(ed.samples)
end


function kolmogorov_smirnov_statistic(ed::EmpiricalDistribution, other)
    sort!(ed.samples)
    sup_diff=0.0
    n=length(ed.samples)
    for i in 1:n
        t=ed.samples[i]
        F_empirical=i/n
        sup_diff=max(sup_diff, abs(F_empirical-cdf(other, t)))
    end
    c_alpha=1.36 # 0.05 confidence interval
    (sup_diff, sup_diff > c_alpha*sqrt(2/n))
end
