
using Distributions
import Distributions: quantile, rand, cdf, logccdf, invlogccdf
import Base: rand, push!, isless, length

export TransitionDistribution, WrappedDistribution, TransitionExponential
export TransitionWeibull, TransitionGamma, TransitionLogLogistic
export rand, test, hazard_integral, implicit_hazard_integral, cdf
export parameters, quantile, EnablingTime!
export EmpiricalDistribution, push!, build!
export NelsonAalenDistribution, multiple_measures
export UniformDistribution, DiracDistribution

"""
These are the distributions of stochastic processes in absolute time.
So it's exponential, starting at some enabling time, or
Weibull, but with an enabling time. They can be sampled at any
time with the assumption that they have not fired up until
some time called "now."
"""
abstract TransitionDistribution


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
    WrappedDistribution(d, e)=new(d, e)
end


function Parameters(distribution::WrappedDistribution)
    params(distribution.relative_distribution)
end


function Parameters!(distribution, params_list)
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

function MeasuredSample(distribution::TransitionExponential, now::Float64, rng)
    u=randexp(rng)
    value=now+u*(1.0/scale(d.dist))
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


function ImplicitHazardIntegral(dist::TransitionDistribution,
        cumulative_hazard, current_time)
    @assert(cumulative_hazard>=0)
    current_time+cumulative_hazard*scale(dist.relative_distribution)
end


function Putative(dist::TransitionDistribution, when,
        interval, consumed_interval)
    ImplicitHazardIntegral(dist, interval-consumed_interval, now)
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



"""
Weibull distribution. lambda is the scale parameter and k the
shape parameter.
F(T)=1-exp(-((T-Te)/lambda)^k)
"""
type TransitionWeibull <: TransitionDistribution
    parameters::Array{Float64,1}
end


function TransitionWeibull(lambda, k)
    TransitionWeibull([lambda, k, 0])
end


function TransitionWeibull(lambda, k, enabling_time)
    TransitionWeibull([lambda, k, enabling_time])
end


Parameters(tw::TransitionWeibull)=tw.parameters
function Parameters!(tw::TransitionWeibull, λ, k, Te)
    parameters[1]=λ
    parameters[2]=k
    parameters[3]=Te
end


function EnablingTime!(tw::TransitionWeibull, t::Float64)
    tw.parameters[3]=t
end


function Sample(dist::TransitionWeibull, now::Float64, rng::MersenneTwister)
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


function MeasuredSample(distribution::TransitionWeibull, now::Float64, rng)
    (λ, k, tₑ)=dist.parameters
    d=now-tₑ
    value=0
    mlogU=randexp(rng)
    if d>0
        value=λ*(mlogU+(d/λ)^k)^(1/k)-d
    else
        value=-d+λ*(mlogU)^(1/k)
    end
    now+value
    (value, mlogU)
end

function HazardIntegral(dist::TransitionWeibull, last, now)
    (λ, k, tₑ)=dist.parameters
    if now-tₑ>eps(Float64)
        return ((now-tₑ)/λ)^k - ((last-tₑ)/λ)^k
    else
        return 0::Float64
    end
end


function ConsumeSample(dist::TransitionWeibull, xa, start, finish)
    if xa<0
        xa=0
    end
    xa+HazardIntegral(dist, start, finish)
end


function CumulativeDistribution(dist::TransitionWeibull, when, now)
    1-exp(-hazard_integral(dist, now, when))
end


function ImplicitHazardIntegral(dist::TransitionWeibull,
        cumulative_hazard, now)
    (λ, k, tₑ)=dist.parameters
    if now-tₑ>eps(Float64)
        return tₑ + λ*(cumulative_hazard + ((now-tₑ)/λ)^k)^(1.0/k)
    else
        @debug("Weibull.implicit $cumulative_hazard")
        return tₑ + λ*(cumulative_hazard)^(1.0/k)
    end
end


function Putative(dist::TransitionWeibull, when,
        interval, consumed_interval)
    ImplicitHazardIntegral(dist, interval-consumed_interval, now)
end


function test(dist::TransitionWeibull)
    rng=MersenneTwister()
    (λ, k, tₑ)=dist.parameters
    ed=EmpiricalDistribution()
    for i in 1:10000
        push!(ed, Sample(dist, 0.0, rng))
    end
    expected_mean=λ*gamma(1+1/k)
    actual_mean=mean(ed)
    @debug("mean expected ", expected_mean, " actual ", actual_mean,
        " diff ", abs(expected_mean-actual_mean))

    expected_variance=λ^2*(gamma(1+2/k)-gamma(1+1/k)^2)
    obs_var=variance(ed)
    @debug("variance expected ", expected_variance, " actual ", obs_var,
        " diff ", abs(expected_variance-obs_var))

    min_value=min(ed)
    mink=min_value^k
    total=0.0
    for i in 1:length(ed)
        total+=ed.samples[i]^k
    end
    λ_estimator=(total/length(ed)-mink)^(1/k)
    @debug("λ expected ", λ, " actual ", λ_estimator,
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
    @debug("k expected ", k, " actual ", k_est,
        " diff ", abs(k-k_est))
end



"""
The Dirac Distribution always fires at the same time after
the enabling time. Using more than one of these in a simulation
makes it possible for two to fire at the same time, which
would violate requirements of the continuous-time model.
"""
type DiracDistribution <: TransitionDistribution
    value::Float64
    enabling_time::Float64
    DiracDistribution(value, enabling_time)=new(value, enabling_time)
end


Parameters(dd::DiracDistribution)=[dd.value, dd.enabling_time]
function Parameters!(dd::DiracDistribution, list)
    dd.value=list[1]
    dd.enabling_time=list[2]
end


EnablingTime(dd::DiracDistribution)=dd.enabling_time
function EnablingTime!(dd::DiracDistribution, when::Float64)
    dd.enabling_time=when
end


Sample(dd::DiracDistribution, now, rng)=dd.enabling_time+dd.value


CumulativeDistribution(dd::DiracDistribution)=0


function HazardIntegral(dd::DiracDistribution, t1, t2)
    absolute_time=dd.enabling_time+dd.value
    if t1 <= absolute_time <= t2
        return 1.0
    else
        return 0
    end
end


function ImplicitHazardIntegral(dd::DiracDistribution, xa, when)
    dd.enabling_time+dd.value
end



"""
UniformDistribution is between a time ta and a time tb,
both relative to the enabling time.
"""
type UniformDistribution <: TransitionDistribution
  ta::Float64
  tb::Float64
  te::Float64
  UniformDistribution(ta, tb, te)=new(ta, tb, te)
end


Parameters(ud::UniformDistribution)=[ud.ta, ud.tb, ud.te]
function Parameters!(ud::UniformDistribution, params)
    ud.ta=params[1]
    ud.tb=params[2]
    ud.te=params[3]
end


EnablingTime(ud::UniformDistribution)=ud.te
function EnablingTime!(ud::UniformDistribution, when::Float64)
    ud.te=when
end


function Sample(ud::UniformDistribution, now, rng)
    left=now<ud.ta ? ud.ta : now
    ud.te+ud.tb+(ud.tb-left)*rand(rng)
end


function HazardIntegral(ud::UniformDistribution, t1, t2)
    S0=1
    if t1-ud.te > ud.ta
        S0=1-(t1 - ud.te - ud.ta)/(ud.tb - ud.ta)
    end
    S1=1
    if t2 - ud.te > ud.ta
        S1=1-(t2 - ud.te - ud.ta)/(ud.tb - ud.ta)
    end
    if t2 > ud.tb + ud.te || t1 > ud.tb + ud.te
        return Inf
    end
    log(S0) - log(S1)
end


function ImplicitHazardIntegral(ud::UniformDistribution, xa, t1)
    if t1 - ud.te < ud.ta
        t1=ud.ta+ud.te
    end
    ud.te+ud.ta+(ud.tb-ud.ta)*(1-exp(-xa)*(1-(t1-ud.te-ud.ta)/(ud.tb-ud.ta)))
end


function CumulativeDistribution(ud::UniformDistribution, t1, now)
    # Let's be dumb and walk through every possibility given the
    # assertions.
    assert(now<=t1)
    assert(ud.ta<ud.tb)

    if now<ud.ta<ud.tb && t1<ud.ta<ud.tb
        return 0
    elseif now<ud.ta<ud.tb && ud.ta<t1<ud.tb
        return (t1-ud.ta)/(ud.tb-ud.ta)
    elseif now<ud.ta<ud.tb && ud.ta<ud.tb<t1
        return 1
    elseif ud.ta<now<ud.tb && ud.ta<t1<ud.tb
        return (t1-now)/(ud.tb-now)
    elseif ud.ta<now<ud.tb && ud.ta<ud.tb<t1
        return 1
    elseif ud.ta<ud.tb<now
        error("UniformDistribution is not defined after end time.")
    end
end


"""
Piecewise linear hazard. The hazard rate is plays connect-the-dots
between the given times.
Whatever is the last point is
treated as a horizontal line to infinity.
"""
type PiecewiseLinearDistribution <: TransitionDistribution
    b::Array{Float64, 1}
    w::Array{Float64, 1}
    te::Float64
end

function PiecewiseLinearDistribution(times, hazards, enabling_time)
    if times[length(times)]<Inf
        b=[times; Inf]
        w=[hazards; w[length[w]]]
    else
        b=times
        w=hazards
    end
    PiecewiseLinearDistribution(b, w, enabling_time)
end

Parameters(p::PiecewiseLinearDistribution)=[b, w, te]
function Parameters!(p::PiecewiseLinearDistribution, params)
    p.b=params[1]
    p.w=params[2]
    p.te=params[3]
end


function Sample(p::PiecewiseLinearDistribution, now, rng)
    error("Can't sample directly yet.")
end


"""
Integrate the hazard, taking into account when the uniform
interval starts and stops.
"""
function HazardIntegral(p::PiecewiseLinearDistribution, t0, t1)
    t0e=t0-p.te
    t1e=t1-p.te
    if t1e<p.b[1]
        return 0
    end
    if t0e<b[1]
        t0e=b[1]
    end

    total=0.0
    for idx = 1:(length(p.b)-1)
        if p.b[idx]<t1e
            db=p.b[idx+1]-p.b[idx]
            total+=0.5*db*(p.w[idx+1]+p.w[idx])
        else
            db=t1e-p.b[idx]
            total+=0.5*db*(p.w[idx+1]+p.w[idx])
            break
        end
    end
    total
end


function ImplicitHazardIntegral(p::PiecewiseLinearDistribution, xa, t0)
    t0e=t0-self.te
    if t0e<b[1]
        t0e=b[1]
    end
    error("Just not finished")
    for i = 1:(length(p.b)-1)
        if t0e>p.b[i+1]
            continue
        elseif t0e>=p.b[i+1]
            continue
        end
    end
    t0e
end


"""
A Gamma distribution.
α - shape parameter
β - inverse scale parameter, also called rate parameter

pdf=(β^α/Γ(α))x^(α-1) e^(-βx)
"""
function TransitionGamma(α::Float64, β::Float64, te::Float64)
    # The supplied version uses θ=1/β.
    relative=Distributions.Gamma(α, 1/β)
    WrappedDistribution(relative, te)
end

##################################
# F(t)=1/(1 + ((t-te)/α)^(-β))
type LogLogistic <: Distributions.ContinuousUnivariateDistribution
    alpha::Float64
    beta::Float64
end


function rand(d::LogLogistic, rng::MersenneTwister)
    quantile(d, rand(rng))
end

function quantile(d::LogLogistic, U::Float64)
    d.alpha*(U/(1-U))^(1/d.beta)
end

function cdf(d::LogLogistic, t::Real)
    1/(1+(t/d.alpha)^(-d.beta))
end

# Survival
function ccdf(d::LogLogistic, t::Real)
    1/(1+(t/d.alpha)^d.beta)
end

function logccdf(d::LogLogistic, t::Real)
    -log( 1 + (t/d.alpha)^d.beta )
end

function invlogccdf(d::LogLogistic, lp::Float64)
    d.alpha*(1-exp(-lp)^(1/d.beta))
end

function TransitionLogLogistic(a::Float64, b::Float64, t::Float64)
    WrappedDistribution(LogLogistic(a,b), t)
end


"""
An empirical distribution is an estimator of a distribution
given a set of samples of times at which it fired.
For this implementation, first make the `EmpiricalDistribution`.
Then use `push!` to add values. Then call `build!` before
sampling from it.
"""
type EmpiricalDistribution
    samples::Array{Float64,1}
    built::Bool
    EmpiricalDistribution()=new(Array(Float64,0), false)
end

function cdf(ed::EmpiricalDistribution, which::Int)
    (ed.samples[which], which/length(ed.samples))
end

function build!(ed::EmpiricalDistribution)
    if !ed.built
        sort!(ed.samples)
        ed.built=true
    end
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

"""
This compares two distributions. It returns the
largest difference between an `EmpiricalDistribution`
and the `other` distribution, along with a `Bool`
that says whether the null hypothesis, that they agree, is rejected
at level 0.05.
"""
function kolmogorov_smirnov_statistic(ed::EmpiricalDistribution, other)
    build!(ed)
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


"""
Given a set of observations of firing times and
cancellation times of an event,
this constructs a Nelson-Aalen estimator for the
distribution.
"""
type NelsonAalenEntry
    when::Float64
    hazard_sum::Float64
end

"""
H(t)=\int \lambda=\sum_{t_i<=t} (d/n)
"""
type NelsonAalenDistribution
    integrated_hazard::Array{NelsonAalenEntry,1}
    NelsonAalenDistribution(cnt::Int)=new(Array(NelsonAalenEntry, cnt))
end

function cdf(dist::NelsonAalenDistribution, when::Float64)
    entry_idx=findfirst(x->(x.when>when), dist.integrated_hazard)
    if entry_idx==0
        entry_idx=length(dist.hazard_sum)+1
    end
    1-exp(-dist.integrated_hazard[entry_idx-1].hazard_sum)
end


function cdf(dist::NelsonAalenDistribution, bypoint::Int)
    entry=dist.integrated_hazard[bypoint]
    (entry.when, 1-exp(-entry.hazard_sum))
end


function kolmogorov_smirnov_statistic(ed::NelsonAalenDistribution, other)
    sup_diff=0.0
    n=length(ed.integrated_hazard)
    for i in 1:n
        t=ed.integrated_hazard[i].when
        F_empirical=ed.integrated_hazard[i].hazard_sum
        sup_diff=max(sup_diff, abs(F_empirical-cdf(other, t)))
    end
    c_alpha=1.36 # 0.05 confidence interval
    (sup_diff, sup_diff > c_alpha*sqrt(2/n))
end


type NelsonAalenSortEntry
    time::Float64
    fired::Int
end

function isless(a::NelsonAalenSortEntry, b::NelsonAalenSortEntry)
    a.time<b.time
end

function NelsonAalenDistribution(fired::Array{Float64,1}, not_fired::Array{Float64,1})
end


"""
Given measurements of several outcomes, construct Nelson-Aalen distributions
"""
function multiple_measures(eds::Array{EmpiricalDistribution,1})
    dist_len=Int64[length(x) for x in eds]
    nad=Array(NelsonAalenDistribution, length(eds))
    for nad_idx=1:length(eds)
        # +1 for the entry with integrated hazard of zero
        nad[nad_idx]=NelsonAalenDistribution(dist_len[nad_idx]+1)
    end
    for ed in eds
        build!(ed)
    end
    total=sum(dist_len)
    fired=Array(NelsonAalenSortEntry, total)
    fired_idx=1
    for add_idx=1:length(eds)
        for s in eds[add_idx].samples
            fired[fired_idx]=NelsonAalenSortEntry(s, add_idx)
            fired_idx+=1
        end
    end
    sort!(fired)
    dist_idx=ones(Int, length(eds))
    for start_zero=1:length(eds)
        nad[start_zero].integrated_hazard[dist_idx[start_zero]]=NelsonAalenEntry(0.0, 0.0)
        dist_idx[start_zero]+=1
    end
    for entry_idx=1:length(fired)
        at_risk=total-entry_idx+1
        who=fired[entry_idx].fired
        when=fired[entry_idx].time
        previous=nad[who].integrated_hazard[dist_idx[who]-1].hazard_sum
        nad[who].integrated_hazard[dist_idx[who]]=NelsonAalenEntry(when, previous+1.0/at_risk)
        dist_idx[who]+=1
    end
    for check_idx=1:length(eds)
        assert(dist_idx[check_idx]==2+dist_len[check_idx])
    end
    nad
end

