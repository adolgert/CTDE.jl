
"""
Weibull distribution. lambda is the scale parameter and k the
shape parameter.
F(T)=1-exp(-((T-Te)/lambda)^k)
"""
struct TransitionWeibull <: TransitionDistribution
    parameters::Array{Float64,1}
    te::Float64
    TransitionWeibull(lambda, k)=new([lambda, k], 0.0)
end


Parameters(tw::TransitionWeibull)=tw.parameters
function Parameters!(tw::TransitionWeibull, λ, k)
    tw.parameters[1]=λ
    tw.parameters[2]=k
end


function EnablingTime!(tw::TransitionWeibull, t::Float64)
    tw.te=t
end


function Sample(dist::TransitionWeibull, now::Float64, rng::MersenneTwister)
    (λ, k)=dist.parameters
    d=now-dist.te
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
    (λ, k)=distribution.parameters
    d=now-distribution.te
    value=0
    mlogU=randexp(rng)
    if d>0
        value=λ*(mlogU+(d/λ)^k)^(1/k)-d
    else
        value=-d+λ*(mlogU)^(1/k)
    end
    (now+value, mlogU)
end

function HazardIntegral(dist::TransitionWeibull, last, now)
    (λ, k)=dist.parameters
    if now-dist.te>eps(Float64)
        return ((now-dist.te)/λ)^k - ((last-dist.te)/λ)^k
    else
        return 0::Float64
    end
end


function ConsumeSample(dist::TransitionWeibull, xa, start, finish)
    xa=(xa<0) ? 0 : xa
    xa+HazardIntegral(dist, start, finish)
end


function CumulativeDistribution(dist::TransitionWeibull, when, now)
    1-exp(-hazard_integral(dist, now, when))
end


function ImplicitHazardIntegral(dist::TransitionWeibull,
        cumulative_hazard, when)
    (λ, k)=dist.parameters
    if when-dist.te>eps(Float64)
        return dist.te + λ*(cumulative_hazard + ((when-dist.te)/λ)^k)^(1.0/k)
    else
        return dist.te + λ*(cumulative_hazard)^(1.0/k)
    end
end


function Putative(dist::TransitionWeibull, when,
        interval, consumed_interval)
    ImplicitHazardIntegral(dist, interval-consumed_interval, when)
end


"""
Given a distribution and an EmpiricalDistribution with
samples from it, compute a set of statistics on those samples,
returning a dictionary where the keys are the names of
the statistics and the entries are an array of
[expected, estimated].
"""
function TestStatistics(d::TransitionWeibull, ed::EmpiricalDistribution)
    statistic=Dict{AbstractString, Array}()
    (λ, k)=d.parameters
    expected_mean=λ*gamma(1+1/k)
    actual_mean=mean(ed)
    statistic["mean"]=[expected_mean, actual_mean]
    statistic["variance"]=[
        (λ^2)*(gamma(1+2/k)-gamma(1+1/k)^2),
        variance(ed)]

    min_value=min(ed)
    mink=min_value^k
    total=0.0
    for i in 1:length(ed)
        total+=ed.samples[i]^k
    end
    statistic["scale"]=[λ,
        (total/length(ed)-mink)^(1/k)]


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

    statistic["exponent"]=[k, k_est]
    statistic
end


function test(dist::TransitionWeibull)
    rng=MersenneTwister()
    (λ, k)=dist.parameters
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

